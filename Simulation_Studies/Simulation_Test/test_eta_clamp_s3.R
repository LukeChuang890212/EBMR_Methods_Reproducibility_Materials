## Test: eta_max clamping to address NA reps in S3, scenario 9-2, model 2, n=500, miss30
## Clamp eta = design_mat %*% alpha to [-eta_max, eta_max] before computing pi
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)
library(parallel)

ETA_MAX <- 20  # plogis(-20) ~ 2e-9, plogis(20) ~ 1-2e-9

n_val <- 500; n_reps <- 1000
ps_spec <- get_ps_spec("9")
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
compute_pi_clamped <- function(eta) plogis(-pmin(pmax(eta, -ETA_MAX), ETA_MAX))
compute_pi_unclamped <- function(eta) plogis(-eta)

mu_true <- get_mu_true("setting3")

## GMM with eta clamping
run_gmm_clamped <- function(design_mat, h_x, r_vec, n) {
  k <- ncol(design_mat); h_dim <- ncol(h_x)
  init <- rep(0, k); tol <- 1e-6; max_outer <- 200L; trust_radius <- 2.0

  g_fn <- function(alpha) {
    eta <- as.vector(design_mat %*% alpha)
    pi_v <- compute_pi_clamped(eta)
    (r_vec / pi_v - 1) * h_x
  }
  Gamma_fn <- function(alpha) {
    eta <- as.vector(design_mat %*% alpha)
    pi_v <- compute_pi_clamped(eta)
    cf <- r_vec * (1 - pi_v) / pi_v
    crossprod(h_x * cf, design_mat) / n
  }

  nr_inner <- function(start, W_hat) {
    alpha <- start
    for (iter in 1:500) {
      g_mat <- g_fn(alpha)
      G_vec <- matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
      Gamma <- Gamma_fn(alpha)
      gv <- 2 * as.vector(crossprod(Gamma, W_hat %*% G_vec))
      if (max(abs(gv)) < 1e-8) break
      H_gn <- crossprod(Gamma, W_hat %*% Gamma)
      direction <- tryCatch(solve(H_gn, -gv), error = function(e) -gv)
      sn <- sqrt(sum(direction^2))
      if (sn > trust_radius) direction <- direction * (trust_radius / sn)
      obj_cur <- as.numeric(crossprod(G_vec, W_hat %*% G_vec))
      step <- 1.0; accepted <- FALSE
      for (ls in 1:20) {
        cand <- alpha + step * direction
        g_cand <- g_fn(cand)
        G_cand <- matrix(.colMeans(g_cand, n, h_dim), h_dim, 1)
        obj_new <- as.numeric(crossprod(G_cand, W_hat %*% G_cand))
        if (is.finite(obj_new) && obj_new < obj_cur - 1e-4 * step * sum(gv * direction)) {
          accepted <- TRUE; break
        }
        step <- step * 0.5
      }
      if (accepted) alpha <- cand else break
    }
    alpha
  }

  estimates <- nr_inner(init, diag(h_dim))
  best_grad_norm <- Inf; best_estimates <- estimates; no_improve <- 0L

  for (t in 1:max_outer) {
    g_mat <- g_fn(estimates)
    G_vec <- matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
    W_hat <- tryCatch(W_fn(g_mat), error = function(e) diag(h_dim))
    cg <- 2 * as.vector(crossprod(Gamma_fn(estimates), W_hat %*% G_vec))
    grad_norm <- max(abs(cg))
    if (grad_norm < tol) break
    if (grad_norm < best_grad_norm) {
      best_grad_norm <- grad_norm; best_estimates <- estimates; no_improve <- 0L
    } else {
      no_improve <- no_improve + 1L
      if (no_improve >= 20L) { estimates <- best_estimates; break }
    }
    estimates <- nr_inner(estimates, W_hat)
  }

  pi_final <- compute_pi_clamped(as.vector(design_mat %*% estimates))
  list(alpha = estimates, pi = pi_final)
}

## SE computation with clamped pi
compute_se_clamped <- function(design_mat, h_x, r_vec, y_vec, alpha_hat, n) {
  eta <- as.vector(design_mat %*% alpha_hat)
  pi_hat <- compute_pi_clamped(eta)
  mu_hat <- mean(r_vec * y_vec / pi_hat)
  g_mat <- (r_vec / pi_hat - 1) * h_x
  W_hat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) NULL)
  if (is.null(W_hat)) return(list(mu = mu_hat, se = NA))
  cf <- r_vec * (1 - pi_hat) / pi_hat
  Gamma_hat <- crossprod(h_x * cf, design_mat) / n
  GtWG <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
  k <- length(alpha_hat); h <- ncol(h_x)
  GtW <- crossprod(Gamma_hat, W_hat)
  G_vec <- colMeans(g_mat)
  dc_X <- cf * design_mat
  R_mat <- matrix(0, h * k, k)
  for (j in 1:k) {
    block <- crossprod(dc_X, design_mat[, j] * h_x) / n
    for (l in 1:h) R_mat[(l-1)*k + j, ] <- block[, l]
  }
  eta_s <- matrix(G_vec, h, 1)
  WG <- as.vector(W_hat %*% eta_s)
  GW_row <- as.vector(t(eta_s) %*% W_hat)
  S_mat <- matrix(0, h^2, k)
  for (j in 1:k) {
    dg_j <- cf * design_mat[, j] * h_x
    cross <- (crossprod(dg_j, g_mat) + crossprod(g_mat, dg_j)) / n
    S_mat[, j] <- as.vector(cross)
  }
  GW_kron_Ik <- kronecker(matrix(GW_row, 1, h), diag(k))
  GW_kron_GtW <- kronecker(matrix(GW_row, 1, h), GtW)
  H_mat <- GtWG + GW_kron_Ik %*% R_mat - GW_kron_GtW %*% S_mat
  H_cond <- tryCatch({
    ev <- eigen(H_mat, symmetric = FALSE, only.values = TRUE)$values
    max(Mod(ev)) / max(min(Mod(ev)), 1e-15)
  }, error = function(e) Inf)
  if (!is.finite(H_cond) || H_cond > 1e15) return(list(mu = mu_hat, se = NA))
  H_inv <- tryCatch(solve(H_mat), error = function(e) NULL)
  if (is.null(H_inv)) return(list(mu = mu_hat, se = NA))
  GtW_g <- GtW %*% t(g_mat)
  g_WG <- as.vector(g_mat %*% WG)
  cf_hxWG <- cf * as.vector(h_x %*% WG)
  GammaT_WG <- t(design_mat * cf_hxWG)
  Q_mat <- GtW_g + GammaT_WG - sweep(GtW_g, 2, g_WG, `*`)
  psi <- -(H_inv %*% Q_mat)
  dot_pi <- -design_mat * (pi_hat * (1 - pi_hat))
  ry_ps_inv2 <- as.vector(r_vec * y_vec * (pi_hat^(-2)))
  H_alpha <- colMeans(dot_pi * ry_ps_inv2)
  mu_iid <- as.vector(t(r_vec/pi_hat*y_vec) - t(H_alpha) %*% psi)
  se <- sqrt(var(mu_iid) / n)
  list(mu = mu_hat, se = se)
}

## Build design matrices manually without calling package constructor
## Model 2 of scenario 9: u1_z2 = r ~ y + u1 + z2
## design_mat = (Intercept, y, u1, z2) with y-terms first after intercept
## h_alpha = full = c("u1", "u2", "z1", "z2") -> h_x = (1, u1, u2, z1, z2)
## (u1 and z1 are binary -> factor expansion adds nothing since 2 levels -> 1 col each)
build_dm <- function(dat_b) {
  # Design matrix: (Intercept, y, u1, z2)
  design_mat <- cbind(1, dat_b$y, dat_b$u1, dat_b$z2)
  colnames(design_mat) <- c("(Intercept)", "y", "u1", "z2")
  # h_x: (1, u1, u2, z1, z2, u1u2, z1z2, u1z1, z1u2, u1z2, u2z2)
  h_x <- cbind(1, dat_b$u1, dat_b$u2, dat_b$z1, dat_b$z2,
               dat_b$u1*dat_b$u2, dat_b$z1*dat_b$z2,
               dat_b$u1*dat_b$z1, dat_b$z1*dat_b$u2,
               dat_b$u1*dat_b$z2, dat_b$u2*dat_b$z2)
  colnames(h_x) <- c("1","u1","u2","z1","z2","u1u2","z1z2","u1z1","z1u2","u1z2","u2z2")
  list(design_mat = design_mat, h_x = h_x)
}

## Load data
all_data <- readRDS("Simulation_Data/Setting3.B2_n500_replicate1000.RDS")

cat(sprintf("=== Testing eta_max = %d clamping on S3, 9-2, M2, n=500, miss30 ===\n\n", ETA_MAX))
cat(sprintf("True mean: %.4f\n\n", mu_true))

## Run all 1000 reps with clamping (parallel)
n_cores <- min(detectCores() - 1, 10)
cl <- makeCluster(n_cores)
clusterExport(cl, c("all_data", "n_val", "compute_pi_clamped", "W_fn",
                     "run_gmm_clamped", "compute_se_clamped", "build_dm",
                     "ETA_MAX", "compute_pi_unclamped"), envir = environment())

results <- parLapply(cl, 1:n_reps, function(rep_i) {
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
  tryCatch({
    di <- build_dm(dat)
    gmm_res <- run_gmm_clamped(di$design_mat, di$h_x, dat$r, n_val)
    se_res <- compute_se_clamped(di$design_mat, di$h_x, dat$r, dat$y, gmm_res$alpha, n_val)
    list(mu = se_res$mu, se = se_res$se,
         pi_min = min(gmm_res$pi), pi_max = max(gmm_res$pi),
         alpha = gmm_res$alpha)
  }, error = function(e) {
    list(mu = NA, se = NA, pi_min = NA, pi_max = NA, alpha = NA)
  })
})
stopCluster(cl)

mu_vals <- sapply(results, `[[`, "mu")
se_vals <- sapply(results, `[[`, "se")

valid <- !is.na(mu_vals) & !is.na(se_vals)
n_na <- sum(!valid)
mu_v <- mu_vals[valid]; se_v <- se_vals[valid]

# Outlier detection
Q1 <- quantile(mu_v, 0.25); Q3 <- quantile(mu_v, 0.75); IQR_val <- Q3 - Q1
is_outlier <- mu_v < (Q1 - 3 * IQR_val) | mu_v > (Q3 + 3 * IQR_val)
max_remove <- floor(0.01 * length(mu_v))
if (sum(is_outlier) > max_remove && max_remove > 0) {
  dist_med <- abs(mu_v - median(mu_v))
  oi <- which(is_outlier)
  keep <- oi[order(dist_med[oi], decreasing = TRUE)[(max_remove+1):length(oi)]]
  is_outlier[keep] <- FALSE
}
n_out <- sum(is_outlier)
mu_c <- mu_v[!is_outlier]; se_c <- se_v[!is_outlier]

cat(sprintf("--- With eta_max = %d clamping ---\n", ETA_MAX))
cat(sprintf("  Total: %d, NA: %d, Valid: %d, Outliers: %d\n", n_reps, n_na, sum(valid), n_out))
cat(sprintf("\n  Raw (n=%d):\n", sum(valid)))
cat(sprintf("    Bias     = %.4f\n", mean(mu_v) - mu_true))
cat(sprintf("    ESD      = %.4f\n", sd(mu_v)))
cat(sprintf("    ESE      = %.4f\n", mean(se_v)))
cat(sprintf("    ESE/ESD  = %.3f\n", mean(se_v) / sd(mu_v)))
ci_lo <- mu_v - 1.96 * se_v; ci_hi <- mu_v + 1.96 * se_v
cat(sprintf("    CP       = %.3f\n", mean(mu_true >= ci_lo & mu_true <= ci_hi)))

cat(sprintf("\n  Cleaned (n=%d):\n", length(mu_c)))
cat(sprintf("    Bias     = %.4f\n", mean(mu_c) - mu_true))
cat(sprintf("    ESD      = %.4f\n", sd(mu_c)))
cat(sprintf("    ESE      = %.4f\n", mean(se_c)))
cat(sprintf("    ESE/ESD  = %.3f\n", mean(se_c) / sd(mu_c)))
ci_lo <- mu_c - 1.96 * se_c; ci_hi <- mu_c + 1.96 * se_c
cat(sprintf("    CP       = %.3f\n", mean(mu_true >= ci_lo & mu_true <= ci_hi)))

## Compare with original (no clamping) from saved file
cat("\n--- Original (no clamping, from saved file) ---\n")
sim_orig <- readRDS("Simulation_Results/EBMR_IPW_setting3-miss30-scenario9-2_2_n500_replicate1000_test59.RDS")
cleaned_orig <- clean_sim_result(sim_orig, multiplier = 3, verbose = FALSE)
sim_c <- cleaned_orig$result
cat(sprintf("  Total: %d, NA: %d, Outliers: %d, Used: %d\n",
    cleaned_orig$n_total, cleaned_orig$n_na, cleaned_orig$n_outliers, cleaned_orig$n_successful))
cat(sprintf("  Bias     = %.4f\n", mean(sim_c[1, ]) - mu_true))
cat(sprintf("  ESD      = %.4f\n", sd(sim_c[1, ])))
cat(sprintf("  ESE      = %.4f\n", mean(sim_c[3, ])))
cat(sprintf("  ESE/ESD  = %.3f\n", mean(sim_c[3, ]) / sd(sim_c[1, ])))
ci_lo <- sim_c[1, ] - 1.96 * sim_c[3, ]; ci_hi <- sim_c[1, ] + 1.96 * sim_c[3, ]
cat(sprintf("  CP       = %.3f\n", mean(mu_true >= ci_lo & mu_true <= ci_hi)))
