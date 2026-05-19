## Test: cond_threshold=1e4 with bootstrap SE comparison (fast version)
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Mental_Health_Data_Application")
library(Matrix)
library(numDeriv)
library(parallel)
library(foreach)
library(doSNOW)
source("MHD_functions.R")

original_data <- read.csv("data_application.csv")
percent <- original_data$Percentage
n_total <- 2486
class_n <- round(n_total * percent / 100)
dat <- gen_data(original_data, class_n, n_total)
dat$y <- dat$teacher_report
n <- nrow(dat)

W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
compute_pi_fn <- function(eta) plogis(-eta)
COND_THRESHOLD <- 1e5

## Build design matrix and h_x manually (no package needed)
## design_mat columns determined by R formula; h_x from named columns + intercept
build_matrices_manual <- function(dat_b, formula_idx) {
  # All 3 formulas share same h_alpha columns
  h_alpha_names <- c("health", "father", "parent_report", "fp", "fh", "hp")

  # Build design matrix via model.matrix
  formulas <- list(
    ~ teacher_report + health + father + teacher_report:health + teacher_report:father,
    ~ teacher_report + health + parent_report + teacher_report:health + teacher_report:parent_report,
    ~ teacher_report + parent_report + father + teacher_report:parent_report + teacher_report:father
  )
  design_mat <- model.matrix(formulas[[formula_idx]], data = dat_b)

  # Build h_x: intercept + discrete vars as factors + continuous vars
  # All h_alpha vars are binary -> treated as factors -> model.matrix adds intercept
  h_x1 <- dat_b[h_alpha_names]
  for (j in 1:ncol(h_x1)) h_x1[, j] <- as.factor(h_x1[, j])
  d <- model.matrix(lm(rep(1, nrow(dat_b)) ~ ., data = h_x1))
  h_x <- d  # no continuous h_x2 vars

  list(design_mat = design_mat, h_x = h_x)
}

## Inline GMM with cond_threshold (simplified: fewer restarts, no eigen at every step)
run_gmm_cond <- function(design_mat, h_x, r_vec, n, cond_threshold = COND_THRESHOLD) {
  k <- ncol(design_mat); h_dim <- ncol(h_x)
  init <- rep(0, k); tol <- 1e-6; max_outer <- 200L; trust_radius <- 2.0

  g_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    (r_vec / pi_v - 1) * h_x
  }
  Gamma_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    cf <- r_vec * (1 - pi_v) / pi_v
    crossprod(h_x * cf, design_mat) / n
  }
  check_cond <- function(alpha, W_hat) {
    Gamma <- Gamma_fn(alpha)
    H <- crossprod(Gamma, W_hat %*% Gamma)
    eig <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    max(eig) / max(min(eig), 1e-15)
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
          # Cond check only on accepted Armijo steps
          cc <- check_cond(cand, W_hat)
          if (cc < cond_threshold) { accepted <- TRUE; break }
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

  pi_final <- compute_pi_fn(as.vector(design_mat %*% estimates))
  list(alpha = estimates, pi = pi_final)
}

## SE computation (CUE SE2)
compute_se2 <- function(design_mat, h_x, r_vec, y_vec, alpha_hat, n) {
  pi_hat <- compute_pi_fn(as.vector(design_mat %*% alpha_hat))
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

#------------------------------------------------------------------------------#
# Point estimates
#------------------------------------------------------------------------------#
cat("=== Point estimates with cond_threshold = 1e4 ===\n\n")

mu_cc <- mean(dat$teacher_report[dat$r == 1])
se_cc <- sd(dat$teacher_report[dat$r == 1]) / sqrt(sum(dat$r))
cat(sprintf("Complete case: mu=%.4f, se=%.4f\n\n", mu_cc, se_cc))

point_results <- list()
for (i in 1:3) {
  di <- build_matrices_manual(dat, i)
  set.seed(42)
  gmm_res <- run_gmm_cond(di$design_mat, di$h_x, dat$r, n)
  se_res <- compute_se2(di$design_mat, di$h_x, dat$r, dat$teacher_report, gmm_res$alpha, n)
  point_results[[i]] <- list(mu = se_res$mu, se = se_res$se)
  ps <- gmm_res$pi
  cat(sprintf("Model %d: mu=%.4f, se=%.4f, PS=[%.4f,%.4f], >0.99:%d, <0.01:%d\n",
      i, se_res$mu, se_res$se, min(ps), max(ps), sum(ps>0.99), sum(ps<0.01)))
}

#------------------------------------------------------------------------------#
# Bootstrap (parallel, no package loading)
#------------------------------------------------------------------------------#
cat("\n=== Bootstrap (B=1000) ===\n\n")

B <- 1000
n_cores <- min(detectCores() - 1, 10)
cat(sprintf("cores=%d\n", n_cores))

cl <- makeCluster(n_cores)
registerDoSNOW(cl)
clusterExport(cl, c("dat", "n", "W_fn", "compute_pi_fn", "COND_THRESHOLD",
                     "run_gmm_cond", "build_matrices_manual"),
              envir = environment())

pb <- txtProgressBar(max = B, style = 3)
progress <- function(nn) setTxtProgressBar(pb, nn)
opts <- list(progress = progress)

boot_results <- foreach(
  b = 1:B, .combine = 'rbind', .options.snow = opts
) %dopar% {
  set.seed(12345 + b)
  idx <- sample(1:n, n, replace = TRUE)
  dat_b <- dat[idx, ]
  n_b <- nrow(dat_b)
  r_b <- dat_b$r

  mu_cc_b <- mean(dat_b$teacher_report[r_b == 1])
  mu_vec <- rep(NA, 3)
  for (i in 1:3) {
    tryCatch({
      di <- build_matrices_manual(dat_b, i)
      gmm_b <- run_gmm_cond(di$design_mat, di$h_x, r_b, n_b)
      mu_vec[i] <- mean(r_b * dat_b$teacher_report / gmm_b$pi)
    }, error = function(e) {})
  }
  c(CC = mu_cc_b, M1 = mu_vec[1], M2 = mu_vec[2], M3 = mu_vec[3])
}

close(pb)
stopCluster(cl)

#------------------------------------------------------------------------------#
# Results
#------------------------------------------------------------------------------#
cat("\n\n=== Results ===\n\n")

compute_trimmed_se <- function(x) {
  x <- x[!is.na(x)]
  n_valid <- length(x)
  Q1 <- quantile(x, 0.25); Q3 <- quantile(x, 0.75); IQR_val <- Q3 - Q1
  outlier <- x < (Q1 - 2 * IQR_val) | x > (Q3 + 2 * IQR_val)
  list(se_raw = sd(x), se_trim = sd(x[!outlier]), n_out = sum(outlier), n_valid = n_valid)
}

labels <- c("CC", "M1", "M2", "M3")
estimates <- c(mu_cc, point_results[[1]]$mu, point_results[[2]]$mu, point_results[[3]]$mu)
analytical_se <- c(se_cc, point_results[[1]]$se, point_results[[2]]$se, point_results[[3]]$se)

cat(sprintf("  %-5s %10s %12s %12s %12s %10s %10s\n",
            "Label", "Estimate", "Analytical", "Boot(raw)", "Boot(trim)", "Outliers", "Valid"))
cat("  ", paste(rep("-", 76), collapse = ""), "\n")

for (j in 1:4) {
  bt <- compute_trimmed_se(boot_results[, j])
  cat(sprintf("  %-5s %10.4f %12.4f %12.4f %12.4f %10d %10d\n",
              labels[j], estimates[j], analytical_se[j], bt$se_raw, bt$se_trim, bt$n_out, bt$n_valid))
}
cat("  ", paste(rep("-", 76), collapse = ""), "\n")
