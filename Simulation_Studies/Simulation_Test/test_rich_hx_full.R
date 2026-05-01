## Full simulation: Enriched h_x = (1, u1, u2, z1, z2, u2^2, z2^2, u1*u2)
## 8 instruments, 4 params, overid=4
## Test S4 M3, S4 M2, S3 M2
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)
library(parallel)

n_val <- 2000; n_reps <- 1000
ps_spec <- get_ps_spec("9-alt1")
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
compute_pi_fn <- function(eta) plogis(-eta)
cond_limit <- 1e8

run_gmm <- function(dat, design_mat, h_x, link_type) {
  r_vec <- dat[["r"]]
  n <- nrow(dat); k <- ncol(design_mat); h_dim <- ncol(h_x)
  init <- rep(0, k)
  tol <- 1e-6; max_outer <- 1000L; trust_radius <- 0.5

  g_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    (r_vec / pi_v - 1) * h_x
  }
  Gamma_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    crossprod(h_x * (r_vec * (1 - pi_v) / pi_v), design_mat) / n
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
      for (ls in 1:30) {
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
  for (t in 1:max_outer) {
    g_mat <- g_fn(estimates)
    G_vec <- matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
    W_hat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(h_dim))
    Gamma <- Gamma_fn(estimates)
    cg <- 2 * as.vector(crossprod(Gamma, W_hat %*% G_vec))
    if (max(abs(cg)) < tol) break
    H_gn <- crossprod(Gamma, W_hat %*% Gamma)
    direction <- tryCatch(solve(H_gn, -cg), error = function(e) -cg)
    sn <- sqrt(sum(direction^2))
    if (sn > trust_radius) direction <- direction * (trust_radius / sn)
    obj_cur <- as.numeric(crossprod(G_vec, W_hat %*% G_vec))
    step <- 1.0; accepted <- FALSE
    for (ls in 1:30) {
      cand <- estimates + step * direction
      g_cand <- g_fn(cand)
      G_cand <- matrix(.colMeans(g_cand, n, h_dim), h_dim, 1)
      obj_new <- as.numeric(crossprod(G_cand, W_hat %*% G_cand))
      if (is.finite(obj_new) && obj_new < obj_cur - 1e-4 * step * sum(cg * direction)) {
        accepted <- TRUE; break
      }
      step <- step * 0.5
    }
    if (accepted) estimates <- cand
  }
  g_mat_f <- g_fn(estimates)
  pi_final <- compute_pi_fn(as.vector(design_mat %*% estimates))
  cond_val <- tryCatch({
    W_f <- solve(crossprod(g_mat_f) / n)
    Gamma_f <- Gamma_fn(estimates)
    H <- crossprod(Gamma_f, W_f %*% Gamma_f)
    ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    max(ev) / max(min(ev), 1e-15)
  }, error = function(e) Inf)
  list(estimates = estimates, pi_hat = pi_final, cond = cond_val)
}

compute_se <- function(dat, alpha_hat, design_mat, h_x, link_type, n) {
  r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
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
  if (!is.null(link_type) && link_type == "logistic_complement") {
    dot_pi <- -design_mat * (pi_hat * (1 - pi_hat))
  } else {
    dot_pi <- design_mat * (pi_hat * (1 - pi_hat))
  }
  ry_ps_inv2 <- as.vector(r_vec * y_vec * (pi_hat^(-2)))
  H_alpha <- colMeans(dot_pi * ry_ps_inv2)
  mu_iid <- as.vector(t(r_vec/pi_hat*y_vec) - t(H_alpha) %*% psi)
  se <- sqrt(var(mu_iid) / n)
  list(mu = mu_hat, se = se)
}

## Build enriched h_x for a dataset
build_h_rich <- function(dat) {
  h <- cbind(1, dat$u1, dat$u2, dat$z1, dat$z2,
             dat$u2^2, dat$z2^2, dat$u1 * dat$u2)
  colnames(h) <- c("1","u1","u2","z1","z2","u2sq","z2sq","u1u2")
  h
}

settings <- list(
  list(name = "S4 M3", setting = "setting4", m_idx = 3),
  list(name = "S4 M2", setting = "setting4", m_idx = 2),
  list(name = "S3 M2", setting = "setting3", m_idx = 2)
)

cat("=== Full simulation: Enriched h_x (8 instruments, overid=4) ===\n\n")

for (cfg in settings) {
  cat(sprintf("========== %s ==========\n", cfg$name))
  data_file <- misspecified_model_all_data_file.list[[cfg$setting]][["miss50"]][[1]]
  all_data <- readRDS(data_file)
  mu_true <- get_mu_true(cfg$setting)
  cat(sprintf("True mean: %.6f\n", mu_true))

  cat("Pre-building design matrices...\n")
  design_info <- vector("list", n_reps)
  for (rep_i in 1:n_reps) {
    dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
    single_ps <- list(
      formula.list = list(ps_spec[["formula.list"]][[cfg$m_idx]]),
      h_alpha.list = list(ps_spec[["h_alpha.list"]][[cfg$m_idx]]),
      inv_link = ps_spec[["inv_link"]], outcome = ps_spec[["outcome"]],
      alpha_init.list = list(NULL), optimizer = "L-BFGS-B"
    )
    tryCatch({
      ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
      dm <- ebmr$ps_fit.list[[1]]$design_matrix
      lt <- ebmr$ps_fit.list[[1]]$link_type
      h_rich <- build_h_rich(dat)
      design_info[[rep_i]] <- list(design_mat = dm, h_x = h_rich, link_type = lt)
    }, error = function(e) NULL)
  }
  cat("Done.\n")

  n_cores <- min(detectCores() - 1, 10)
  cl <- makeCluster(n_cores)
  clusterExport(cl, c("run_gmm", "compute_se", "compute_pi_fn", "W_fn",
                       "all_data", "design_info", "n_val"),
                envir = environment())

  results <- parLapply(cl, 1:n_reps, function(rep_i) {
    if (is.null(design_info[[rep_i]])) return(list(mu = NA, se = NA, cond = NA))
    dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
    di <- design_info[[rep_i]]
    tryCatch({
      res <- run_gmm(dat, di$design_mat, di$h_x, di$link_type)
      se_res <- compute_se(dat, res$estimates, di$design_mat, di$h_x, di$link_type, n_val)
      list(mu = se_res$mu, se = se_res$se, cond = res$cond)
    }, error = function(e) list(mu = NA, se = NA, cond = NA))
  })
  stopCluster(cl)

  mu_vals <- sapply(results, `[[`, "mu")
  se_vals <- sapply(results, `[[`, "se")
  cond_vals <- sapply(results, `[[`, "cond")

  valid <- !is.na(mu_vals) & !is.na(se_vals)
  n_degen <- sum(cond_vals[valid] >= cond_limit, na.rm = TRUE)

  cat(sprintf("  Valid: %d/%d, Degen: %d\n", sum(valid), n_reps, n_degen))
  cat(sprintf("  Bias     = %.4f\n", mean(mu_vals[valid]) - mu_true))
  cat(sprintf("  ESD      = %.4f\n", sd(mu_vals[valid])))
  cat(sprintf("  ESE      = %.4f\n", mean(se_vals[valid])))
  cat(sprintf("  ESE/ESD  = %.3f\n", mean(se_vals[valid]) / sd(mu_vals[valid])))
  ci_lo <- mu_vals[valid] - 1.96 * se_vals[valid]
  ci_hi <- mu_vals[valid] + 1.96 * se_vals[valid]
  cat(sprintf("  CP       = %.3f\n", mean(mu_true >= ci_lo & mu_true <= ci_hi)))
  cat("\n")
}
