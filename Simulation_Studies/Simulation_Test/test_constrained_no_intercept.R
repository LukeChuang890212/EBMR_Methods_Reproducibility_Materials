## Constrained GMM: minimize Q(alpha) s.t. E[r/pi(alpha)-1]=0
## h_alpha drops the intercept (since E[r/pi-1]=0 is enforced as constraint)
## Uses nloptr::slsqp for constrained optimization
## Iterative GMM: fix W, solve constrained problem, update W, repeat

setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)
library(nloptr)

n_val <- 2000
n_reps <- 1000
ps_spec <- get_ps_spec("9-alt1")
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
compute_pi_fn <- function(eta) plogis(-eta)
tol <- 1e-6
max_outer <- 100L
cond_limit <- 1e8

run_constrained_gmm <- function(dat, design_mat, h_x_full, link_type) {
  r_vec <- dat[["r"]]
  n <- nrow(dat); k <- ncol(design_mat)

  # Drop intercept column from h_x (first column is all 1s)
  h_x <- h_x_full[, -1, drop = FALSE]
  h_dim <- ncol(h_x)

  init <- rep(0, k)

  g_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    (r_vec / pi_v - 1) * h_x
  }
  G_bar_fn <- function(alpha) {
    g_mat <- g_fn(alpha)
    matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
  }
  Gamma_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    crossprod(h_x * (r_vec * (1 - pi_v) / pi_v), design_mat) / n
  }
  forward_cond <- function(alpha) {
    g_c <- g_fn(alpha)
    W_c <- tryCatch(solve(crossprod(g_c) / n), error = function(e) diag(h_dim))
    Gamma <- Gamma_fn(alpha)
    H <- crossprod(Gamma, W_c %*% Gamma)
    eig <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    max(eig) / max(min(eig), 1e-15)
  }

  # Equality constraint: E[r/pi(alpha) - 1] = 0
  constraint_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    mean(r_vec / pi_v - 1)
  }
  constraint_jac <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    as.vector(colMeans(r_vec * (1 - pi_v) / pi_v * design_mat))
  }

  # SLSQP constrained optimization for fixed W
  solve_constrained <- function(start, W_hat) {
    obj_fn <- function(alpha) {
      G_v <- G_bar_fn(alpha)
      as.numeric(crossprod(G_v, W_hat %*% G_v))
    }
    grad_fn <- function(alpha) {
      G_v <- G_bar_fn(alpha)
      Gamma <- Gamma_fn(alpha)
      2 * as.vector(crossprod(Gamma, W_hat %*% G_v))
    }
    res <- nloptr::slsqp(
      x0 = start,
      fn = obj_fn,
      gr = grad_fn,
      heq = constraint_fn,
      heqjac = constraint_jac,
      control = list(maxeval = 2000, xtol_rel = 1e-10, ftol_rel = 1e-12)
    )
    res$par
  }

  # Step 1: W = I
  estimates <- solve_constrained(init, diag(h_dim))

  # Step 2: Iterative W updates
  for (t in 1:max_outer) {
    g_mat <- g_fn(estimates)
    G_vec <- G_bar_fn(estimates)
    W_hat <- tryCatch(W_fn(g_mat), error = function(e) diag(h_dim))
    Gamma <- Gamma_fn(estimates)
    cg <- 2 * as.vector(crossprod(Gamma, W_hat %*% G_vec))
    grad_norm <- max(abs(cg))
    if (grad_norm < tol) break
    estimates <- solve_constrained(estimates, W_hat)
  }

  # Final diagnostics
  g_mat_f <- g_fn(estimates)
  G_f <- G_bar_fn(estimates)
  W_f <- tryCatch(W_fn(g_mat_f), error = function(e) diag(h_dim))
  final_grad_norm <- max(abs(2 * as.vector(crossprod(Gamma_fn(estimates), W_f %*% G_f))))
  final_cond <- forward_cond(estimates)
  pi_final <- compute_pi_fn(as.vector(design_mat %*% estimates))
  constraint_viol <- mean(r_vec / pi_final - 1)

  list(estimates = estimates, grad_norm = final_grad_norm, cond = final_cond,
       pi_hat = pi_final, constraint_viol = constraint_viol, h_x = h_x)
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

# ============================================================
# Test across 3 settings
# ============================================================
settings <- list(
  list(name = "S4 M3", setting = "setting4", m_idx = 3),
  list(name = "S4 M2", setting = "setting4", m_idx = 2),
  list(name = "S3 M2", setting = "setting3", m_idx = 2)
)

cat("=== Constrained GMM via SLSQP (no intercept in h_x): E[r/pi-1]=0 ===\n\n")

for (cfg in settings) {
  cat(sprintf("========== %s ==========\n", cfg$name))
  data_file <- misspecified_model_all_data_file.list[[cfg$setting]][["miss50"]][[1]]
  all_data <- readRDS(data_file)
  mu_true <- get_mu_true(cfg$setting)

  cat("Pre-building design matrices...\n")
  design_info <- vector("list", n_reps)
  for (rep_i in 1:n_reps) {
    dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
    single_ps <- list(
      formula.list = list(ps_spec[["formula.list"]][[cfg$m_idx]]),
      h_alpha.list = list(ps_spec[["h_alpha.list"]][[cfg$m_idx]]),
      inv_link = ps_spec[["inv_link"]],
      outcome = ps_spec[["outcome"]],
      alpha_init.list = list(NULL),
      optimizer = "L-BFGS-B"
    )
    tryCatch({
      ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
      design_info[[rep_i]] <- list(
        design_mat = ebmr$ps_fit.list[[1]]$design_matrix,
        h_x = ebmr$ps_fit.list[[1]]$h_x,
        link_type = ebmr$ps_fit.list[[1]]$link_type
      )
    }, error = function(e) NULL)
  }
  cat("Done.\n")

  mu_vals <- se_vals <- cond_vals <- grad_vals <- cv_vals <- rep(NA_real_, n_reps)
  # Track pi distribution for y=1 in degenerate reps
  degen_count <- 0

  for (rep_i in 1:n_reps) {
    if (rep_i %% 200 == 0) cat(sprintf("  rep %d...\n", rep_i))
    if (is.null(design_info[[rep_i]])) next
    dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
    di <- design_info[[rep_i]]
    tryCatch({
      res <- run_constrained_gmm(dat, di$design_mat, di$h_x, di$link_type)
      cond_vals[rep_i] <- res$cond
      grad_vals[rep_i] <- res$grad_norm
      cv_vals[rep_i] <- res$constraint_viol

      # Check pi for y=1 in first few degenerate reps
      if (res$cond >= cond_limit && degen_count < 3) {
        degen_count <- degen_count + 1
        y <- dat$y; pi_h <- res$pi_hat
        cat(sprintf("  [Degen rep %d] alpha=%s\n", rep_i,
            paste(round(res$estimates, 3), collapse=", ")))
        cat(sprintf("    pi(y=1): min=%.6f, mean=%.6f, max=%.6f\n",
            min(pi_h[y==1]), mean(pi_h[y==1]), max(pi_h[y==1])))
        cat(sprintf("    pi(y=0): min=%.6f, mean=%.6f, max=%.6f\n",
            min(pi_h[y==0]), mean(pi_h[y==0]), max(pi_h[y==0])))
        cat(sprintf("    constraint viol = %.2e\n", res$constraint_viol))
      }

      # SE uses h_x WITHOUT intercept (same as optimization)
      se_res <- compute_se(dat, res$estimates, di$design_mat, res$h_x, di$link_type, n_val)
      mu_vals[rep_i] <- se_res$mu
      se_vals[rep_i] <- se_res$se
    }, error = function(e) {
      cat(sprintf("  Rep %d ERROR: %s\n", rep_i, e$message))
    })
  }

  valid <- !is.na(mu_vals) & !is.na(se_vals)
  n_degen <- sum(cond_vals[valid] >= cond_limit, na.rm = TRUE)
  n_conv <- sum(grad_vals[valid] < tol, na.rm = TRUE)
  cat(sprintf("  Valid: %d/%d, Degen: %d, Conv: %d\n",
      sum(valid), n_reps, n_degen, n_conv))
  cat(sprintf("  Bias     = %.4f\n", mean(mu_vals[valid]) - mu_true))
  cat(sprintf("  ESD      = %.4f\n", sd(mu_vals[valid])))
  cat(sprintf("  ESE      = %.4f\n", mean(se_vals[valid])))
  cat(sprintf("  ESE/ESD  = %.3f\n", mean(se_vals[valid]) / sd(mu_vals[valid])))
  ci_lo <- mu_vals[valid] - 1.96 * se_vals[valid]
  ci_hi <- mu_vals[valid] + 1.96 * se_vals[valid]
  cat(sprintf("  CP       = %.3f\n", mean(mu_true >= ci_lo & mu_true <= ci_hi)))
  cat(sprintf("  Constraint |E[r/pi-1]|: median=%.2e, max=%.2e\n",
      median(abs(cv_vals[valid])), max(abs(cv_vals[valid]))))
  cat("\n")
}
