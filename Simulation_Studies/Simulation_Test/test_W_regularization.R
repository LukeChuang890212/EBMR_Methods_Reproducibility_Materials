setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

n_val <- 2000
n_reps <- 200
ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
mu_true <- get_mu_true("setting4")
model_idx <- 3
compute_pi_fn <- function(eta) plogis(-eta)
tol <- 1e-6
max_outer <- 100

# SE2 computation
compute_se2_at_alpha <- function(dat, alpha_hat, design_mat, h_x, link_type, n, W_fn_se) {
  r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
  pi_hat <- compute_pi_fn(as.vector(design_mat %*% alpha_hat))
  mu_hat <- mean(r_vec * y_vec / pi_hat)
  g_mat <- (r_vec / pi_hat - 1) * h_x
  W_hat <- tryCatch(W_fn_se(g_mat), error = function(e) NULL)
  if (is.null(W_hat)) return(list(mu = mu_hat, se2 = NA, cond = NA))
  cf <- r_vec * (1 - pi_hat) / pi_hat
  Gamma_hat <- crossprod(h_x * cf, design_mat) / n
  GtWG <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
  eig <- eigen(GtWG, symmetric = TRUE, only.values = TRUE)$values
  cond_val <- max(eig) / max(min(eig), 1e-15)
  k <- length(alpha_hat); h <- ncol(h_x)
  GtW <- crossprod(Gamma_hat, W_hat)
  G_vec <- colMeans(g_mat)
  dc_deta <- cf; dc_X <- dc_deta * design_mat
  R_mat <- matrix(0, h * k, k)
  for (j in 1:k) {
    Xj_h <- design_mat[, j] * h_x
    block <- crossprod(dc_X, Xj_h) / n
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
  H2_mat <- GtWG + GW_kron_Ik %*% R_mat - GW_kron_GtW %*% S_mat
  H2_cond <- tryCatch({
    ev <- eigen(H2_mat, symmetric = FALSE, only.values = TRUE)$values
    max(Mod(ev)) / max(min(Mod(ev)), 1e-15)
  }, error = function(e) Inf)
  if (!is.finite(H2_cond) || H2_cond > 1e15)
    return(list(mu = mu_hat, se2 = NA, cond = cond_val))
  H2_inv <- tryCatch(solve(H2_mat), error = function(e) NULL)
  if (is.null(H2_inv)) return(list(mu = mu_hat, se2 = NA, cond = cond_val))
  GtW_g <- GtW %*% t(g_mat)
  g_WG <- as.vector(g_mat %*% WG)
  cf_hxWG <- cf * as.vector(h_x %*% WG)
  GammaT_WG <- t(design_mat * cf_hxWG)
  Q_mat <- GtW_g + GammaT_WG - sweep(GtW_g, 2, g_WG, `*`)
  psi2 <- -(H2_inv %*% Q_mat)
  if (!is.null(link_type) && link_type == "logistic_complement") {
    dot_pi <- -design_mat * (pi_hat * (1 - pi_hat))
  } else {
    dot_pi <- design_mat * (pi_hat * (1 - pi_hat))
  }
  ry_ps_inv2 <- as.vector(r_vec * y_vec * (pi_hat^(-2)))
  H_alpha <- colMeans(dot_pi * ry_ps_inv2)
  mu_iid <- as.vector(t(r_vec/pi_hat*y_vec) - t(H_alpha) %*% psi2)
  se2 <- sqrt(var(mu_iid) / n)
  list(mu = mu_hat, se2 = se2, cond = cond_val)
}

# Regularized W: W = solve(g'g/n + lambda * I)
run_reg_gmm <- function(dat, design_mat, h_x, link_type, lambda) {
  r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
  n <- nrow(dat); k <- ncol(design_mat); h_dim <- ncol(h_x)

  g_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    (r_vec / pi_v - 1) * h_x
  }
  Gamma_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    crossprod(h_x * (r_vec * (1 - pi_v) / pi_v), design_mat) / n
  }
  W_reg <- function(g_mat) solve(crossprod(g_mat) / n + lambda * diag(h_dim))

  # NR inner with cond guard
  nr_inner <- function(alpha, W_hat, cond_limit = 1e8) {
    for (iter in 1:50) {
      g_i <- g_fn(alpha)
      G_i <- matrix(.colMeans(g_i, n, h_dim), h_dim, 1)
      Gamma_i <- Gamma_fn(alpha)
      gv <- 2 * as.vector(crossprod(Gamma_i, W_hat %*% G_i))
      if (max(abs(gv)) < tol) break
      GtWG <- crossprod(Gamma_i, W_hat %*% Gamma_i)
      H_mat <- GtWG + GtWG
      direction <- tryCatch(-solve(H_mat, gv),
                            error = function(e) -gv * 2.0 / max(abs(gv)))
      if (sqrt(sum(direction^2)) > 2.0)
        direction <- direction * 2.0 / sqrt(sum(direction^2))
      obj_cur <- as.numeric(crossprod(G_i, W_hat %*% G_i))
      accepted <- FALSE; step <- 1.0
      for (ls in 1:20) {
        cand <- alpha + step * direction
        g_cand <- g_fn(cand)
        G_cand <- matrix(.colMeans(g_cand, n, h_dim), h_dim, 1)
        obj_new <- as.numeric(crossprod(G_cand, W_hat %*% G_cand))
        if (is.finite(obj_new) && obj_new < obj_cur - 1e-4 * step * sum(gv * direction)) {
          Gamma_c <- Gamma_fn(cand)
          H_c <- crossprod(Gamma_c, W_hat %*% Gamma_c)
          eig_c <- eigen(H_c, symmetric = TRUE, only.values = TRUE)$values
          cc <- max(eig_c) / max(min(eig_c), 1e-15)
          if (cc < cond_limit) { accepted <- TRUE; break }
        }
        step <- step * 0.5
      }
      if (accepted) alpha <- cand else break
    }
    alpha
  }

  # Step 1: W = I
  est <- nr_inner(rep(0, k), diag(h_dim))

  # Outer iterations with regularized W
  best_grad <- Inf; best_est <- est
  for (t in 1:max_outer) {
    g_mat <- g_fn(est)
    W_hat <- tryCatch(W_reg(g_mat), error = function(e) diag(h_dim))
    Gamma <- Gamma_fn(est)
    G_vec <- matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
    grad <- max(abs(2 * as.vector(crossprod(Gamma, W_hat %*% G_vec))))
    if (grad < best_grad) { best_grad <- grad; best_est <- est }
    if (grad < tol) break
    est <- nr_inner(est, W_hat)
  }

  # Compute cond at best solution using unregularized W for fair comparison
  g_f <- g_fn(best_est)
  W_unreg <- tryCatch(solve(crossprod(g_f) / n), error = function(e) diag(h_dim))
  Gamma_f <- Gamma_fn(best_est)
  GtWG <- crossprod(Gamma_f, W_unreg %*% Gamma_f)
  eig <- eigen(GtWG, symmetric = TRUE, only.values = TRUE)$values
  cond_unreg <- max(eig) / max(min(eig), 1e-15)

  list(estimates = best_est, grad_norm = best_grad, cond = cond_unreg)
}

# Test multiple lambda values
for (lambda in c(0.001, 0.01, 0.1, 1.0)) {
  cat(sprintf("\n=== Lambda = %.3f ===\n", lambda))
  mu_vals <- se2_vals <- cond_vals <- rep(NA_real_, n_reps)
  n_degen <- 0

  W_fn_unreg <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

  for (rep_i in 1:n_reps) {
    dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
    single_ps <- list(
      formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
      h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
      inv_link = ps_spec[["inv_link"]],
      outcome = ps_spec[["outcome"]],
      alpha_init.list = list(NULL),
      optimizer = "L-BFGS-B"
    )
    tryCatch({
      ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn_unreg)
      design_mat <- ebmr$ps_fit.list[[1]]$design_matrix
      h_x <- ebmr$ps_fit.list[[1]]$h_x
      link_type <- ebmr$ps_fit.list[[1]]$link_type

      res <- run_reg_gmm(dat, design_mat, h_x, link_type, lambda)
      cond_vals[rep_i] <- res$cond
      if (res$cond >= 1e8) n_degen <- n_degen + 1

      # SE2 with unregularized W (standard formula)
      se_res <- compute_se2_at_alpha(dat, res$estimates, design_mat, h_x, link_type, n_val, W_fn_unreg)
      mu_vals[rep_i] <- se_res$mu
      se2_vals[rep_i] <- se_res$se2
    }, error = function(e) NULL)
  }

  valid <- !is.na(mu_vals) & !is.na(se2_vals)
  well_cond <- valid & cond_vals < 1e8
  cat(sprintf("  Valid: %d, Degenerate (cond>=1e8): %d\n", sum(valid), n_degen))
  cat(sprintf("  All valid:  Bias=%.4f, ESD=%.4f, ESE2=%.4f, ESE2/ESD=%.3f\n",
      mean(mu_vals[valid]) - mu_true, sd(mu_vals[valid]),
      mean(se2_vals[valid]), mean(se2_vals[valid]) / sd(mu_vals[valid])))
  if (sum(well_cond) > 0 && sum(well_cond) < sum(valid)) {
    cat(sprintf("  Cond<1e8:   Bias=%.4f, ESD=%.4f, ESE2=%.4f, ESE2/ESD=%.3f (n=%d)\n",
        mean(mu_vals[well_cond]) - mu_true, sd(mu_vals[well_cond]),
        mean(se2_vals[well_cond]), mean(se2_vals[well_cond]) / sd(mu_vals[well_cond]),
        sum(well_cond)))
  }
}
