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
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
compute_pi_fn <- function(eta) plogis(-eta)
tol <- 1e-6
max_outer <- 100
cond_limit <- 1e8

# SGD-based iterative GMM with cond guard
# Uses mini-batch stochastic gradient for the inner loop.
# The stochasticity helps escape degenerate regions because the noisy gradient
# doesn't perfectly align with the direction toward degeneracy.
run_sgd_gmm <- function(dat, design_mat, h_x, link_type) {
  r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
  n <- nrow(dat); k <- ncol(design_mat); h_dim <- ncol(h_x)

  # Full-batch g and Gamma
  g_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    (r_vec / pi_v - 1) * h_x
  }
  Gamma_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    crossprod(h_x * (r_vec * (1 - pi_v) / pi_v), design_mat) / n
  }
  get_cond <- function(alpha, W_hat) {
    Gamma <- Gamma_fn(alpha)
    H <- crossprod(Gamma, W_hat %*% Gamma)
    eig <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    max(eig) / max(min(eig), 1e-15)
  }

  # Mini-batch g and Gamma
  g_batch <- function(alpha, idx) {
    pi_v <- compute_pi_fn(as.vector(design_mat[idx, , drop = FALSE] %*% alpha))
    (r_vec[idx] / pi_v - 1) * h_x[idx, , drop = FALSE]
  }
  Gamma_batch <- function(alpha, idx) {
    pi_v <- compute_pi_fn(as.vector(design_mat[idx, , drop = FALSE] %*% alpha))
    nb <- length(idx)
    crossprod(h_x[idx, , drop = FALSE] * (r_vec[idx] * (1 - pi_v) / pi_v),
              design_mat[idx, , drop = FALSE]) / nb
  }

  # SGD inner loop: minimize Q(alpha) = G'W G for fixed W
  # Uses mini-batch gradient with decreasing step size
  sgd_inner <- function(start, W_hat, batch_size = 256, max_iter = 500, lr0 = 0.1) {
    alpha <- start
    best_alpha <- alpha
    best_obj <- Inf

    for (iter in 1:max_iter) {
      # Mini-batch gradient
      idx <- sample.int(n, min(batch_size, n))
      g_b <- g_batch(alpha, idx)
      G_b <- matrix(colMeans(g_b), h_dim, 1)
      Gamma_b <- Gamma_batch(alpha, idx)
      grad <- 2 * as.vector(crossprod(Gamma_b, W_hat %*% G_b))

      # Learning rate schedule: lr0 / (1 + iter/50)
      lr <- lr0 / (1 + iter / 50)

      # Gradient step with norm clamp
      step_vec <- -lr * grad
      sn <- sqrt(sum(step_vec^2))
      if (sn > 1.0) step_vec <- step_vec * (1.0 / sn)

      cand <- alpha + step_vec

      # Cond guard on full batch
      if (iter %% 10 == 0 || iter <= 5) {
        cond_c <- get_cond(cand, W_hat)
        if (cond_c >= cond_limit) {
          # Reject step, try smaller
          step_vec <- step_vec * 0.1
          cand <- alpha + step_vec
          cond_c <- get_cond(cand, W_hat)
          if (cond_c >= cond_limit) next  # skip this step entirely
        }
      }

      alpha <- cand

      # Track best on full batch every 20 iters
      if (iter %% 20 == 0) {
        g_full <- g_fn(alpha)
        G_full <- matrix(.colMeans(g_full, n, h_dim), h_dim, 1)
        obj <- as.numeric(crossprod(G_full, W_hat %*% G_full))
        if (is.finite(obj) && obj < best_obj) {
          best_obj <- obj
          best_alpha <- alpha
        }
        # Check convergence
        Gamma_full <- Gamma_fn(alpha)
        grad_full <- max(abs(2 * as.vector(crossprod(Gamma_full, W_hat %*% G_full))))
        if (grad_full < tol) break
      }
    }
    best_alpha
  }

  # Full-batch GD inner loop (deterministic, for comparison/refinement)
  gd_inner <- function(start, W_hat, max_iter = 200, lr0 = 0.05) {
    alpha <- start
    best_alpha <- alpha
    best_obj <- Inf

    for (iter in 1:max_iter) {
      g_mat <- g_fn(alpha)
      G_vec <- matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
      Gamma <- Gamma_fn(alpha)
      grad <- 2 * as.vector(crossprod(Gamma, W_hat %*% G_vec))

      if (max(abs(grad)) < tol) break

      lr <- lr0 / (1 + iter / 100)
      step_vec <- -lr * grad
      sn <- sqrt(sum(step_vec^2))
      if (sn > 1.0) step_vec <- step_vec * (1.0 / sn)

      cand <- alpha + step_vec

      # Cond guard
      cond_c <- get_cond(cand, W_hat)
      if (cond_c >= cond_limit) {
        # Backtrack step size
        for (bt in 1:10) {
          step_vec <- step_vec * 0.5
          cand <- alpha + step_vec
          cond_c <- get_cond(cand, W_hat)
          if (cond_c < cond_limit) break
        }
        if (cond_c >= cond_limit) next
      }

      alpha <- cand
      obj <- as.numeric(crossprod(G_vec, W_hat %*% G_vec))
      if (is.finite(obj) && obj < best_obj) {
        best_obj <- obj
        best_alpha <- alpha
      }
    }
    best_alpha
  }

  # Step 1: W = I, use GD
  est <- gd_inner(rep(0, k), diag(h_dim), max_iter = 500, lr0 = 0.1)

  # Check if solution is good enough, refine with SGD
  g_mat <- g_fn(est)
  G_vec <- matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
  grad0 <- max(abs(2 * as.vector(crossprod(Gamma_fn(est), G_vec))))
  if (grad0 >= tol) {
    est2 <- sgd_inner(est, diag(h_dim), batch_size = 256, max_iter = 1000, lr0 = 0.05)
    g2 <- g_fn(est2)
    G2 <- matrix(.colMeans(g2, n, h_dim), h_dim, 1)
    grad2 <- max(abs(2 * as.vector(crossprod(Gamma_fn(est2), G2))))
    if (grad2 < grad0) est <- est2
  }

  # Outer loop: iterative W updates
  best_grad <- Inf; best_est <- est

  for (t in 1:max_outer) {
    g_mat <- g_fn(est)
    W_hat <- tryCatch(W_fn(g_mat), error = function(e) diag(h_dim))

    # Check cond before proceeding
    cond_cur <- get_cond(est, W_hat)
    if (cond_cur >= cond_limit) {
      # W update causes degeneracy — try GD with smaller steps
      est_retry <- gd_inner(est, W_hat, max_iter = 500, lr0 = 0.01)
      cond_retry <- get_cond(est_retry, W_hat)
      if (cond_retry >= cond_limit) {
        # Also try SGD
        est_retry <- sgd_inner(est, W_hat, batch_size = 128, max_iter = 500, lr0 = 0.01)
        cond_retry <- get_cond(est_retry, W_hat)
      }
      if (cond_retry >= cond_limit) {
        est <- best_est
        break
      }
      est <- est_retry
      next
    }

    Gamma <- Gamma_fn(est)
    G_vec <- matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
    grad <- max(abs(2 * as.vector(crossprod(Gamma, W_hat %*% G_vec))))

    if (grad < best_grad) { best_grad <- grad; best_est <- est }
    if (grad < tol) break

    # Inner optimization with GD (safer than NR for degenerate landscape)
    new_est <- gd_inner(est, W_hat, max_iter = 300, lr0 = 0.05)

    # Check cond after
    g_new <- g_fn(new_est)
    W_new <- tryCatch(W_fn(g_new), error = function(e) diag(h_dim))
    new_cond <- get_cond(new_est, W_new)
    new_grad <- max(abs(2 * as.vector(crossprod(Gamma_fn(new_est), W_new %*%
                    matrix(.colMeans(g_new, n, h_dim), h_dim, 1)))))

    if (new_cond < cond_limit) {
      est <- new_est
      if (new_grad < best_grad) { best_grad <- new_grad; best_est <- new_est }
    } else {
      est <- best_est
      break
    }
  }

  # Final solution cond
  g_f <- g_fn(best_est)
  W_f <- tryCatch(W_fn(g_f), error = function(e) diag(h_dim))
  final_cond <- get_cond(best_est, W_f)

  list(estimates = best_est, grad_norm = best_grad, cond = final_cond)
}

# SE2 computation
compute_se2 <- function(dat, alpha_hat, design_mat, h_x, link_type, n) {
  r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
  pi_hat <- compute_pi_fn(as.vector(design_mat %*% alpha_hat))
  mu_hat <- mean(r_vec * y_vec / pi_hat)
  g_mat <- (r_vec / pi_hat - 1) * h_x
  W_hat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) NULL)
  if (is.null(W_hat)) return(list(mu = mu_hat, se2 = NA, cond = NA))
  cf <- r_vec * (1 - pi_hat) / pi_hat
  Gamma_hat <- crossprod(h_x * cf, design_mat) / n
  GtWG <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
  eig <- eigen(GtWG, symmetric = TRUE, only.values = TRUE)$values
  cond_val <- max(eig) / max(min(eig), 1e-15)
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
  H2_mat <- GtWG + GW_kron_Ik %*% R_mat - GW_kron_GtW %*% S_mat
  H2_cond <- tryCatch({
    ev <- eigen(H2_mat, symmetric = FALSE, only.values = TRUE)$values
    max(Mod(ev)) / max(min(Mod(ev)), 1e-15)
  }, error = function(e) Inf)
  if (!is.finite(H2_cond) || H2_cond > 1e15) return(list(mu = mu_hat, se2 = NA, cond = cond_val))
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

# ========== Test M3 ==========
cat("=== Model 3: SGD-based GMM ===\n")
m_idx <- 3
mu_vals <- se2_vals <- cond_vals <- grad_vals <- rep(NA_real_, n_reps)

for (rep_i in 1:n_reps) {
  if (rep_i %% 50 == 0) cat(sprintf("  rep %d...\n", rep_i))
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]

  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[m_idx]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[m_idx]]),
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL),
    optimizer = "L-BFGS-B"
  )
  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
    design_mat <- ebmr$ps_fit.list[[1]]$design_matrix
    h_x <- ebmr$ps_fit.list[[1]]$h_x
    link_type <- ebmr$ps_fit.list[[1]]$link_type

    res <- run_sgd_gmm(dat, design_mat, h_x, link_type)
    cond_vals[rep_i] <- res$cond
    grad_vals[rep_i] <- res$grad_norm

    se_res <- compute_se2(dat, res$estimates, design_mat, h_x, link_type, n_val)
    mu_vals[rep_i] <- se_res$mu
    se2_vals[rep_i] <- se_res$se2
  }, error = function(e) {
    cat(sprintf("  Rep %d ERROR: %s\n", rep_i, e$message))
  })
}

valid <- !is.na(mu_vals) & !is.na(se2_vals)
n_degen <- sum(cond_vals >= cond_limit, na.rm = TRUE)
n_conv <- sum(grad_vals < 1e-6, na.rm = TRUE)
cat(sprintf("\n  Valid: %d/%d, Degenerate: %d, Converged: %d\n",
    sum(valid), n_reps, n_degen, n_conv))
cat(sprintf("  Bias  = %.4f\n", mean(mu_vals[valid]) - mu_true))
cat(sprintf("  ESD   = %.4f\n", sd(mu_vals[valid])))
cat(sprintf("  ESE2  = %.4f\n", mean(se2_vals[valid])))
cat(sprintf("  ESE2/ESD = %.3f\n", mean(se2_vals[valid]) / sd(mu_vals[valid])))

well_cond <- valid & cond_vals < cond_limit
if (sum(well_cond) > 1 && sum(well_cond) < sum(valid)) {
  cat(sprintf("  [Cond<1e8] n=%d, Bias=%.4f, ESD=%.4f, ESE2=%.4f, ESE2/ESD=%.3f\n",
      sum(well_cond), mean(mu_vals[well_cond]) - mu_true, sd(mu_vals[well_cond]),
      mean(se2_vals[well_cond]), mean(se2_vals[well_cond]) / sd(mu_vals[well_cond])))
}
