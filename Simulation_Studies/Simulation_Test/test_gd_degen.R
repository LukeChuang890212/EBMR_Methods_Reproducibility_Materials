setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

n_val <- 2000
ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
mu_true <- get_mu_true("setting4")
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
compute_pi_fn <- function(eta) plogis(-eta)
tol <- 1e-6
cond_limit <- 1e8

# Known degenerate reps from previous analysis
degen_reps <- c(34, 35, 38, 68, 109, 112, 133, 143, 155, 156, 161, 166, 177, 188, 190, 194)
m_idx <- 3

# Full-batch gradient descent with cond guard
# Key difference from NR: step direction = -gradient (steepest descent),
# much more conservative than Newton direction
run_gd_gmm <- function(dat, design_mat, h_x, link_type) {
  r_vec <- dat[["r"]]
  n <- nrow(dat); k <- ncol(design_mat); h_dim <- ncol(h_x)

  g_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    (r_vec / pi_v - 1) * h_x
  }
  Gamma_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    crossprod(h_x * (r_vec * (1 - pi_v) / pi_v), design_mat) / n
  }
  check_cond <- function(alpha, W_hat) {
    Gamma <- Gamma_fn(alpha)
    H <- crossprod(Gamma, W_hat %*% Gamma)
    eig <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    max(eig) / max(min(eig), 1e-15)
  }
  # Also check "forward cond": cond with W computed from g at candidate
  forward_cond <- function(alpha) {
    g_c <- g_fn(alpha)
    W_c <- tryCatch(solve(crossprod(g_c) / n), error = function(e) diag(h_dim))
    check_cond(alpha, W_c)
  }

  # GD inner: minimize Q(alpha) = G'WG for fixed W
  # Uses gradient descent with Barzilai-Borwein step size + cond guard
  gd_inner <- function(start, W_hat, max_iter = 1000) {
    alpha <- start
    best_alpha <- alpha; best_obj <- Inf
    prev_grad <- NULL; prev_alpha <- NULL

    for (iter in 1:max_iter) {
      g_mat <- g_fn(alpha)
      G_vec <- matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
      Gamma <- Gamma_fn(alpha)
      grad <- 2 * as.vector(crossprod(Gamma, W_hat %*% G_vec))
      obj <- as.numeric(crossprod(G_vec, W_hat %*% G_vec))

      if (is.finite(obj) && obj < best_obj) { best_obj <- obj; best_alpha <- alpha }
      if (max(abs(grad)) < tol) break

      # Barzilai-Borwein step size (adaptive)
      if (!is.null(prev_grad)) {
        s <- alpha - prev_alpha
        y <- grad - prev_grad
        sy <- sum(s * y)
        if (abs(sy) > 1e-15) {
          lr <- abs(sum(s * s) / sy)  # BB1
          lr <- min(lr, 1.0)  # cap
        } else {
          lr <- 0.01
        }
      } else {
        lr <- 0.01
      }

      prev_grad <- grad; prev_alpha <- alpha
      step_vec <- -lr * grad

      # Clamp step norm
      sn <- sqrt(sum(step_vec^2))
      if (sn > 2.0) step_vec <- step_vec * (2.0 / sn)

      # Backtracking with cond guard
      accepted <- FALSE
      for (bt in 1:20) {
        cand <- alpha + step_vec
        g_cand <- g_fn(cand)
        G_cand <- matrix(.colMeans(g_cand, n, h_dim), h_dim, 1)
        obj_new <- as.numeric(crossprod(G_cand, W_hat %*% G_cand))
        if (is.finite(obj_new) && obj_new < obj - 1e-4 * sum(grad * step_vec)) {
          cc <- check_cond(cand, W_hat)
          if (cc < cond_limit) { accepted <- TRUE; break }
        }
        step_vec <- step_vec * 0.5
      }
      if (accepted) alpha <- cand else break
    }
    best_alpha
  }

  # Step 1: W = I
  est <- gd_inner(rep(0, k), diag(h_dim))

  # Outer loop: W updates with forward cond guard
  best_grad <- Inf; best_est <- est

  for (t in 1:100) {
    g_mat <- g_fn(est)
    W_hat <- tryCatch(W_fn(g_mat), error = function(e) diag(h_dim))
    cond_cur <- check_cond(est, W_hat)

    if (cond_cur >= cond_limit) {
      # Don't give up — try GD from current est with smaller steps
      est_retry <- gd_inner(est, W_hat, max_iter = 500)
      cond_retry <- check_cond(est_retry, W_hat)
      if (cond_retry >= cond_limit) { est <- best_est; break }
      est <- est_retry
      next
    }

    Gamma <- Gamma_fn(est)
    G_vec <- matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
    grad <- max(abs(2 * as.vector(crossprod(Gamma, W_hat %*% G_vec))))

    if (grad < best_grad) { best_grad <- grad; best_est <- est }
    if (grad < tol) break

    est <- gd_inner(est, W_hat)
  }

  final_fwd_cond <- forward_cond(best_est)
  list(estimates = best_est, grad_norm = best_grad, cond = final_fwd_cond)
}

# Test on degenerate reps
cat("=== GD-based GMM on known degenerate reps ===\n")
for (rep_i in degen_reps) {
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

    res <- run_gd_gmm(dat, design_mat, h_x, link_type)
    pi_hat <- compute_pi_fn(as.vector(design_mat %*% res$estimates))
    r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
    mu <- mean(r_vec * y_vec / pi_hat)

    cat(sprintf("  Rep %3d: grad=%.2e, cond=%.2e, mu=%.4f, alpha[2]=%.3f, #pi>0.99=%d\n",
        rep_i, res$grad_norm, res$cond, mu, res$estimates[2], sum(pi_hat > 0.99)))
  }, error = function(e) {
    cat(sprintf("  Rep %3d: ERROR %s\n", rep_i, e$message))
  })
}
