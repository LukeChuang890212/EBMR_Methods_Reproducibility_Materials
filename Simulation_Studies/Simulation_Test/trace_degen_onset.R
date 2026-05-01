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
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
model_idx <- 3
compute_pi_fn <- function(eta) plogis(-eta)
tol <- 1e-6

# Known degenerate reps
degen_reps <- c(34, 35, 38, 68, 109, 112, 133, 143, 155, 156, 161, 166, 177, 188, 190, 194)

for (rep_i in degen_reps[1:6]) {
  cat(sprintf("\n========== Rep %d ==========\n", rep_i))
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
  r_vec <- dat[["r"]]; n <- n_val

  # Get design_mat and h_x
  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL),
    optimizer = "L-BFGS-B"
  )
  ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
  design_mat <- ebmr$ps_fit.list[[1]]$design_matrix
  h_x <- ebmr$ps_fit.list[[1]]$h_x
  k <- ncol(design_mat); h_dim <- ncol(h_x)

  g_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    (r_vec / pi_v - 1) * h_x
  }
  Gamma_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    crossprod(h_x * (r_vec * (1 - pi_v) / pi_v), design_mat) / n
  }
  cond_fn <- function(alpha) {
    g_mat <- g_fn(alpha)
    W_hat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(h_dim))
    Gamma <- Gamma_fn(alpha)
    H <- crossprod(Gamma, W_hat %*% Gamma)
    eig <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    max(eig) / max(min(eig), 1e-15)
  }

  # Manual CNR: trace each outer iteration
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

  # Outer iteration 0: W = I
  estimates <- rep(0, k)
  cat(sprintf("  Outer 0 (W=I): alpha=%s\n", paste(round(estimates, 3), collapse=", ")))

  estimates <- nr_inner(estimates, diag(h_dim))
  pi_v <- compute_pi_fn(as.vector(design_mat %*% estimates))
  g_mat <- g_fn(estimates)
  G_vec <- matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
  W_hat <- diag(h_dim)
  Gamma <- Gamma_fn(estimates)
  grad <- max(abs(2 * as.vector(crossprod(Gamma, W_hat %*% G_vec))))
  cn <- cond_fn(estimates)
  cat(sprintf("  After NR(W=I): alpha=%s\n", paste(round(estimates, 3), collapse=", ")))
  cat(sprintf("    grad=%.2e, cond=%.2e, pi:[%.4f,%.4f], mean=%.4f, #pi>0.99=%d\n",
      grad, cn, min(pi_v), max(pi_v), mean(pi_v), sum(pi_v > 0.99)))

  # Outer iterations with W updates
  best_grad <- grad; best_est <- estimates
  for (t in 1:20) {
    g_mat <- g_fn(estimates)
    W_hat <- tryCatch(W_fn(g_mat), error = function(e) diag(h_dim))
    Gamma <- Gamma_fn(estimates)
    G_vec <- matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
    grad <- max(abs(2 * as.vector(crossprod(Gamma, W_hat %*% G_vec))))

    if (grad < best_grad) { best_grad <- grad; best_est <- estimates }
    if (grad < tol) {
      cat(sprintf("  Outer %d: CONVERGED, grad=%.2e\n", t, grad))
      break
    }

    old_est <- estimates
    estimates <- nr_inner(estimates, W_hat)
    pi_v <- compute_pi_fn(as.vector(design_mat %*% estimates))
    cn <- cond_fn(estimates)

    # Check if NR made any progress
    moved <- max(abs(estimates - old_est))

    cat(sprintf("  Outer %2d: alpha[2]=%.3f, grad=%.2e, cond=%.2e, pi_max=%.6f, #pi>0.99=%d, moved=%.4f\n",
        t, estimates[2], grad, cn, max(pi_v), sum(pi_v > 0.99), moved))

    if (moved < 1e-10) {
      cat(sprintf("    STALLED (no movement)\n"))
      break
    }
  }

  cat(sprintf("  Final best: grad=%.2e, alpha=%s\n",
      best_grad, paste(round(best_est, 3), collapse=", ")))
}
