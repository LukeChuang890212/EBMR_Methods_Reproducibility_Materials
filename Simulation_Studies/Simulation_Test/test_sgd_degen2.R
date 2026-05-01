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

degen_reps <- c(34, 112, 161)  # 3 reps to trace
m_idx <- 3

run_sgd_gmm_traced <- function(dat, design_mat, h_x, link_type,
                                batch_size = 256, max_inner = 2000, lr0 = 0.05) {
  r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
  n <- nrow(dat); k <- ncol(design_mat); h_dim <- ncol(h_x)

  g_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    (r_v / pi_v - 1) * h_x
  }
  # Fix: use r_vec not r_v
  g_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    (r_vec / pi_v - 1) * h_x
  }
  Gamma_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    crossprod(h_x * (r_vec * (1 - pi_v) / pi_v), design_mat) / n
  }
  obj_fn <- function(alpha, W_hat) {
    g_mat <- g_fn(alpha)
    G <- matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
    as.numeric(crossprod(G, W_hat %*% G))
  }
  grad_fn <- function(alpha, W_hat) {
    g_mat <- g_fn(alpha)
    G <- matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
    Gamma <- Gamma_fn(alpha)
    2 * as.vector(crossprod(Gamma, W_hat %*% G))
  }
  check_cond <- function(alpha) {
    g_c <- g_fn(alpha)
    W_c <- tryCatch(solve(crossprod(g_c) / n), error = function(e) diag(h_dim))
    Gamma_c <- Gamma_fn(alpha)
    H_c <- crossprod(Gamma_c, W_c %*% Gamma_c)
    eig_c <- eigen(H_c, symmetric = TRUE, only.values = TRUE)$values
    max(eig_c) / max(min(eig_c), 1e-15)
  }

  sgd_grad <- function(alpha, W_hat, idx) {
    nb <- length(idx)
    dm_b <- design_mat[idx, , drop = FALSE]
    hx_b <- h_x[idx, , drop = FALSE]
    r_b <- r_vec[idx]
    pi_v <- compute_pi_fn(as.vector(dm_b %*% alpha))
    g_b <- (r_b / pi_v - 1) * hx_b
    G_b <- matrix(colMeans(g_b), h_dim, 1)
    cf_b <- r_b * (1 - pi_v) / pi_v
    Gamma_b <- crossprod(hx_b * cf_b, dm_b) / nb
    2 * as.vector(crossprod(Gamma_b, W_hat %*% G_b))
  }

  sgd_inner <- function(start, W_hat, label = "") {
    alpha <- start
    best_alpha <- alpha
    best_obj <- obj_fn(alpha, W_hat)

    for (iter in 1:max_inner) {
      idx <- sample.int(n, min(batch_size, n))
      sg <- sgd_grad(alpha, W_hat, idx)
      lr <- lr0 / sqrt(iter)
      step_vec <- -lr * sg
      sn <- sqrt(sum(step_vec^2))
      if (sn > 1.0) step_vec <- step_vec * (1.0 / sn)
      cand <- alpha + step_vec

      if (iter %% 50 == 0) {
        Gamma_c <- Gamma_fn(cand)
        H_c <- crossprod(Gamma_c, W_hat %*% Gamma_c)
        eig_c <- eigen(H_c, symmetric = TRUE, only.values = TRUE)$values
        cc <- max(eig_c) / max(min(eig_c), 1e-15)
        if (cc >= cond_limit) {
          step_vec <- step_vec * 0.1
          cand <- alpha + step_vec
        }
      }
      alpha <- cand

      if (iter %% 100 == 0) {
        obj <- obj_fn(alpha, W_hat)
        if (is.finite(obj) && obj < best_obj) {
          best_obj <- obj; best_alpha <- alpha
        }
        gf <- max(abs(grad_fn(alpha, W_hat)))
        if (gf < tol) break
      }
    }

    # Polish
    alpha <- best_alpha
    for (iter in 1:100) {
      gf <- grad_fn(alpha, W_hat)
      if (max(abs(gf)) < tol) break
      lr <- 0.01 / (1 + iter / 20)
      step_vec <- -lr * gf
      sn <- sqrt(sum(step_vec^2))
      if (sn > 0.5) step_vec <- step_vec * (0.5 / sn)
      cand <- alpha + step_vec
      obj_new <- obj_fn(cand, W_hat)
      if (is.finite(obj_new) && obj_new < obj_fn(alpha, W_hat)) {
        alpha <- cand
      } else break
    }
    obj_final <- obj_fn(alpha, W_hat)
    if (obj_final < best_obj) best_alpha <- alpha

    gf <- max(abs(grad_fn(best_alpha, W_hat)))
    cat(sprintf("    %s: grad=%.2e, obj=%.2e, alpha[2]=%.3f\n",
        label, gf, obj_fn(best_alpha, W_hat), best_alpha[2]))
    best_alpha
  }

  # Step 1: W = I
  cat("  Step1 (W=I):\n")
  est <- sgd_inner(rep(0, k), diag(h_dim), "SGD(W=I)")
  cat(sprintf("    forward_cond=%.2e\n", check_cond(est)))

  # Outer W updates with trace
  best_grad <- Inf; best_est <- est
  for (t in 1:50) {
    g_mat <- g_fn(est)
    W_hat <- tryCatch(W_fn(g_mat), error = function(e) diag(h_dim))
    grad <- max(abs(grad_fn(est, W_hat)))
    if (grad < best_grad) { best_grad <- grad; best_est <- est }

    fc <- check_cond(est)
    cat(sprintf("  Outer %2d: grad(W_new)=%.2e, fwd_cond=%.2e", t, grad, fc))

    if (grad < tol) { cat(" CONVERGED\n"); break }

    if (fc >= cond_limit) {
      cat(" DEGEN -> retry SGD\n")
      est_retry <- sgd_inner(est, W_hat, sprintf("retry_t%d", t))
      fc_retry <- check_cond(est_retry)
      cat(sprintf("    retry fwd_cond=%.2e\n", fc_retry))
      if (fc_retry >= cond_limit) {
        cat("    retry FAILED, stop\n")
        est <- best_est; break
      }
      est <- est_retry
      next
    }

    cat("\n")
    est <- sgd_inner(est, W_hat, sprintf("SGD_t%d", t))
  }

  final_cond <- check_cond(best_est)
  final_grad <- max(abs(grad_fn(best_est,
    tryCatch(W_fn(g_fn(best_est)), error = function(e) diag(h_dim)))))
  list(estimates = best_est, grad_norm = final_grad, cond = final_cond)
}

for (rep_i in degen_reps) {
  cat(sprintf("\n========== Rep %d ==========\n", rep_i))
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
    set.seed(rep_i)
    res <- run_sgd_gmm_traced(dat, design_mat, h_x, link_type)
    cat(sprintf("  FINAL: grad=%.2e, cond=%.2e, alpha=%s\n",
        res$grad_norm, res$cond, paste(round(res$estimates, 3), collapse=", ")))
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
  })
}
