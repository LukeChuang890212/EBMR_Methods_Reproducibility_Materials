setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})

n_val <- 2000
n_reps <- 200
ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
mu_true <- get_mu_true("setting4")

formula_j <- ps_spec[["formula.list"]][[3]]
h_alpha_vars <- ps_spec[["h_alpha.list"]][[3]]

run_method <- function(method_name, inner_solver) {
  mu_vals <- rep(NA_real_, n_reps)
  se1_vals <- rep(NA_real_, n_reps)
  alpha_mat <- matrix(NA, n_reps, 4)
  conv_vals <- rep(NA, n_reps)
  grad_vals <- rep(NA_real_, n_reps)

  for (rep_i in 1:n_reps) {
    dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
    tryCatch({
      r_vec <- dat[["r"]]
      y_vec <- dat[["y"]]
      design_mat <- model.matrix(formula_j, data=dat)
      h_alpha_mat <- as.matrix(dat[, h_alpha_vars, drop=FALSE])
      h_full <- cbind(1, h_alpha_mat)
      n <- n_val
      p <- ncol(design_mat)
      h_dim <- ncol(h_full)

      model_fn <- function(alpha) {
        eta <- as.vector(design_mat %*% alpha)
        1 / (1 + exp(-eta))
      }

      Phi_alpha <- function(alpha) {
        pi_hat <- model_fn(alpha)
        (r_vec / pi_hat - 1) * h_full
      }

      Gamma_fn <- function(alpha) {
        pi_hat <- model_fn(alpha)
        common <- -r_vec * (1 - pi_hat) / pi_hat
        crossprod(h_full * common, design_mat) / n
      }

      W_func <- function(g_mat) solve(crossprod(g_mat) / nrow(g_mat))

      result <- inner_solver(rep(0, p), p, h_dim, model_fn, Phi_alpha, Gamma_fn, W_func,
                             design_mat, r_vec, y_vec, n)

      estimates <- result$estimates
      pi_hat <- model_fn(estimates)
      mu_hat <- mean(r_vec * y_vec / pi_hat)

      g_mat <- Phi_alpha(estimates)
      W_hat <- tryCatch(W_func(g_mat), error = function(e) diag(h_dim))
      Gamma_hat <- Gamma_fn(estimates)
      H_mat <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
      H_inv <- tryCatch(solve(H_mat), error = function(e) NULL)

      if (!is.null(H_inv)) {
        psi1 <- H_inv %*% crossprod(Gamma_hat, W_hat) %*% t(g_mat)
        dmu_dalpha <- colMeans(-r_vec * y_vec * (1 - pi_hat) / pi_hat * design_mat)
        iid1 <- r_vec * y_vec / pi_hat - mu_hat - as.vector(t(dmu_dalpha) %*% psi1)
        se1_vals[rep_i] <- sqrt(mean(iid1^2) / n)
      }

      mu_vals[rep_i] <- mu_hat
      alpha_mat[rep_i, ] <- estimates
      conv_vals[rep_i] <- result$converged
      grad_vals[rep_i] <- result$grad
    }, error = function(e) NULL)
  }

  valid <- !is.na(mu_vals) & !is.na(se1_vals)
  conv <- conv_vals[valid]
  esd <- sd(mu_vals[valid])
  ese1 <- mean(se1_vals[valid])

  cat(sprintf("\n===== %s =====\n", method_name))
  cat(sprintf("Valid: %d, Converged: %d\n", sum(valid), sum(conv)))
  cat(sprintf("Bias=%.4f, ESD=%.4f, ESE1=%.4f, ESE1/ESD=%.3f\n",
      mean(mu_vals[valid]) - mu_true, esd, ese1, ese1/esd))
  cat(sprintf("alpha[2]: mean=%.3f, sd=%.3f, range=[%.3f, %.3f]\n",
      mean(alpha_mat[valid, 2]), sd(alpha_mat[valid, 2]),
      min(alpha_mat[valid, 2]), max(alpha_mat[valid, 2])))
}

# ============================================================
# Method 1: Ridge penalty THROUGHOUT (not just step 1)
# ============================================================
solver_ridge_all <- function(init, p, h_dim, model_fn, Phi_alpha, Gamma_fn, W_func,
                             design_mat, r_vec, y_vec, n) {
  lambda <- 0.01
  obj_I <- function(alpha) {
    G <- colMeans(Phi_alpha(alpha)); sum(G^2) + lambda * sum(alpha^2)
  }
  opt0 <- optim(init, obj_I, method = "L-BFGS-B", control = list(maxit = 5000))
  estimates <- opt0$par

  best_grad <- Inf; best_est <- estimates; no_improve <- 0L; converged <- FALSE
  for (t in 1:200) {
    g_mat <- Phi_alpha(estimates)
    G_vec <- colMeans(g_mat)
    W_hat <- tryCatch(W_func(g_mat), error = function(e) diag(h_dim))
    Gamma_hat <- Gamma_fn(estimates)
    conv_grad <- 2 * as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))
    grad_norm <- max(abs(conv_grad))
    if (grad_norm < 1e-6) { converged <- TRUE; break }
    if (grad_norm < best_grad) { best_grad <- grad_norm; best_est <- estimates; no_improve <- 0L
    } else { no_improve <- no_improve + 1L; if (no_improve >= 20L) { estimates <- best_est; break } }

    obj_fn <- function(alpha) {
      G <- colMeans(Phi_alpha(alpha))
      as.numeric(t(G) %*% W_hat %*% G) + lambda * sum(alpha^2)
    }
    opt <- optim(estimates, obj_fn, method = "L-BFGS-B", control = list(maxit = 1000))
    estimates <- opt$par
  }
  list(estimates = estimates, converged = converged, grad = best_grad)
}

# ============================================================
# Method 2: One NR step per W update (the best known method)
# ============================================================
solver_one_nr <- function(init, p, h_dim, model_fn, Phi_alpha, Gamma_fn, W_func,
                          design_mat, r_vec, y_vec, n) {
  estimates <- init
  best_grad <- Inf; best_est <- estimates
  best_local_grad <- Inf; best_local_est <- estimates
  no_improve <- 0L; converged <- FALSE

  for (t in 1:1000) {
    g_mat <- Phi_alpha(estimates)
    G_vec <- colMeans(g_mat)
    W_hat <- tryCatch(W_func(g_mat), error = function(e) diag(h_dim))
    Gamma_hat <- Gamma_fn(estimates)
    g_vec <- 2 * as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))
    grad_norm <- max(abs(g_vec))

    d <- sqrt(sum((estimates - init)^2))
    if (d < 5 && grad_norm < best_local_grad) {
      best_local_grad <- grad_norm; best_local_est <- estimates
    }

    if (grad_norm < 1e-6) { converged <- TRUE; break }
    if (grad_norm < best_grad) { best_grad <- grad_norm; best_est <- estimates; no_improve <- 0L
    } else { no_improve <- no_improve + 1L
      if (no_improve >= 20L) {
        estimates <- if (best_local_grad < 0.1) best_local_est else best_est; break
      }
    }

    H_gn <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
    direction <- tryCatch(solve(H_gn, -g_vec), error = function(e) -g_vec)
    step_norm <- sqrt(sum(direction^2))
    if (step_norm > 0.5) direction <- direction * (0.5 / step_norm)

    obj_cur <- as.numeric(t(G_vec) %*% W_hat %*% G_vec)
    step <- 1.0
    for (ls in 1:20) {
      candidate <- estimates + step * direction
      G_new <- colMeans(Phi_alpha(candidate))
      obj_new <- as.numeric(t(G_new) %*% W_hat %*% G_new)
      if (is.finite(obj_new) && obj_new < obj_cur - 1e-4 * step * sum(g_vec * direction)) break
      step <- step * 0.5
    }
    estimates <- estimates + step * direction
  }

  d_final <- sqrt(sum((estimates - init)^2))
  if (d_final > 5 && best_local_grad < 0.1) estimates <- best_local_est
  list(estimates = estimates, converged = converged, grad = best_grad)
}

# ============================================================
# Method 3: nlminb with strict step control (rel.tol, step.min)
# ============================================================
solver_nlminb <- function(init, p, h_dim, model_fn, Phi_alpha, Gamma_fn, W_func,
                          design_mat, r_vec, y_vec, n) {
  obj_I <- function(alpha) { G <- colMeans(Phi_alpha(alpha)); sum(G^2) }
  opt0 <- nlminb(init, obj_I, control = list(iter.max = 5000, step.min = 0.01))
  estimates <- opt0$par

  best_grad <- Inf; best_est <- estimates; no_improve <- 0L; converged <- FALSE
  for (t in 1:200) {
    g_mat <- Phi_alpha(estimates)
    G_vec <- colMeans(g_mat)
    W_hat <- tryCatch(W_func(g_mat), error = function(e) diag(h_dim))
    Gamma_hat <- Gamma_fn(estimates)
    conv_grad <- 2 * as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))
    grad_norm <- max(abs(conv_grad))
    if (grad_norm < 1e-6) { converged <- TRUE; break }
    if (grad_norm < best_grad) { best_grad <- grad_norm; best_est <- estimates; no_improve <- 0L
    } else { no_improve <- no_improve + 1L; if (no_improve >= 20L) { estimates <- best_est; break } }

    obj_fn <- function(alpha) {
      G <- colMeans(Phi_alpha(alpha))
      as.numeric(t(G) %*% W_hat %*% G)
    }
    opt <- nlminb(estimates, obj_fn, control = list(iter.max = 1000, step.min = 0.01))
    estimates <- opt$par
  }
  list(estimates = estimates, converged = converged, grad = best_grad)
}

# ============================================================
# Method 4: Nelder-Mead (derivative-free, less aggressive steps)
# ============================================================
solver_nm <- function(init, p, h_dim, model_fn, Phi_alpha, Gamma_fn, W_func,
                      design_mat, r_vec, y_vec, n) {
  obj_I <- function(alpha) { G <- colMeans(Phi_alpha(alpha)); sum(G^2) }
  opt0 <- optim(init, obj_I, method = "Nelder-Mead", control = list(maxit = 10000))
  estimates <- opt0$par

  best_grad <- Inf; best_est <- estimates; no_improve <- 0L; converged <- FALSE
  for (t in 1:200) {
    g_mat <- Phi_alpha(estimates)
    G_vec <- colMeans(g_mat)
    W_hat <- tryCatch(W_func(g_mat), error = function(e) diag(h_dim))
    Gamma_hat <- Gamma_fn(estimates)
    conv_grad <- 2 * as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))
    grad_norm <- max(abs(conv_grad))
    if (grad_norm < 1e-6) { converged <- TRUE; break }
    if (grad_norm < best_grad) { best_grad <- grad_norm; best_est <- estimates; no_improve <- 0L
    } else { no_improve <- no_improve + 1L; if (no_improve >= 20L) { estimates <- best_est; break } }

    obj_fn <- function(alpha) {
      G <- colMeans(Phi_alpha(alpha))
      as.numeric(t(G) %*% W_hat %*% G)
    }
    opt <- optim(estimates, obj_fn, method = "Nelder-Mead", control = list(maxit = 10000))
    estimates <- opt$par
  }
  list(estimates = estimates, converged = converged, grad = best_grad)
}

# ============================================================
# Method 5: Multi-start Step 1 (pick smallest ||alpha|| among good solutions)
# then standard iterative GMM
# ============================================================
solver_multistart <- function(init, p, h_dim, model_fn, Phi_alpha, Gamma_fn, W_func,
                              design_mat, r_vec, y_vec, n) {
  obj_I <- function(alpha) { G <- colMeans(Phi_alpha(alpha)); sum(G^2) }

  # Run from zero + 9 random starts
  candidates <- list()
  starts <- c(list(init), lapply(1:9, function(s) rnorm(p, sd = 0.5)))
  for (s in seq_along(starts)) {
    opt <- optim(starts[[s]], obj_I, method = "L-BFGS-B", control = list(maxit = 5000))
    candidates[[s]] <- list(par = opt$par, value = opt$value, norm = sqrt(sum(opt$par^2)))
  }

  # Among solutions with obj < 2 * best_obj, pick smallest norm
  best_obj <- min(sapply(candidates, `[[`, "value"))
  good <- which(sapply(candidates, `[[`, "value") < best_obj * 2 + 1e-6)
  norms <- sapply(candidates[good], `[[`, "norm")
  best_idx <- good[which.min(norms)]
  estimates <- candidates[[best_idx]]$par

  # Standard iterative GMM
  best_grad <- Inf; best_est <- estimates; no_improve <- 0L; converged <- FALSE
  for (t in 1:200) {
    g_mat <- Phi_alpha(estimates)
    G_vec <- colMeans(g_mat)
    W_hat <- tryCatch(W_func(g_mat), error = function(e) diag(h_dim))
    Gamma_hat <- Gamma_fn(estimates)
    conv_grad <- 2 * as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))
    grad_norm <- max(abs(conv_grad))
    if (grad_norm < 1e-6) { converged <- TRUE; break }
    if (grad_norm < best_grad) { best_grad <- grad_norm; best_est <- estimates; no_improve <- 0L
    } else { no_improve <- no_improve + 1L; if (no_improve >= 20L) { estimates <- best_est; break } }

    obj_fn <- function(alpha) {
      G <- colMeans(Phi_alpha(alpha))
      as.numeric(t(G) %*% W_hat %*% G)
    }
    opt <- optim(estimates, obj_fn, method = "L-BFGS-B", control = list(maxit = 1000))
    estimates <- opt$par
  }
  list(estimates = estimates, converged = converged, grad = best_grad)
}

# ============================================================
# Method 6: Ridge penalty throughout with LARGER lambda
# ============================================================
solver_ridge_strong <- function(init, p, h_dim, model_fn, Phi_alpha, Gamma_fn, W_func,
                                design_mat, r_vec, y_vec, n) {
  lambda <- 0.1
  obj_I <- function(alpha) {
    G <- colMeans(Phi_alpha(alpha)); sum(G^2) + lambda * sum(alpha^2)
  }
  opt0 <- optim(init, obj_I, method = "L-BFGS-B", control = list(maxit = 5000))
  estimates <- opt0$par

  best_grad <- Inf; best_est <- estimates; no_improve <- 0L; converged <- FALSE
  for (t in 1:200) {
    g_mat <- Phi_alpha(estimates)
    G_vec <- colMeans(g_mat)
    W_hat <- tryCatch(W_func(g_mat), error = function(e) diag(h_dim))
    Gamma_hat <- Gamma_fn(estimates)
    conv_grad <- 2 * as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))
    grad_norm <- max(abs(conv_grad))
    if (grad_norm < 1e-6) { converged <- TRUE; break }
    if (grad_norm < best_grad) { best_grad <- grad_norm; best_est <- estimates; no_improve <- 0L
    } else { no_improve <- no_improve + 1L; if (no_improve >= 20L) { estimates <- best_est; break } }

    obj_fn <- function(alpha) {
      G <- colMeans(Phi_alpha(alpha))
      as.numeric(t(G) %*% W_hat %*% G) + lambda * sum(alpha^2)
    }
    opt <- optim(estimates, obj_fn, method = "L-BFGS-B", control = list(maxit = 1000))
    estimates <- opt$par
  }
  list(estimates = estimates, converged = converged, grad = best_grad)
}

# ============================================================
# Method 7: One NR step WITHOUT local preference (pure one-step)
# ============================================================
solver_one_nr_pure <- function(init, p, h_dim, model_fn, Phi_alpha, Gamma_fn, W_func,
                               design_mat, r_vec, y_vec, n) {
  estimates <- init
  best_grad <- Inf; best_est <- estimates; no_improve <- 0L; converged <- FALSE

  for (t in 1:1000) {
    g_mat <- Phi_alpha(estimates)
    G_vec <- colMeans(g_mat)
    W_hat <- tryCatch(W_func(g_mat), error = function(e) diag(h_dim))
    Gamma_hat <- Gamma_fn(estimates)
    g_vec <- 2 * as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))
    grad_norm <- max(abs(g_vec))

    if (grad_norm < 1e-6) { converged <- TRUE; break }
    if (grad_norm < best_grad) { best_grad <- grad_norm; best_est <- estimates; no_improve <- 0L
    } else { no_improve <- no_improve + 1L; if (no_improve >= 20L) { estimates <- best_est; break } }

    H_gn <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
    direction <- tryCatch(solve(H_gn, -g_vec), error = function(e) -g_vec)
    step_norm <- sqrt(sum(direction^2))
    if (step_norm > 0.5) direction <- direction * (0.5 / step_norm)

    obj_cur <- as.numeric(t(G_vec) %*% W_hat %*% G_vec)
    step <- 1.0
    for (ls in 1:20) {
      candidate <- estimates + step * direction
      G_new <- colMeans(Phi_alpha(candidate))
      obj_new <- as.numeric(t(G_new) %*% W_hat %*% G_new)
      if (is.finite(obj_new) && obj_new < obj_cur - 1e-4 * step * sum(g_vec * direction)) break
      step <- step * 0.5
    }
    estimates <- estimates + step * direction
  }
  list(estimates = estimates, converged = converged, grad = best_grad)
}

# Run all methods
run_method("1. Ridge lambda=0.01 throughout", solver_ridge_all)
run_method("2. One NR step + trust region + local pref", solver_one_nr)
run_method("3. nlminb with step.min=0.01", solver_nlminb)
run_method("4. Nelder-Mead", solver_nm)
run_method("5. Multi-start (10 starts, pick min norm)", solver_multistart)
run_method("6. Ridge lambda=0.1 throughout", solver_ridge_strong)
run_method("7. One NR step + trust region (no local pref)", solver_one_nr_pure)
