setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})

n_val <- 2000
n_reps <- 200

# ============================================================
# Jacobian-guarded inner optimizer
# After L-BFGS-B converges, check if Gamma'W*Gamma is degenerate.
# If so, reject the solution and use the last non-degenerate iterate.
# ============================================================
run_jacobian_guard <- function(setting, miss, ps_spec_id, model_idx) {
  ps_spec <- get_ps_spec(ps_spec_id)
  data_file <- misspecified_model_all_data_file.list[[setting]][[miss]][[1]]
  all_data <- readRDS(data_file)
  mu_true <- get_mu_true(setting)

  formula_j <- ps_spec[["formula.list"]][[model_idx]]
  h_alpha_vars <- ps_spec[["h_alpha.list"]][[model_idx]]

  cat(sprintf("Setting: %s, PS spec: %s, Model: %d\n", setting, ps_spec_id, model_idx))
  cat(sprintf("Formula: %s\n", deparse(formula_j)))
  cat(sprintf("h_alpha: %s\n", paste(h_alpha_vars, collapse=", ")))
  cat(sprintf("mu.true = %.4f\n\n", mu_true))

  mu_vals <- rep(NA_real_, n_reps)
  se1_vals <- rep(NA_real_, n_reps)
  alpha_mat <- matrix(NA, n_reps, length(all.vars(formula_j)))
  conv_vals <- rep(NA, n_reps)
  grad_vals <- rep(NA_real_, n_reps)
  guard_triggered <- rep(FALSE, n_reps)

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

      # Compute Jacobian norm (Frobenius) as health indicator
      gamma_norm <- function(alpha) {
        Gamma_hat <- Gamma_fn(alpha)
        sqrt(sum(Gamma_hat^2))
      }

      # Jacobian-guarded inner optimizer:
      # Run L-BFGS-B, but after convergence check if Gamma is degenerate.
      # If degenerate, use the pre-optimization iterate instead.
      inner_opt_guarded <- function(start, W_hat, gamma_ref) {
        obj_fn <- function(alpha) {
          G <- colMeans(Phi_alpha(alpha))
          as.numeric(t(G) %*% W_hat %*% G)
        }

        opt <- optim(start, obj_fn, method = "L-BFGS-B", control = list(maxit = 1000))
        candidate <- opt$par

        # Check: is the Jacobian at the candidate degenerate?
        gamma_cand <- gamma_norm(candidate)

        # Guard: if Jacobian norm dropped below 10% of reference, reject
        if (gamma_cand < 0.1 * gamma_ref) {
          return(list(par = start, guarded = TRUE))
        }

        list(par = candidate, guarded = FALSE)
      }

      # === Step 1: W=I ===
      obj_I <- function(alpha) { G <- colMeans(Phi_alpha(alpha)); sum(G^2) }
      opt0 <- optim(rep(0, p), obj_I, method = "L-BFGS-B", control = list(maxit = 5000))
      estimates <- opt0$par

      # Compute reference Jacobian norm at step 1 result
      gamma_ref <- gamma_norm(estimates)
      # If step 1 itself is already degenerate, use zero init's Jacobian as ref
      gamma_zero <- gamma_norm(rep(0, p))
      if (gamma_ref < 0.1 * gamma_zero) {
        estimates <- rep(0, p)  # reject step 1
        gamma_ref <- gamma_zero
        guard_triggered[rep_i] <- TRUE
      }

      # === Step 2+: iterative GMM with Jacobian guard ===
      best_grad <- Inf
      best_est <- estimates
      no_improve <- 0L
      converged <- FALSE

      for (t in 1:200) {
        g_mat <- Phi_alpha(estimates)
        G_vec <- colMeans(g_mat)
        W_hat <- tryCatch(W_func(g_mat), error = function(e) diag(h_dim))

        Gamma_hat <- Gamma_fn(estimates)
        conv_grad <- 2 * as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))
        grad_norm <- max(abs(conv_grad))

        if (grad_norm < 1e-6) { converged <- TRUE; break }

        if (grad_norm < best_grad) {
          best_grad <- grad_norm
          best_est <- estimates
          no_improve <- 0L
        } else {
          no_improve <- no_improve + 1L
          if (no_improve >= 20L) { estimates <- best_est; break }
        }

        result <- inner_opt_guarded(estimates, W_hat, gamma_ref)
        if (result$guarded) {
          guard_triggered[rep_i] <- TRUE
          # Don't update estimates — stay at current iterate
        } else {
          estimates <- result$par
          # Update reference: use running max to avoid ratcheting down
          gamma_ref <- max(gamma_ref, gamma_norm(estimates))
        }
      }

      # Final check: is the solution degenerate?
      gamma_final <- gamma_norm(estimates)
      if (gamma_final < 0.1 * gamma_zero) {
        guard_triggered[rep_i] <- TRUE
        estimates <- best_est  # use best non-degenerate
      }

      pi_hat <- model_fn(estimates)
      mu_hat <- mean(r_vec * y_vec / pi_hat)

      # SE1
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
      conv_vals[rep_i] <- converged
      grad_vals[rep_i] <- best_grad
    }, error = function(e) NULL)
  }

  valid <- !is.na(mu_vals) & !is.na(se1_vals)
  conv <- conv_vals[valid]
  esd <- sd(mu_vals[valid])
  ese1 <- mean(se1_vals[valid])

  cat(sprintf("Valid: %d, Converged: %d, Guard triggered: %d\n", sum(valid), sum(conv), sum(guard_triggered)))
  cat(sprintf("Bias=%.4f, ESD=%.4f\n", mean(mu_vals[valid]) - mu_true, esd))
  cat(sprintf("ESE1=%.4f, ESE1/ESD=%.3f\n", ese1, ese1/esd))
  cat(sprintf("alpha[2]: mean=%.3f, sd=%.3f, range=[%.3f, %.3f]\n",
      mean(alpha_mat[valid, 2]), sd(alpha_mat[valid, 2]),
      min(alpha_mat[valid, 2]), max(alpha_mat[valid, 2])))

  pi_means <- rep(NA, n_reps)
  for (i in which(valid)) {
    dat <- all_data[((i-1)*n_val + 1):(i*n_val), ]
    dm <- model.matrix(formula_j, data=dat)
    eta <- as.vector(dm %*% alpha_mat[i,])
    pi_means[i] <- mean(1/(1+exp(-eta)))
  }
  cat(sprintf("mean(pi): mean=%.4f, range=[%.4f, %.4f]\n\n",
      mean(pi_means[valid]), min(pi_means[valid]), max(pi_means[valid])))
}

# Test on setting4 model3 (the known problem case)
cat("============================================================\n")
cat("Setting4, Model 3 (r ~ y + u2 + z2, h_alpha = full)\n")
cat("============================================================\n")
run_jacobian_guard("setting4", "miss50", "9-alt1", 3)

# Test on setting3 model2
cat("============================================================\n")
cat("Setting3, Model 2 (r ~ y + u1 + z2, h_alpha = full)\n")
cat("============================================================\n")
run_jacobian_guard("setting3", "miss50", "9-alt1", 2)
