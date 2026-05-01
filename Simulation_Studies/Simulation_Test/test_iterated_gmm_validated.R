setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})

n_val <- 2000
n_reps <- 200

run_test <- function(setting, miss, ps_spec_id, model_idx, cond_threshold = 1e6) {
  ps_spec <- get_ps_spec(ps_spec_id)
  data_file <- misspecified_model_all_data_file.list[[setting]][[miss]][[1]]
  all_data <- readRDS(data_file)
  mu_true <- get_mu_true(setting)

  formula_j <- ps_spec[["formula.list"]][[model_idx]]
  h_alpha_vars <- ps_spec[["h_alpha.list"]][[model_idx]]

  mu_vals <- rep(NA_real_, n_reps)
  se1_vals <- rep(NA_real_, n_reps)
  alpha2_vals <- rep(NA_real_, n_reps)
  conv_vals <- rep(NA, n_reps)
  grad_vals <- rep(NA_real_, n_reps)
  reject_count <- 0L

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

      # Check if solution is on a degenerate flat plateau
      is_degenerate <- function(alpha, W_hat) {
        Gamma_hat <- Gamma_fn(alpha)
        H <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
        eig <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
        cond_num <- max(eig) / max(min(eig), 1e-15)
        cond_num > cond_threshold
      }

      # === Step 1: W = I, full minimization ===
      obj_I <- function(alpha) {
        G <- colMeans(Phi_alpha(alpha))
        sum(G^2)
      }
      opt0 <- optim(rep(0, p), obj_I, method = "L-BFGS-B", control = list(maxit = 5000))
      estimates <- opt0$par

      # Validate step 1: if degenerate, reject and use zero init
      if (is_degenerate(estimates, diag(h_dim))) {
        reject_count <- reject_count + 1L
        estimates <- rep(0, p)
      }

      # === Step 2+: standard iterated GMM with full inner minimization ===
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

        # Full inner minimization
        obj_fn <- function(alpha) {
          G <- colMeans(Phi_alpha(alpha))
          as.numeric(t(G) %*% W_hat %*% G)
        }
        opt <- optim(estimates, obj_fn, method = "L-BFGS-B", control = list(maxit = 1000))
        candidate <- opt$par

        # Validate: reject if candidate landed on flat plateau
        if (is_degenerate(candidate, W_hat)) {
          reject_count <- reject_count + 1L
          # Don't update estimates — keep current iterate, W will update next round
        } else {
          estimates <- candidate
        }
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
      alpha2_vals[rep_i] <- estimates[2]
      conv_vals[rep_i] <- converged
      grad_vals[rep_i] <- best_grad
    }, error = function(e) NULL)
  }

  valid <- !is.na(mu_vals) & !is.na(se1_vals)
  conv <- conv_vals[valid]
  esd <- sd(mu_vals[valid])
  ese1 <- mean(se1_vals[valid])

  cat(sprintf("\n=== %s, model %d (cond_threshold=%.0e) ===\n", setting, model_idx, cond_threshold))
  cat(sprintf("Valid: %d, Converged: %d, Rejections: %d\n", sum(valid), sum(conv), reject_count))
  cat(sprintf("Bias=%.4f, ESD=%.4f\n", mean(mu_vals[valid]) - mu_true, esd))
  cat(sprintf("ESE1=%.4f, ESE1/ESD=%.3f\n", ese1, ese1/esd))
  cat(sprintf("alpha[2]: mean=%.3f, sd=%.3f, range=[%.3f, %.3f]\n",
      mean(alpha2_vals[valid]), sd(alpha2_vals[valid]),
      min(alpha2_vals[valid]), max(alpha2_vals[valid])))
}

# Test with different thresholds
for (ct in c(1e4, 1e5, 1e6, 1e8)) {
  run_test("setting4", "miss50", "9-alt1", 3, cond_threshold = ct)
}

# Also test setting3 model2 with best threshold
cat("\n\n========== Setting3 Model2 ==========\n")
run_test("setting3", "miss50", "9-alt1", 2, cond_threshold = 1e6)
