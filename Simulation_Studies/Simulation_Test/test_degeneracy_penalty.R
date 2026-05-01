setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})

n_val <- 2000
n_reps <- 200

# Idea: add a penalty to the GMM objective that becomes large when
# the Jacobian Gamma degenerates (i.e., when we enter the flat plateau).
#
# The penalty: lambda * max(0, log(cond(H)) - log(cond_ok))^2
# This is zero when cond(H) < cond_ok, and grows when cond(H) exceeds it.
#
# But computing cond(H) at every function evaluation is expensive.
# Cheaper proxy: ||1 - pi||^2 should be bounded away from 0 at good solutions.
# When pi -> 1, (1-pi) -> 0, which is exactly the degeneracy condition.
#
# Even simpler: min(1 - pi_hat) is a direct measure.
# Or: mean(pi*(1-pi)) — this is proportional to the "information" in alpha.
# At good solutions: mean(pi*(1-pi)) ~ 0.25 (balanced).
# At degenerate solutions: mean(pi*(1-pi)) -> 0 (pi -> 0 or 1).
#
# Penalty: -lambda * log(mean(pi*(1-pi)))
# This is small when pi is balanced, large when pi -> 0 or 1.

run_test <- function(setting, miss, ps_spec_id, model_idx, lambda) {
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

      # Degeneracy penalty: penalize when pi*(1-pi) is too small
      # This prevents the optimizer from entering the flat plateau
      deg_penalty <- function(alpha) {
        pi_hat <- model_fn(alpha)
        info <- mean(pi_hat * (1 - pi_hat))
        if (info < 1e-10) return(1e10)
        -lambda * log(info)
      }

      # === Step 1: W=I with degeneracy penalty ===
      obj_I <- function(alpha) {
        G <- colMeans(Phi_alpha(alpha))
        sum(G^2) + deg_penalty(alpha)
      }
      opt0 <- optim(rep(0, p), obj_I, method = "L-BFGS-B", control = list(maxit = 5000))
      estimates <- opt0$par

      # === Step 2+: standard iterated GMM with degeneracy penalty in inner loop ===
      best_grad <- Inf
      best_est <- estimates
      no_improve <- 0L
      converged <- FALSE

      for (t in 1:200) {
        g_mat <- Phi_alpha(estimates)
        G_vec <- colMeans(g_mat)
        W_hat <- tryCatch(W_func(g_mat), error = function(e) diag(h_dim))

        # Convergence check on UNPENALIZED gradient
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

        # Full inner minimization with penalty
        obj_fn <- function(alpha) {
          G <- colMeans(Phi_alpha(alpha))
          as.numeric(t(G) %*% W_hat %*% G) + deg_penalty(alpha)
        }
        opt <- optim(estimates, obj_fn, method = "L-BFGS-B", control = list(maxit = 1000))
        estimates <- opt$par
      }

      pi_hat <- model_fn(estimates)
      mu_hat <- mean(r_vec * y_vec / pi_hat)

      # SE1 (based on unpenalized GMM)
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

  cat(sprintf("\n=== %s model%d, lambda=%.4f ===\n", setting, model_idx, lambda))
  cat(sprintf("Valid: %d, Converged: %d\n", sum(valid), sum(conv)))
  cat(sprintf("Bias=%.4f, ESD=%.4f\n", mean(mu_vals[valid]) - mu_true, esd))
  cat(sprintf("ESE1=%.4f, ESE1/ESD=%.3f\n", ese1, ese1/esd))
  cat(sprintf("alpha[2]: mean=%.3f, sd=%.3f, range=[%.3f, %.3f]\n",
      mean(alpha2_vals[valid]), sd(alpha2_vals[valid]),
      min(alpha2_vals[valid]), max(alpha2_vals[valid])))
  cat(sprintf("mean(pi): "))
  pi_means <- sapply(which(valid), function(i) {
    dat <- all_data[((i-1)*n_val + 1):(i*n_val), ]
    dm <- model.matrix(formula_j, data=dat)
    mean(1/(1+exp(-as.vector(dm %*% c(mu_vals[i], alpha2_vals[i])))))
  })
  # Just report alpha range as proxy
  cat(sprintf("(see alpha range)\n"))
}

# Setting4 Model3: test different lambda values
for (lam in c(0.001, 0.01, 0.1, 1.0)) {
  run_test("setting4", "miss50", "9-alt1", 3, lambda = lam)
}

# Setting3 Model2 with best lambda
cat("\n\n========== Setting3 Model2 ==========\n")
run_test("setting3", "miss50", "9-alt1", 2, lambda = 0.01)
