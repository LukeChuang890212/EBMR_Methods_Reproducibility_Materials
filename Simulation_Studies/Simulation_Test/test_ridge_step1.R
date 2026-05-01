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

# Try different lambda values for ridge penalty at Step 1
for (lambda in c(0.001, 0.01, 0.1)) {
  cat(sprintf("\n============ lambda = %.3f ============\n", lambda))

  mu_vals <- rep(NA_real_, n_reps)
  se1_vals <- rep(NA_real_, n_reps)
  se2_vals <- rep(NA_real_, n_reps)
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

      # Step 1: W=I with ridge penalty: minimize G'G + lambda * ||alpha||^2
      obj_ridge <- function(alpha) {
        G <- colMeans(Phi_alpha(alpha))
        sum(G^2) + lambda * sum(alpha^2)
      }
      opt0 <- optim(rep(0, p), obj_ridge, method = "L-BFGS-B",
                    control = list(maxit = 5000))
      estimates <- opt0$par

      # Step 2+: standard iterative GMM, NO penalty, NO bounds
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

        obj_fn <- function(alpha) {
          G <- colMeans(Phi_alpha(alpha))
          as.numeric(t(G) %*% W_hat %*% G)
        }
        opt <- optim(estimates, obj_fn, method = "L-BFGS-B",
                     control = list(maxit = 1000))
        estimates <- opt$par
      }

      # Compute mu and SE1
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
      conv_vals[rep_i] <- converged
      grad_vals[rep_i] <- best_grad
    }, error = function(e) NULL)
  }

  valid <- !is.na(mu_vals) & !is.na(se1_vals)
  conv <- conv_vals[valid]

  esd <- sd(mu_vals[valid])
  ese1 <- mean(se1_vals[valid])

  cat(sprintf("Valid: %d, Converged: %d, Not converged: %d\n", sum(valid), sum(conv), sum(!conv)))
  cat(sprintf("Bias=%.4f, ESD=%.4f\n", mean(mu_vals[valid]) - mu_true, esd))
  cat(sprintf("ESE1=%.4f, ESE1/ESD=%.3f\n", ese1, ese1/esd))
  cat(sprintf("alpha[2]: mean=%.4f, sd=%.4f, range=[%.4f, %.4f]\n",
      mean(alpha_mat[valid, 2]), sd(alpha_mat[valid, 2]),
      min(alpha_mat[valid, 2]), max(alpha_mat[valid, 2])))
}
