setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})

n_val <- 2000
n_reps <- 200
ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting3"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
mu_true <- get_mu_true("setting3")

# Model 2: formula = r ~ y + u1 + z2, h_alpha = (u1, u2, z1, z2)
formula_j <- ps_spec[["formula.list"]][[2]]
h_alpha_vars <- ps_spec[["h_alpha.list"]][[2]]

cat(sprintf("Formula: %s\n", deparse(formula_j)))
cat(sprintf("h_alpha: %s\n", paste(h_alpha_vars, collapse=", ")))

run_method <- function(method_name, use_one_nr) {
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

      if (use_one_nr) {
        # === One NR step + trust region + local preference ===
        estimates <- rep(0, p)
        init <- estimates
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

      } else {
        # === Baseline: L-BFGS-B iterative GMM ===
        obj_I <- function(alpha) { G <- colMeans(Phi_alpha(alpha)); sum(G^2) }
        opt0 <- optim(rep(0, p), obj_I, method = "L-BFGS-B", control = list(maxit = 5000))
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
          opt <- optim(estimates, obj_fn, method = "L-BFGS-B", control = list(maxit = 1000))
          estimates <- opt$par
        }
      }

      # === Compute mu, SE1, SE2 ===
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

        # SE2: CUE correction
        G_vec <- colMeans(g_mat)
        WG <- W_hat %*% G_vec
        gi_dot_WG <- as.vector(g_mat %*% WG)
        Omega <- crossprod(g_mat) / n
        OmegaWG <- as.vector(Omega %*% WG)
        GtW <- crossprod(Gamma_hat, W_hat)
        GtW_g <- GtW %*% t(g_mat)
        GtW_OmegaWG <- as.vector(GtW %*% OmegaWG)
        correction <- sweep(GtW_g * rep(gi_dot_WG, each = p), 1, GtW_OmegaWG, "-")
        psi2 <- psi1 - H_inv %*% correction

        iid2 <- r_vec * y_vec / pi_hat - mu_hat - as.vector(t(dmu_dalpha) %*% psi2)
        se2_vals[rep_i] <- sqrt(mean(iid2^2) / n)
      }

      mu_vals[rep_i] <- mu_hat
      alpha_mat[rep_i, ] <- estimates
      conv_vals[rep_i] <- converged
      grad_vals[rep_i] <- if (use_one_nr) best_grad else best_grad
    }, error = function(e) NULL)
  }

  valid <- !is.na(mu_vals) & !is.na(se1_vals)
  conv <- conv_vals[valid]
  esd <- sd(mu_vals[valid])
  ese1 <- mean(se1_vals[valid])
  ese2 <- mean(se2_vals[valid], na.rm=TRUE)

  cat(sprintf("\n===== %s =====\n", method_name))
  cat(sprintf("Valid: %d, Converged: %d\n", sum(valid), sum(conv)))
  cat(sprintf("Bias=%.4f, ESD=%.4f\n", mean(mu_vals[valid]) - mu_true, esd))
  cat(sprintf("ESE1=%.4f, ESE1/ESD=%.3f\n", ese1, ese1/esd))
  cat(sprintf("ESE2=%.4f, ESE2/ESD=%.3f\n", ese2, ese2/esd))
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
  cat(sprintf("mean(pi): mean=%.4f, range=[%.4f, %.4f]\n",
      mean(pi_means[valid]), min(pi_means[valid]), max(pi_means[valid])))
}

cat(sprintf("mu.true = %.4f\n\n", mu_true))

run_method("Baseline (L-BFGS-B iterative GMM)", use_one_nr = FALSE)
run_method("One NR step + trust region + local pref", use_one_nr = TRUE)
