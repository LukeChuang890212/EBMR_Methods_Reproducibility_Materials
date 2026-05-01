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
  guard_count <- 0L

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

      # Compute cond(Gamma' W Gamma) at alpha, using given or computed W
      get_condH <- function(alpha, W_hat = NULL) {
        if (is.null(W_hat)) {
          g_mat <- Phi_alpha(alpha)
          W_hat <- tryCatch(W_func(g_mat), error = function(e) diag(h_dim))
        }
        Gamma_hat <- Gamma_fn(alpha)
        H <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
        eig <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
        max(eig) / max(min(eig), 1e-15)
      }

      # === Step 1: W=I ===
      # Use NR with step rejection based on cond(H)
      # Start from zero, take steps but reject any that make cond(H) too large
      alpha <- rep(0, p)
      W_hat <- diag(h_dim)

      for (k in 1:200) {
        g_mat <- Phi_alpha(alpha)
        G_vec <- colMeans(g_mat)
        Gamma_hat <- Gamma_fn(alpha)
        g_vec <- 2 * as.vector(crossprod(Gamma_hat, G_vec))  # W=I
        inner_grad <- max(abs(g_vec))
        if (inner_grad < 1e-8) break

        H_gn <- crossprod(Gamma_hat, Gamma_hat)  # W=I
        direction <- tryCatch(solve(H_gn, -g_vec), error = function(e) -g_vec)

        # Trust region
        step_norm <- sqrt(sum(direction^2))
        if (step_norm > 0.5) direction <- direction * (0.5 / step_norm)

        # Line search with degeneracy check
        obj_cur <- sum(G_vec^2)
        step <- 1.0
        accepted <- FALSE
        for (ls in 1:20) {
          candidate <- alpha + step * direction
          G_new <- colMeans(Phi_alpha(candidate))
          obj_new <- sum(G_new^2)

          # Check Armijo AND non-degeneracy
          if (is.finite(obj_new) && obj_new < obj_cur - 1e-4 * step * sum(g_vec * direction)) {
            # Check cond(H) at candidate
            cond_cand <- get_condH(candidate, diag(h_dim))
            if (cond_cand < cond_threshold) {
              accepted <- TRUE
              break
            }
          }
          step <- step * 0.5
        }

        if (accepted) {
          alpha <- candidate
        } else {
          guard_count <- guard_count + 1L
          break  # can't make progress without entering degenerate region
        }
      }
      estimates <- alpha

      # === Step 2+: iterative GMM ===
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

        # Inner: NR with cond(H) guard
        alpha_inner <- estimates
        for (k in 1:200) {
          g_mat_k <- Phi_alpha(alpha_inner)
          G_k <- colMeans(g_mat_k)
          Gamma_k <- Gamma_fn(alpha_inner)
          g_k <- 2 * as.vector(crossprod(Gamma_k, W_hat %*% G_k))
          ig <- max(abs(g_k))
          if (ig < 1e-8) break

          H_k <- crossprod(Gamma_k, W_hat %*% Gamma_k)
          dir_k <- tryCatch(solve(H_k, -g_k), error = function(e) -g_k)
          sn <- sqrt(sum(dir_k^2))
          if (sn > 0.5) dir_k <- dir_k * (0.5 / sn)

          obj_k <- as.numeric(t(G_k) %*% W_hat %*% G_k)
          step <- 1.0
          accepted <- FALSE
          for (ls in 1:20) {
            cand <- alpha_inner + step * dir_k
            G_new <- colMeans(Phi_alpha(cand))
            obj_new <- as.numeric(t(G_new) %*% W_hat %*% G_new)
            if (is.finite(obj_new) && obj_new < obj_k - 1e-4 * step * sum(g_k * dir_k)) {
              cond_cand <- get_condH(cand, W_hat)
              if (cond_cand < cond_threshold) {
                accepted <- TRUE
                break
              }
            }
            step <- step * 0.5
          }

          if (accepted) {
            alpha_inner <- cand
          } else {
            guard_count <- guard_count + 1L
            break
          }
        }
        estimates <- alpha_inner
      }

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
  cat(sprintf("Valid: %d, Converged: %d, Guard triggered: %d\n", sum(valid), sum(conv), guard_count))
  cat(sprintf("Bias=%.4f, ESD=%.4f\n", mean(mu_vals[valid]) - mu_true, esd))
  cat(sprintf("ESE1=%.4f, ESE1/ESD=%.3f\n", ese1, ese1/esd))
  cat(sprintf("alpha[2]: mean=%.3f, sd=%.3f, range=[%.3f, %.3f]\n",
      mean(alpha2_vals[valid]), sd(alpha2_vals[valid]),
      min(alpha2_vals[valid]), max(alpha2_vals[valid])))
}

run_test("setting4", "miss50", "9-alt1", 3, cond_threshold = 1e6)
run_test("setting3", "miss50", "9-alt1", 2, cond_threshold = 1e6)
