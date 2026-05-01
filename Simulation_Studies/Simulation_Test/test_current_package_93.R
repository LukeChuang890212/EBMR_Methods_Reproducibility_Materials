setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

n_val <- 2000
n_reps <- 200
ps_spec <- get_ps_spec("9-alt1")
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
compute_pi_fn <- function(eta) plogis(-eta)
cond_limit <- 1e8

run_sim <- function(setting, m_idx, label) {
  cat(sprintf("\n=== %s ===\n", label))
  data_file <- misspecified_model_all_data_file.list[[setting]][["miss50"]][[1]]
  all_data <- readRDS(data_file)
  mu_true <- get_mu_true(setting)

  mu_vals <- se2_vals <- cond_vals <- grad_vals <- rep(NA_real_, n_reps)

  for (rep_i in 1:n_reps) {
    if (rep_i %% 50 == 0) cat(sprintf("  rep %d...\n", rep_i))
    dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
    single_ps <- list(
      formula.list = list(ps_spec[["formula.list"]][[m_idx]]),
      h_alpha.list = list(ps_spec[["h_alpha.list"]][[m_idx]]),
      inv_link = ps_spec[["inv_link"]],
      outcome = ps_spec[["outcome"]],
      alpha_init.list = list(NULL),
      optimizer = "constrained_nr"
    )
    tryCatch({
      ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
      gmm_fit <- ebmr$ps_fit.list[[1]]$gmm_fit
      cond_vals[rep_i] <- gmm_fit$opt$solution_cond
      grad_vals[rep_i] <- gmm_fit$opt$final_grad_norm

      pi_hat <- ebmr$ps_fit.list[[1]]$fitted.values
      design_mat <- ebmr$ps_fit.list[[1]]$design_matrix
      link_type <- ebmr$ps_fit.list[[1]]$link_type
      r_vec <- dat[["r"]]; y_vec <- dat[["y"]]

      psi2_mat <- gmm_fit$psi2
      if (!is.null(psi2_mat) && is.matrix(psi2_mat) && !any(is.na(psi2_mat))) {
        if (!is.null(link_type) && link_type == "logistic_complement") {
          dot_pi <- -design_mat * (pi_hat * (1 - pi_hat))
        } else {
          dot_pi <- design_mat * (pi_hat * (1 - pi_hat))
        }
        ry_ps_inv2 <- as.vector(r_vec * y_vec * (pi_hat^(-2)))
        H_alpha <- colMeans(dot_pi * ry_ps_inv2)
        mu_iid <- as.vector(t(r_vec/pi_hat*y_vec) - t(H_alpha) %*% psi2_mat)
        mu_vals[rep_i] <- mean(r_vec * y_vec / pi_hat)
        se2_vals[rep_i] <- sqrt(var(mu_iid) / n_val)
      }
    }, error = function(e) {
      cat(sprintf("  Rep %d ERROR: %s\n", rep_i, e$message))
    })
  }

  valid <- !is.na(mu_vals) & !is.na(se2_vals)
  n_degen <- sum(cond_vals[valid] >= cond_limit, na.rm = TRUE)
  n_conv <- sum(grad_vals[valid] < 1e-6, na.rm = TRUE)
  cat(sprintf("  Valid: %d/%d, Degenerate: %d, Converged: %d\n",
      sum(valid), n_reps, n_degen, n_conv))
  cat(sprintf("  Bias     = %.4f\n", mean(mu_vals[valid]) - mu_true))
  cat(sprintf("  ESD      = %.4f\n", sd(mu_vals[valid])))
  cat(sprintf("  ESE2     = %.4f\n", mean(se2_vals[valid])))
  cat(sprintf("  ESE2/ESD = %.3f\n", mean(se2_vals[valid]) / sd(mu_vals[valid])))

  # Coverage
  ci_lo <- mu_vals[valid] - 1.96 * se2_vals[valid]
  ci_hi <- mu_vals[valid] + 1.96 * se2_vals[valid]
  cp <- mean(mu_true >= ci_lo & mu_true <= ci_hi)
  cat(sprintf("  CP       = %.3f\n", cp))
}

run_sim("setting4", 3, "Setting 4, M3 (r ~ y + u2 + z2), n=2000, miss50, scenario 9-3")
run_sim("setting4", 2, "Setting 4, M2 (r ~ y + u1 + z2), n=2000, miss50, scenario 9-3")
run_sim("setting3", 3, "Setting 3, M3 (r ~ y + u2 + z2), n=2000, miss50, scenario 9-3")
