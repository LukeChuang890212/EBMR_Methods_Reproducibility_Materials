setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})

library(EBMRalgorithmFast4)

n_val <- 2000
n_reps <- 200

run_package_single_model <- function(setting, miss, ps_spec_id, model_idx) {
  ps_spec <- get_ps_spec(ps_spec_id)
  data_file <- misspecified_model_all_data_file.list[[setting]][[miss]][[1]]
  all_data <- readRDS(data_file)
  mu_true <- get_mu_true(setting)

  # Build single-model ps_specifications in the format the package expects
  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL)
  )

  W <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

  mu_vals <- rep(NA_real_, n_reps)
  se1_vals <- rep(NA_real_, n_reps)
  se2_vals <- rep(NA_real_, n_reps)
  alpha2_vals <- rep(NA_real_, n_reps)

  for (rep_i in 1:n_reps) {
    dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
    tryCatch({
      ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W)

      r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
      pi_hat <- ebmr$ps_fit.list[[1]]$fitted.values
      mu_hat <- mean(r_vec * y_vec / pi_hat)
      n <- n_val

      gmm_fit <- ebmr$ps_fit.list[[1]]$gmm_fit
      alpha <- gmm_fit$estimates
      psi1 <- gmm_fit$psi   # k x n
      psi2_mat <- gmm_fit$psi2  # k x n (may be NULL or have NAs)

      # SE using psi1
      design_mat <- ebmr$ps_fit.list[[1]]$design_matrix
      link_type <- ebmr$ps_fit.list[[1]]$link_type
      if (!is.null(link_type) && link_type == "logistic_complement") {
        dot_pi <- -design_mat * (pi_hat * (1 - pi_hat))
      } else {
        dot_pi <- design_mat * (pi_hat * (1 - pi_hat))
      }
      ry_ps_inv2 <- as.vector(r_vec * y_vec * (pi_hat^(-2)))
      H_alpha <- colMeans(dot_pi * ry_ps_inv2)

      mu_iid1 <- as.vector(t(r_vec/pi_hat*y_vec) - t(H_alpha) %*% psi1)
      se1_vals[rep_i] <- sqrt(var(mu_iid1)/n)

      # SE using psi2
      if (!is.null(psi2_mat) && is.matrix(psi2_mat) && !any(is.na(psi2_mat))) {
        mu_iid2 <- as.vector(t(r_vec/pi_hat*y_vec) - t(H_alpha) %*% psi2_mat)
        se2_vals[rep_i] <- sqrt(var(mu_iid2)/n)
      }

      mu_vals[rep_i] <- mu_hat
      alpha2_vals[rep_i] <- alpha[2]
    }, error = function(e) {
      if (rep_i <= 3) cat(sprintf("  Rep %d error: %s\n", rep_i, conditionMessage(e)))
    })
  }

  valid1 <- !is.na(mu_vals) & !is.na(se1_vals)
  valid2 <- !is.na(mu_vals) & !is.na(se2_vals)

  cat(sprintf("\n=== Package baseline: %s, model %d ===\n", setting, model_idx))
  cat(sprintf("Valid (SE1): %d, Valid (SE2): %d\n", sum(valid1), sum(valid2)))
  if (sum(valid1) > 5) {
    esd <- sd(mu_vals[valid1]); ese1 <- mean(se1_vals[valid1])
    cat(sprintf("Bias=%.4f, ESD=%.4f\n", mean(mu_vals[valid1]) - mu_true, esd))
    cat(sprintf("ESE1=%.4f, ESE1/ESD=%.3f\n", ese1, ese1/esd))
  }
  if (sum(valid2) > 5) {
    esd2 <- sd(mu_vals[valid2]); ese2 <- mean(se2_vals[valid2])
    cat(sprintf("ESE2=%.4f, ESE2/ESD=%.3f\n", ese2, ese2/esd2))
  } else {
    cat(sprintf("SE2: only %d valid (not enough)\n", sum(valid2)))
  }
  if (any(!is.na(alpha2_vals))) {
    v <- !is.na(alpha2_vals)
    cat(sprintf("alpha[2]: mean=%.3f, sd=%.3f, range=[%.3f, %.3f]\n",
        mean(alpha2_vals[v]), sd(alpha2_vals[v]),
        min(alpha2_vals[v]), max(alpha2_vals[v])))
  }
}

# Setting4 model 3 (the case we fixed)
run_package_single_model("setting4", "miss50", "9-alt1", 3)

# Setting4 model 2
run_package_single_model("setting4", "miss50", "9-alt1", 2)

# Setting3 model 2
run_package_single_model("setting3", "miss50", "9-alt1", 2)
