setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})

# Reload the modified package
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

n_val <- 2000
n_reps <- 200

run_test <- function(setting, miss, ps_spec_id, model_idx, optimizer_name) {
  ps_spec <- get_ps_spec(ps_spec_id)
  data_file <- misspecified_model_all_data_file.list[[setting]][[miss]][[1]]
  all_data <- readRDS(data_file)
  mu_true <- get_mu_true(setting)

  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL),
    optimizer = optimizer_name
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
      psi1 <- gmm_fit$psi
      psi2_mat <- gmm_fit$psi2

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
      alpha2_vals[rep_i] <- gmm_fit$estimates[2]
    }, error = function(e) NULL)
  }

  valid1 <- !is.na(mu_vals) & !is.na(se1_vals)
  valid2 <- !is.na(mu_vals) & !is.na(se2_vals)
  esd <- sd(mu_vals[valid1])
  ese1 <- if(sum(valid1)>0) mean(se1_vals[valid1]) else NA
  ese2 <- if(sum(valid2)>0) mean(se2_vals[valid2]) else NA

  cat(sprintf("  %-15s: V1=%d, V2=%d, Bias=%.4f, ESD=%.4f, ESE1/ESD=%.3f, ESE2/ESD=%.3f\n",
      optimizer_name, sum(valid1), sum(valid2),
      mean(mu_vals[valid1]) - mu_true, esd,
      ese1/esd, if(sum(valid2)>0) ese2/sd(mu_vals[valid2]) else NA))
}

cat("=== Setting4, model 3 ===\n")
run_test("setting4", "miss50", "9-alt1", 3, "L-BFGS-B")
run_test("setting4", "miss50", "9-alt1", 3, "constrained_nr")

cat("\n=== Setting4, model 2 ===\n")
run_test("setting4", "miss50", "9-alt1", 2, "L-BFGS-B")
run_test("setting4", "miss50", "9-alt1", 2, "constrained_nr")

cat("\n=== Setting3, model 2 ===\n")
run_test("setting3", "miss50", "9-alt1", 2, "L-BFGS-B")
run_test("setting3", "miss50", "9-alt1", 2, "constrained_nr")
