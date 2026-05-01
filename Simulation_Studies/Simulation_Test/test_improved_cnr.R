setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

n_val <- 2000
n_reps <- 200
ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
mu_true <- get_mu_true("setting4")
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

run_test <- function(model_idx, optimizer_name) {
  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL),
    optimizer = optimizer_name
  )

  mu_vals <- se2_mu <- foc_vals <- rep(NA_real_, n_reps)
  converged_flag <- rep(FALSE, n_reps)

  for (rep_i in 1:n_reps) {
    dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
    tryCatch({
      ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
      r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
      pi_hat <- ebmr$ps_fit.list[[1]]$fitted.values
      mu_hat <- mean(r_vec * y_vec / pi_hat)

      gmm_fit <- ebmr$ps_fit.list[[1]]$gmm_fit
      psi2_mat <- gmm_fit$psi2
      design_mat <- ebmr$ps_fit.list[[1]]$design_matrix
      link_type <- ebmr$ps_fit.list[[1]]$link_type

      GtW <- crossprod(gmm_fit$Gamma.hat, gmm_fit$W.hat)
      foc_vals[rep_i] <- max(abs(GtW %*% gmm_fit$eta_s))
      converged_flag[rep_i] <- gmm_fit$opt$converged

      if (!is.null(link_type) && link_type == "logistic_complement") {
        dot_pi <- -design_mat * (pi_hat * (1 - pi_hat))
      } else {
        dot_pi <- design_mat * (pi_hat * (1 - pi_hat))
      }
      ry_ps_inv2 <- as.vector(r_vec * y_vec * (pi_hat^(-2)))
      H_alpha <- colMeans(dot_pi * ry_ps_inv2)

      if (!is.null(psi2_mat) && is.matrix(psi2_mat) && !any(is.na(psi2_mat))) {
        mu_iid <- as.vector(t(r_vec/pi_hat*y_vec) - t(H_alpha) %*% psi2_mat)
        se2_mu[rep_i] <- sqrt(var(mu_iid)/n_val)
      }
      mu_vals[rep_i] <- mu_hat
    }, error = function(e) NULL)
  }

  valid <- !is.na(mu_vals) & !is.na(se2_mu)
  conv <- valid & converged_flag

  cat(sprintf("\n  M%d %-15s:\n", model_idx, optimizer_name))
  cat(sprintf("    All valid:  n=%3d, Bias=%.4f, ESD=%.4f, ESE2/ESD=%.3f\n",
      sum(valid), mean(mu_vals[valid]) - mu_true, sd(mu_vals[valid]),
      mean(se2_mu[valid]) / sd(mu_vals[valid])))
  if (sum(conv) > 10) {
    cat(sprintf("    Converged:  n=%3d, Bias=%.4f, ESD=%.4f, ESE2/ESD=%.3f\n",
        sum(conv), mean(mu_vals[conv]) - mu_true, sd(mu_vals[conv]),
        mean(se2_mu[conv]) / sd(mu_vals[conv])))
  } else {
    cat(sprintf("    Converged:  n=%3d (too few)\n", sum(conv)))
  }
  cat(sprintf("    FOC: median=%.2e, mean=%.2e\n",
      median(foc_vals[valid]), mean(foc_vals[valid])))
}

cat("=== Improved constrained NR ===\n")
for (m in 1:3) {
  for (opt in c("L-BFGS-B", "constrained_nr")) {
    run_test(m, opt)
  }
}
