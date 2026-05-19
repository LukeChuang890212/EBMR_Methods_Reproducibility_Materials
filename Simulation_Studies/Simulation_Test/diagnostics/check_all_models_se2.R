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

cat("PS spec 9-alt1:\n")
for (m in 1:3) {
  f <- ps_spec[["formula.list"]][[m]]
  cat(sprintf("  Model %d: formula=%s\n", m, deparse(f)))
}

# Check SE2 for each model
for (m in 1:3) {
  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[m]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[m]]),
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL),
    optimizer = "L-BFGS-B"
  )

  mu_vals <- rep(NA_real_, n_reps)
  se2_mu <- rep(NA_real_, n_reps)
  eta_max_vals <- rep(NA_real_, n_reps)
  alpha2_vals <- rep(NA_real_, n_reps)

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
      alpha_hat <- gmm_fit$estimates
      alpha2_vals[rep_i] <- alpha_hat[2]

      # Track max |eta|
      eta_vals <- as.vector(design_mat %*% alpha_hat)
      eta_max_vals[rep_i] <- max(abs(eta_vals))
    }, error = function(e) NULL)
  }

  valid <- !is.na(mu_vals) & !is.na(se2_mu)
  esd <- sd(mu_vals[valid])
  ese <- mean(se2_mu[valid])

  cat(sprintf("\nModel %d (L-BFGS-B): valid=%d, Bias=%.4f, ESD=%.4f, ESE2/ESD=%.3f\n",
      m, sum(valid), mean(mu_vals[valid]) - mu_true, esd, ese/esd))
  cat(sprintf("  alpha2: mean=%.3f, sd=%.3f\n", mean(alpha2_vals[valid]), sd(alpha2_vals[valid])))
  cat(sprintf("  max|eta|: mean=%.1f, max=%.1f, pct>20=%.1f%%\n",
      mean(eta_max_vals[valid]), max(eta_max_vals[valid]),
      100*mean(eta_max_vals[valid] > 20)))
}
