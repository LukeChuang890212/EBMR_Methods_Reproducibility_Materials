setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

n_val <- 2000
ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
model_idx <- 3

W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

# ============================
# Jackknife SE for mu_hat
# ============================
# Do leave-one-out on 3 reps and compare jackknife SE with package SE

for (rep_i in 1:3) {
  cat(sprintf("\n=== Rep %d ===\n", rep_i))
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]

  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL),
    optimizer = "constrained_nr"
  )

  # Full-sample fit
  ebmr_full <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
  gmm_fit <- ebmr_full$ps_fit.list[[1]]$gmm_fit
  alpha_full <- gmm_fit$estimates
  pi_full <- ebmr_full$ps_fit.list[[1]]$fitted.values
  r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
  mu_full <- mean(r_vec * y_vec / pi_full)
  n <- n_val

  # Package SE2 for mu
  psi2_mat <- gmm_fit$psi2
  design_mat <- ebmr_full$ps_fit.list[[1]]$design_matrix
  link_type <- ebmr_full$ps_fit.list[[1]]$link_type
  if (!is.null(link_type) && link_type == "logistic_complement") {
    dot_pi <- -design_mat * (pi_full * (1 - pi_full))
  } else {
    dot_pi <- design_mat * (pi_full * (1 - pi_full))
  }
  ry_ps_inv2 <- as.vector(r_vec * y_vec * (pi_full^(-2)))
  H_alpha <- colMeans(dot_pi * ry_ps_inv2)

  if (is.matrix(psi2_mat) && !any(is.na(psi2_mat))) {
    mu_iid_pkg <- as.vector(t(r_vec/pi_full*y_vec) - t(H_alpha) %*% psi2_mat)
    se_pkg <- sqrt(var(mu_iid_pkg)/n)
  } else {
    se_pkg <- NA
  }

  # Simple SE for mu
  Gamma_hat <- gmm_fit$Gamma.hat
  g_mat <- gmm_fit$g.matrix
  W_hat <- gmm_fit$W.hat
  GtW <- crossprod(Gamma_hat, W_hat)
  GtWG <- GtW %*% Gamma_hat
  H_inv <- solve(GtWG)
  psi1_simple <- H_inv %*% GtW %*% t(g_mat)
  dmu <- colMeans(-r_vec*y_vec*(1-pi_full)/pi_full * design_mat)
  iid_simple <- r_vec*y_vec/pi_full - mu_full - as.vector(t(dmu) %*% psi1_simple)
  se_simple <- sqrt(mean(iid_simple^2)/n)

  cat(sprintf("  Package SE2 for mu: %.6f\n", se_pkg))
  cat(sprintf("  Simple SE for mu:   %.6f\n", se_simple))

  # Jackknife: delete 100 random observations
  set.seed(42 + rep_i)
  jack_idx <- sample(1:n, 100)
  mu_jack <- rep(NA, 100)
  alpha_jack <- matrix(NA, 100, length(alpha_full))

  for (k in seq_along(jack_idx)) {
    i_del <- jack_idx[k]
    dat_del <- dat[-i_del, ]

    single_ps_del <- list(
      formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
      h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
      inv_link = ps_spec[["inv_link"]],
      outcome = ps_spec[["outcome"]],
      alpha_init.list = list(alpha_full),  # warm start
      optimizer = "constrained_nr"
    )

    tryCatch({
      ebmr_del <- EBMRAlgorithmFast4$new("y", single_ps_del, dat_del, W_fn)
      pi_del <- ebmr_del$ps_fit.list[[1]]$fitted.values
      r_del <- dat_del[["r"]]; y_del <- dat_del[["y"]]
      mu_jack[k] <- mean(r_del * y_del / pi_del)
      alpha_jack[k, ] <- ebmr_del$ps_fit.list[[1]]$gmm_fit$estimates
    }, error = function(e) NULL)
  }

  valid_jack <- !is.na(mu_jack)
  cat(sprintf("  Jackknife valid: %d/100\n", sum(valid_jack)))

  if (sum(valid_jack) > 10) {
    # Jackknife SE for mu
    jack_se_mu <- sqrt((n-1) * var(mu_jack[valid_jack]))
    cat(sprintf("  Jackknife SE for mu: %.6f\n", jack_se_mu))
    cat(sprintf("  Ratio pkg/jack: %.3f\n", se_pkg / jack_se_mu))
    cat(sprintf("  Ratio simple/jack: %.3f\n", se_simple / jack_se_mu))

    # Jackknife SE for alpha
    jack_se_alpha <- sqrt((n-1) * apply(alpha_jack[valid_jack, ], 2, var))
    cat(sprintf("  Jackknife SE alpha: %s\n",
        paste(sprintf("%.4f", jack_se_alpha), collapse=", ")))
    cat(sprintf("  Package SE2 alpha:  %s\n",
        paste(sprintf("%.4f", gmm_fit$se2), collapse=", ")))
  }
}
