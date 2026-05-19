setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EPS", quiet = TRUE)

n_val <- 2000
n_reps <- 200

ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
mu_true <- get_mu_true("setting4")

W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

diagnose <- function(model_idx, optimizer_name) {
  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL),
    optimizer = optimizer_name
  )

  mu_vals <- se2_vals <- rep(NA_real_, n_reps)
  foc_norm <- G_norm <- obj_vals <- rep(NA_real_, n_reps)
  H2_extra_norm <- H2_cond <- rep(NA_real_, n_reps)
  Q_extra_var <- rep(NA_real_, n_reps)
  converged <- rep(NA, n_reps)

  for (rep_i in 1:n_reps) {
    dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
    tryCatch({
      ebmr <- EPS$new("y", single_ps, dat, W_fn)
      r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
      pi_hat <- ebmr$ps_fit.list[[1]]$fitted.values
      mu_hat <- mean(r_vec * y_vec / pi_hat)

      gmm_fit <- ebmr$ps_fit.list[[1]]$gmm_fit
      design_mat <- ebmr$ps_fit.list[[1]]$design_matrix
      link_type <- ebmr$ps_fit.list[[1]]$link_type

      Gamma_hat <- gmm_fit$Gamma.hat
      g_mat <- gmm_fit$g.matrix
      W_hat <- gmm_fit$W.hat
      eta_s <- gmm_fit$eta_s
      GtW <- crossprod(Gamma_hat, W_hat)
      GtWG <- GtW %*% Gamma_hat

      # Convergence diagnostics
      foc_norm[rep_i] <- max(abs(GtW %*% eta_s))
      G_norm[rep_i] <- sqrt(sum(eta_s^2))
      obj_vals[rep_i] <- as.numeric(crossprod(eta_s, W_hat %*% eta_s))
      converged[rep_i] <- gmm_fit$opt$converged

      # SE2 components
      if (!is.null(link_type) && link_type == "logistic_complement") {
        dot_pi <- -design_mat * (pi_hat * (1 - pi_hat))
      } else {
        dot_pi <- design_mat * (pi_hat * (1 - pi_hat))
      }
      ry_ps_inv2 <- as.vector(r_vec * y_vec * (pi_hat^(-2)))
      H_alpha <- colMeans(dot_pi * ry_ps_inv2)

      psi2_mat <- gmm_fit$psi2
      if (!is.null(psi2_mat) && is.matrix(psi2_mat) && !any(is.na(psi2_mat))) {
        mu_iid <- as.vector(t(r_vec/pi_hat*y_vec) - t(H_alpha) %*% psi2_mat)
        se2_vals[rep_i] <- sqrt(var(mu_iid)/n_val)
      }

      # Measure the "extra" in H2 beyond GtWG
      # H2 = GtWG + (G'W⊗I)R - (G'W⊗Γ'W)S
      # The extra terms magnitude relative to GtWG
      GW_row <- as.vector(t(eta_s) %*% W_hat)
      H2_extra_norm[rep_i] <- sqrt(sum(GW_row^2))  # ~ |G'W|, scales the extra terms

      # Measure variance contribution from Q_full extra terms
      # Q_i = Γ'Wg_i + [Γ_i'WG - (Γ'Wg_i)(g_i'WG)]
      # The extra part is: Γ_i'WG - (Γ'Wg_i)(g_i'WG)
      WG <- as.vector(W_hat %*% eta_s)
      g_WG <- as.vector(g_mat %*% WG)  # g_i'WG for each i
      GtW_g <- GtW %*% t(g_mat)  # k x n

      # The extra terms introduce variance proportional to var(Γ_i'WG) and var((Γ'Wg_i)(g_i'WG))
      # Simple proxy: ratio of var of extra vs base
      base_var <- mean(apply(GtW_g, 1, var))
      extra_var <- var(g_WG) * mean(apply(GtW_g, 1, var))  # approximate
      Q_extra_var[rep_i] <- sqrt(extra_var / max(base_var, 1e-20))

      mu_vals[rep_i] <- mu_hat
    }, error = function(e) NULL)
  }

  valid <- !is.na(mu_vals) & !is.na(se2_vals)
  cat(sprintf("\n=== M%d %s (valid=%d) ===\n", model_idx, optimizer_name, sum(valid)))
  cat(sprintf("  Bias=%.4f, ESD=%.4f, ESE2/ESD=%.3f\n",
      mean(mu_vals[valid]) - mu_true, sd(mu_vals[valid]),
      mean(se2_vals[valid]) / sd(mu_vals[valid])))
  cat(sprintf("  FOC norm: median=%.2e, mean=%.2e, max=%.2e\n",
      median(foc_norm[valid]), mean(foc_norm[valid]), max(foc_norm[valid])))
  cat(sprintf("  |G| norm: median=%.4f, mean=%.4f, max=%.4f\n",
      median(G_norm[valid]), mean(G_norm[valid]), max(G_norm[valid])))
  cat(sprintf("  Objective: median=%.2e, mean=%.2e\n",
      median(obj_vals[valid]), mean(obj_vals[valid])))
  cat(sprintf("  |G'W| (H2 extra scale): median=%.4f, mean=%.4f\n",
      median(H2_extra_norm[valid]), mean(H2_extra_norm[valid])))
  cat(sprintf("  Q extra/base var ratio: median=%.4f, mean=%.4f\n",
      median(Q_extra_var[valid]), mean(Q_extra_var[valid])))
  cat(sprintf("  Converged: %d/%d\n", sum(converged[valid], na.rm=TRUE), sum(valid)))
}

for (m in 1:3) {
  for (opt in c("L-BFGS-B", "constrained_nr")) {
    diagnose(m, opt)
  }
}
