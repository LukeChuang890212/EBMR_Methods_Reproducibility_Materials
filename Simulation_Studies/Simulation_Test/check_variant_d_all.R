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

run_variant_d <- function(model_idx, optimizer_name) {
  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL),
    optimizer = optimizer_name
  )

  mu_vals <- rep(NA_real_, n_reps)
  se_pkg <- rep(NA_real_, n_reps)
  se_D <- rep(NA_real_, n_reps)

  for (rep_i in 1:n_reps) {
    dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
    tryCatch({
      ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
      r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
      pi_hat <- ebmr$ps_fit.list[[1]]$fitted.values
      mu_hat <- mean(r_vec * y_vec / pi_hat)
      n <- n_val

      gmm_fit <- ebmr$ps_fit.list[[1]]$gmm_fit
      design_mat <- ebmr$ps_fit.list[[1]]$design_matrix
      link_type <- ebmr$ps_fit.list[[1]]$link_type

      if (!is.null(link_type) && link_type == "logistic_complement") {
        dot_pi <- -design_mat * (pi_hat * (1 - pi_hat))
      } else {
        dot_pi <- design_mat * (pi_hat * (1 - pi_hat))
      }
      ry_ps_inv2 <- as.vector(r_vec * y_vec * (pi_hat^(-2)))
      H_alpha <- colMeans(dot_pi * ry_ps_inv2)

      # Package SE2
      psi2_mat <- gmm_fit$psi2
      if (!is.null(psi2_mat) && is.matrix(psi2_mat) && !any(is.na(psi2_mat))) {
        mu_iid_pkg <- as.vector(t(r_vec/pi_hat*y_vec) - t(H_alpha) %*% psi2_mat)
        se_pkg[rep_i] <- sqrt(var(mu_iid_pkg)/n)
      }

      # Variant D: simple psi = -(GtWG)^{-1} GtW g_i
      Gamma_hat <- gmm_fit$Gamma.hat
      g_mat <- gmm_fit$g.matrix
      W_hat <- gmm_fit$W.hat
      GtW <- crossprod(Gamma_hat, W_hat)
      GtWG <- GtW %*% Gamma_hat
      GtWG_inv <- tryCatch(solve(GtWG), error = function(e) NULL)
      if (!is.null(GtWG_inv)) {
        psi_D <- -(GtWG_inv %*% GtW %*% t(g_mat))
        mu_iid_D <- as.vector(t(r_vec/pi_hat*y_vec) - t(H_alpha) %*% psi_D)
        se_D[rep_i] <- sqrt(var(mu_iid_D)/n)
      }

      mu_vals[rep_i] <- mu_hat
    }, error = function(e) NULL)
  }

  v_pkg <- !is.na(mu_vals) & !is.na(se_pkg)
  v_D <- !is.na(mu_vals) & !is.na(se_D)

  cat(sprintf("  M%d %-15s: pkg ESE2/ESD=%.3f (v=%d), D ESE/ESD=%.3f (v=%d)\n",
      model_idx, optimizer_name,
      if(sum(v_pkg)>10) mean(se_pkg[v_pkg])/sd(mu_vals[v_pkg]) else NA,
      sum(v_pkg),
      if(sum(v_D)>10) mean(se_D[v_D])/sd(mu_vals[v_D]) else NA,
      sum(v_D)))
}

cat("=== Package SE2 vs Variant D across all models ===\n\n")
for (m in 1:3) {
  for (opt in c("L-BFGS-B", "constrained_nr")) {
    run_variant_d(m, opt)
  }
}
