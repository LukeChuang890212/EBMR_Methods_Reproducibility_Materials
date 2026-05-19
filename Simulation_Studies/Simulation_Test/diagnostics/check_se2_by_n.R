setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EPS", quiet = TRUE)

ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
mu_true <- get_mu_true("setting4")
total_rows <- nrow(all_data)
cat(sprintf("Total rows: %d\n\n", total_rows))

W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
model_idx <- 3

# Test different sample sizes
for (n_val in c(500, 1000, 2000, 4000)) {
  n_reps <- min(total_rows %/% n_val, 200)

  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL),
    optimizer = "L-BFGS-B"
  )

  mu_vals <- rep(NA_real_, n_reps)
  se2_mu <- rep(NA_real_, n_reps)
  se_simple_mu <- rep(NA_real_, n_reps)
  eta_max_vals <- rep(NA_real_, n_reps)

  for (rep_i in 1:n_reps) {
    dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
    tryCatch({
      ebmr <- EPS$new("y", single_ps, dat, W_fn)
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

      # Package SE2
      if (!is.null(psi2_mat) && is.matrix(psi2_mat) && !any(is.na(psi2_mat))) {
        mu_iid <- as.vector(t(r_vec/pi_hat*y_vec) - t(H_alpha) %*% psi2_mat)
        se2_mu[rep_i] <- sqrt(var(mu_iid)/n_val)
      }

      # Simple SE
      Gamma_hat <- gmm_fit$Gamma.hat
      g_mat <- gmm_fit$g.matrix
      W_hat <- gmm_fit$W.hat
      p <- ncol(design_mat)
      GtW <- crossprod(Gamma_hat, W_hat)
      GtWG <- GtW %*% Gamma_hat
      H_inv <- tryCatch(solve(GtWG), error = function(e) NULL)
      if (!is.null(H_inv)) {
        psi1_s <- H_inv %*% GtW %*% t(g_mat)
        dmu <- colMeans(-r_vec*y_vec*(1-pi_hat)/pi_hat * design_mat)
        iid_s <- r_vec*y_vec/pi_hat - mu_hat - as.vector(t(dmu) %*% psi1_s)
        se_simple_mu[rep_i] <- sqrt(mean(iid_s^2)/n_val)
      }

      mu_vals[rep_i] <- mu_hat
      eta_vals <- as.vector(design_mat %*% gmm_fit$estimates)
      eta_max_vals[rep_i] <- max(abs(eta_vals))
    }, error = function(e) NULL)
  }

  v2 <- !is.na(mu_vals) & !is.na(se2_mu)
  vs <- !is.na(mu_vals) & !is.na(se_simple_mu)

  cat(sprintf("n=%4d, reps=%3d: pkg_valid=%3d, simple_valid=%3d\n", n_val, n_reps, sum(v2), sum(vs)))
  if (sum(v2) > 10) {
    cat(sprintf("  ESE2/ESD=%.3f, Simple_ESE/ESD=%.3f, pct_eta>20=%.1f%%\n",
        mean(se2_mu[v2])/sd(mu_vals[v2]),
        mean(se_simple_mu[vs])/sd(mu_vals[vs]),
        100*mean(eta_max_vals[!is.na(mu_vals)] > 20)))
  }
}
