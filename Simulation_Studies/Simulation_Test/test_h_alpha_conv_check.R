setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

n_val <- 2000; n_reps <- 200; m_idx <- 3
ps_spec_orig <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
mu_true <- get_mu_true("setting4")
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

h_alpha_enrich <- function(dat) {
  cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2,
        u2_sq = dat$u2^2, z2_sq = dat$z2^2,
        u1_u2 = dat$u1*dat$u2, u1_z2 = dat$u1*dat$z2,
        z1_u2 = dat$z1*dat$u2, z1_z2 = dat$z1*dat$z2,
        u2_z2 = dat$u2*dat$z2)
}

mu_vals <- se2_vals <- grad_vals <- rep(NA_real_, n_reps)
for (rep_i in 1:n_reps) {
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
  single_ps <- list(
    formula.list = list(ps_spec_orig[["formula.list"]][[m_idx]]),
    h_alpha.list = list(h_alpha_enrich),
    inv_link = ps_spec_orig[["inv_link"]],
    outcome = ps_spec_orig[["outcome"]],
    alpha_init.list = list(NULL),
    optimizer = "constrained_nr"
  )
  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
    gmm_fit <- ebmr$ps_fit.list[[1]]$gmm_fit
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
  }, error = function(e) {})
}

valid <- !is.na(mu_vals) & !is.na(se2_vals)
conv <- valid & grad_vals < 1e-6

cat(sprintf("Converged (grad<1e-6): n=%d\n", sum(conv)))
cat(sprintf("  Bias     = %.4f\n", mean(mu_vals[conv]) - mu_true))
cat(sprintf("  ESD      = %.4f\n", sd(mu_vals[conv])))
cat(sprintf("  ESE2     = %.4f\n", mean(se2_vals[conv])))
cat(sprintf("  ESE2/ESD = %.3f\n", mean(se2_vals[conv]) / sd(mu_vals[conv])))

nonconv <- valid & grad_vals >= 1e-6
cat(sprintf("\nNon-converged: n=%d\n", sum(nonconv)))
if (sum(nonconv) > 1) {
  cat(sprintf("  Bias     = %.4f\n", mean(mu_vals[nonconv]) - mu_true))
  cat(sprintf("  ESD      = %.4f\n", sd(mu_vals[nonconv])))
  cat(sprintf("  ESE2     = %.4f\n", mean(se2_vals[nonconv])))
  cat(sprintf("  ESE2/ESD = %.3f\n", mean(se2_vals[nonconv]) / sd(mu_vals[nonconv])))
}
