setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Data_Generation.r")
source("config/scenarios.R")
library(EBMRalgorithmFast4)
source("Basic_setup.r")

data_file <- misspecified_model_all_data_file.list$setting3$miss50[[1]]
all_data <- readRDS(data_file)
n_val <- 2000

ps_spec <- get_ps_spec("9-alt1")
subset_ps_spec <- list(
  formula.list = ps_spec$formula.list[2],
  h_alpha.list = ps_spec$h_alpha.list[2],
  inv_link = ps_spec$inv_link,
  outcome = ps_spec$outcome
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

for (rep_i in c(437, 1)) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]
  ebmr <- EBMRAlgorithmFast4$new("y", subset_ps_spec, dat, W_func)
  ps_fit <- ebmr$ps_fit.list[[1]]
  ps_fitted <- ps_fit$fitted.values
  r <- as.numeric(dat$r)
  y <- dat$y
  n <- length(y)
  gmm_fit <- ps_fit$gmm_fit

  h_alpha_mat <- as.matrix(ps_fit$design_matrix)
  dot_pi <- h_alpha_mat * ps_fitted * (1 - ps_fitted)

  cat(sprintf("========== Rep %d ==========\n", rep_i))

  # H_alpha_w = mean(dot_pi_i * r_i * y_i / ps_i^2)
  # = mean(h_alpha_i * ps_i*(1-ps_i) * r_i * y_i / ps_i^2)
  # = mean(h_alpha_i * (1-ps_i)/ps_i * r_i * y_i)
  # The factor (1-ps_i)/ps_i * r_i * y_i blows up when ps_i is small

  factor_i <- (1 - ps_fitted) / ps_fitted * r * y  # length n
  cat(sprintf("Factor (1-ps)/ps * r * y:\n"))
  cat(sprintf("  mean=%.4f, sd=%.4f, max=%.2f\n", mean(factor_i), sd(factor_i), max(factor_i)))

  # How many obs contribute > 50%% of H_alpha_w?
  # Each obs contribution to H_alpha_w is: dot_pi_i * r_i*y_i/ps_i^2 / n
  contrib <- dot_pi * as.vector(r * y * ps_fitted^(-2)) / n  # n x alpha_dim
  H_alpha_w <- colSums(contrib)
  cat(sprintf("H_alpha_w = [%s]\n", paste(round(H_alpha_w, 4), collapse=", ")))

  # Focus on H_alpha_w[2] since it's the largest component
  j_max <- which.max(abs(H_alpha_w))
  cat(sprintf("Largest component: H_alpha_w[%d] = %.4f\n", j_max, H_alpha_w[j_max]))

  contrib_j <- contrib[, j_max]
  top_contrib <- order(abs(contrib_j), decreasing = TRUE)

  cat(sprintf("\nTop 10 obs contributions to H_alpha_w[%d]:\n", j_max))
  cumsum_pct <- 0
  for (k in 1:10) {
    i <- top_contrib[k]
    pct <- contrib_j[i] / H_alpha_w[j_max] * 100
    cumsum_pct <- cumsum_pct + pct
    cat(sprintf("  obs %d: contrib=%.4f (%.1f%%, cum=%.1f%%), ps=%.5f, r=%d, y=%.3f, r*y/ps^2=%.1f\n",
                i, contrib_j[i], pct, cumsum_pct,
                ps_fitted[i], r[i], y[i], r[i]*y[i]/ps_fitted[i]^2))
  }

  # Now look at psi_alpha for these observations
  # psi_alpha_i = (GtWG)^{-1} GtW h_i where h_i = (r_i/ps_i - 1) * h_alpha_i
  psi_alpha <- gmm_fit$psi  # h_dim x n
  Gamma_hat <- gmm_fit$Gamma.hat
  W_hat <- gmm_fit$W.hat
  GtWG_inv <- solve(crossprod(Gamma_hat, W_hat) %*% Gamma_hat)

  cat(sprintf("\nGtWG_inv:\n"))
  print(round(GtWG_inv, 4))

  # Check the actual alpha_adj for obs 774 (rep 437) step by step
  obs_check <- top_contrib[1]
  h_i <- as.vector((r[obs_check]/ps_fitted[obs_check] - 1) * h_alpha_mat[obs_check,])
  psi_i <- as.vector(psi_alpha[, obs_check])

  cat(sprintf("\n--- Obs %d detailed ---\n", obs_check))
  cat(sprintf("  ps = %.6f, r/ps = %.2f\n", ps_fitted[obs_check], r[obs_check]/ps_fitted[obs_check]))
  cat(sprintf("  h_alpha_i = [%s]\n", paste(round(h_alpha_mat[obs_check,], 4), collapse=", ")))
  cat(sprintf("  h_i (moment cond) = [%s]\n", paste(round(h_i, 2), collapse=", ")))
  cat(sprintf("  ||h_i|| = %.2f\n", sqrt(sum(h_i^2))))
  cat(sprintf("  psi_i = [%s]\n", paste(round(psi_i, 2), collapse=", ")))
  cat(sprintf("  ||psi_i|| = %.2f\n", sqrt(sum(psi_i^2))))
  cat(sprintf("  alpha_adj_i = H' * psi_i = %.2f\n", sum(H_alpha_w * psi_i)))
  cat(sprintf("  base_i = r/ps * y = %.2f\n", r[obs_check]/ps_fitted[obs_check]*y[obs_check]))
  cat(sprintf("  iid_i = base - adj = %.2f\n",
              r[obs_check]/ps_fitted[obs_check]*y[obs_check] - sum(H_alpha_w * psi_i)))

  cat("\n\n")
}

cat("Done!\n")
