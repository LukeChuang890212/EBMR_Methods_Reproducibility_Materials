setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
source("Simulation.r")
library(EBMRalgorithmFast4)

data_file <- misspecified_model_all_data_file.list$setting4$miss50[[1]]
all_data <- readRDS(data_file)
n_val <- 2000

ps_spec <- get_ps_spec("9-alt1")
subset_ps_spec <- list(
  formula.list = ps_spec$formula.list[c(1, 3)],
  h_alpha.list = ps_spec$h_alpha.list[c(1, 3)],
  inv_link = ps_spec$inv_link,
  outcome = ps_spec$outcome
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

# Check a few outlier reps and normal reps
outlier_reps <- c(8, 28, 83, 92, 133, 155)
normal_reps <- c(1, 2, 3, 4, 5, 10)

cat("=== DIAGNOSING WRONG WEIGHT ASSIGNMENT ===\n\n")

for (rep_i in c(normal_reps, outlier_reps)) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]
  is_outlier <- rep_i %in% outlier_reps

  ebmr <- EBMRAlgorithmFast4$new("y", subset_ps_spec, dat, W_func)
  res <- ebmr$EBMR_IPW(
    h_nu = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2),
    se.fit = TRUE, type = "HT"
  )

  # Individual PS from each model
  ps1 <- ebmr$ps_fit.list[[1]]$fitted.values  # model 1 (correct)
  ps2 <- ebmr$ps_fit.list[[2]]$fitted.values  # model 3 (misspecified)

  # Alpha norms
  alpha1_norm <- sqrt(sum(ebmr$ps_fit.list[[1]]$coefficients^2))
  alpha2_norm <- sqrt(sum(ebmr$ps_fit.list[[2]]$coefficients^2))

  cat(sprintf("--- Rep %d [%s] ---\n", rep_i, ifelse(is_outlier, "OUTLIER", "normal")))
  cat(sprintf("  nu = (%.4f, %.4f), w = (%.4f, %.4f), mu = %.4f\n",
      res$nu.hat[1], res$nu.hat[2], res$w.hat[1], res$w.hat[2], res$mu_ipw))
  cat(sprintf("  ||alpha1|| = %.4f, ||alpha2|| = %.4f\n", alpha1_norm, alpha2_norm))

  cat(sprintf("  Model 1 (correct) PS: min=%.6f, q01=%.6f, q05=%.6f, median=%.4f, mean=%.4f\n",
      min(ps1), quantile(ps1, 0.01), quantile(ps1, 0.05), median(ps1), mean(ps1)))
  cat(sprintf("  Model 3 (misspec) PS: min=%.6f, q01=%.6f, q05=%.6f, median=%.4f, mean=%.4f\n",
      min(ps2), quantile(ps2, 0.01), quantile(ps2, 0.05), median(ps2), mean(ps2)))

  # Ensemble PS for both weight assignments
  ens_correct <- ps1  # w = (1, 0)
  ens_wrong <- ps2    # w = (0, 1)

  cat(sprintf("  If w=(1,0): min(ens_ps)=%.6f, mean(r/ens_ps)=%.4f\n",
      min(ens_correct), mean(dat$r / ens_correct)))
  cat(sprintf("  If w=(0,1): min(ens_ps)=%.6f, mean(r/ens_ps)=%.4f\n",
      min(ens_wrong), mean(dat$r / ens_wrong)))

  # GMM objective for both weight assignments
  r_vec <- as.vector(dat$r)
  h_x_raw <- cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2)
  h_x <- cbind(1, h_x_raw)

  # Objective with w=(1,0) using model 1 PS
  g1 <- (r_vec / ps1 - 1) * h_x
  g_bar1 <- colMeans(g1)
  W1 <- tryCatch(solve(t(g1) %*% g1 / nrow(g1)), error = function(e) diag(ncol(h_x)))
  obj1 <- as.numeric(t(g_bar1) %*% W1 %*% g_bar1) * nrow(g1)

  # Objective with w=(0,1) using model 3 PS
  g2 <- (r_vec / ps2 - 1) * h_x
  g_bar2 <- colMeans(g2)
  W2 <- tryCatch(solve(t(g2) %*% g2 / nrow(g2)), error = function(e) diag(ncol(h_x)))
  obj2 <- as.numeric(t(g_bar2) %*% W2 %*% g_bar2) * nrow(g2)

  cat(sprintf("  GMM obj w=(1,0): %.6f, w=(0,1): %.6f  [winner: %s]\n",
      obj1, obj2, ifelse(obj1 < obj2, "correct", "WRONG")))
  cat(sprintf("  Actual GMM obj: %.6f\n", res$ensemble_fit$gmm_fit$opt$objective))

  # Check moment conditions: E[(r/pi - 1) * h_x] should be ~0 for correct model
  cat(sprintf("  |g_bar| with model 1: %s\n", paste(sprintf("%.4f", g_bar1), collapse=", ")))
  cat(sprintf("  |g_bar| with model 3: %s\n", paste(sprintf("%.4f", g_bar2), collapse=", ")))
  cat("\n")
}
