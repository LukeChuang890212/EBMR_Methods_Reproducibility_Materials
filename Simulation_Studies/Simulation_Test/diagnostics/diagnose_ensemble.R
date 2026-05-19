setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
source("Simulation.r")
library(EPS)

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

# Look at a few outlier reps (from inspect_test59: reps 8, 27, 28, 40, 41)
outlier_reps <- c(8, 27, 28, 40, 41, 83, 92)
normal_reps <- c(1, 2, 3, 4, 5)

cat("=== DIAGNOSING ENSEMBLE GMM ===\n\n")

for (rep_i in c(normal_reps, outlier_reps)) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]

  ebmr <- EPS$new("y", subset_ps_spec, dat, W_func)
  res <- ebmr$EBMR_IPW(
    h_nu = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2),
    se.fit = TRUE, type = "HT"
  )

  is_outlier <- rep_i %in% outlier_reps
  cat(sprintf("Rep %d [%s]: nu=(%.4f, %.4f), w=(%.4f, %.4f), mu=%.4f, obj=%.6f, converged=%s\n",
      rep_i, ifelse(is_outlier, "OUTLIER", "normal"),
      res$nu.hat[1], res$nu.hat[2],
      res$w.hat[1], res$w.hat[2],
      res$mu_ipw,
      res$ensemble_fit$gmm_fit$opt$objective,
      res$ensemble_fit$gmm_fit$opt$converged))

  # Also try manually with different inits to see if there's a better solution
  # Access the ps.matrix
  ps.matrix <- do.call(cbind, lapply(ebmr$ps_fit.list[c(1,3)], function(f) f$fitted.values))

  # Compute GMM objective for both solutions: nu=(1,0) and nu=(0,1)
  r_vec <- as.vector(dat$r)
  h_x_raw <- cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2)
  h_x <- cbind(1, h_x_raw)

  for (nu_test in list(c(1, 0), c(0, 1), res$nu.hat)) {
    w_test <- nu_test^2 / sum(nu_test^2)
    ps_nu <- as.vector(ps.matrix %*% w_test)
    g_mat <- (r_vec / ps_nu - 1) * h_x
    g_bar <- colMeans(g_mat)
    W_mat <- solve(t(g_mat) %*% g_mat / nrow(g_mat))
    obj <- as.numeric(t(g_bar) %*% W_mat %*% g_bar) * nrow(g_mat)
    cat(sprintf("  nu=(%.2f,%.2f) -> w=(%.4f,%.4f): Q=%.6f\n",
        nu_test[1], nu_test[2], w_test[1], w_test[2], obj))
  }
  cat("\n")
}
