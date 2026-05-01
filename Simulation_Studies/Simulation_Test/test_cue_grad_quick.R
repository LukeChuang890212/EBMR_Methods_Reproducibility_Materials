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

# Test on outlier reps (from earlier diagnosis) and a few normal reps
test_reps <- c(1, 2, 3, 8, 28, 83, 92, 133, 155)

cat("=== QUICK TEST: Full CUE gradient convergence ===\n\n")
for (rep_i in test_reps) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]
  is_outlier <- rep_i %in% c(8, 28, 83, 92, 133, 155)

  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", subset_ps_spec, dat, W_func)
    res <- ebmr$EBMR_IPW(
      h_nu = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2),
      se.fit = TRUE, type = "HT"
    )

    alpha2_norm <- sqrt(sum(ebmr$ps_fit.list[[2]]$coefficients^2))
    cat(sprintf("Rep %3d [%s]: nu=(%.4f,%.4f) w=(%.4f,%.4f) mu=%.4f ||a2||=%.2f conv=%s iter=%d grad=%.2e\n",
        rep_i, ifelse(is_outlier, "OUTLIER", "normal"),
        res$nu.hat[1], res$nu.hat[2],
        res$w.hat[1], res$w.hat[2],
        res$mu_ipw, alpha2_norm,
        ebmr$ps_fit.list[[2]]$opt$converged,
        ebmr$ps_fit.list[[2]]$opt$iterations,
        ebmr$ps_fit.list[[2]]$opt$final_grad_norm))
  }, error = function(e) {
    cat(sprintf("Rep %3d [%s]: ERROR - %s\n", rep_i, ifelse(is_outlier, "OUTLIER", "normal"), e$message))
  })
}
