setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Data_Generation.r")
source("config/scenarios.R")
library(EBMRalgorithmFast4)
source("Basic_setup.r")

# Use a simple scenario to test
data_file <- misspecified_model_all_data_file.list$setting4$miss50[[1]]
all_data <- readRDS(data_file)
n_val <- 2000

ps_spec <- get_ps_spec("9-alt1")

# Test all J values
for (J_test in 1:3) {
  model_set <- 1:J_test
  cat("\n=== Testing J =", J_test, "===\n")
  subset_ps_spec <- list(
    formula.list = ps_spec$formula.list[model_set],
    h_alpha.list = ps_spec$h_alpha.list[model_set],
    inv_link = ps_spec$inv_link,
    outcome = ps_spec$outcome
  )

# Run a single replicate manually (not parallel)
i <- 1
dat <- all_data[((i - 1) * n_val + 1):(i * n_val), ]

W <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

tryCatch({
  ebmr <- EBMRAlgorithmFast4$new("y", subset_ps_spec, dat, W)
  cat("Number of PS models:", length(ebmr$ps_fit.list), "\n")

  result <- ebmr$EBMR_IPW(
    h_nu = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2),
    se.fit = TRUE
  )

  cat("mu_ipw:", result$mu_ipw, "\n")
  cat("se_ipw:", result$se_ipw, "\n")
  cat("nu.hat:", result$nu.hat, "\n")
  cat("w.hat:", result$w.hat, "\n")

  estimates <- unlist(result[1:4])
  cat("estimates:", estimates, "\n")

  alpha_out <- unlist(lapply(ebmr$ps_fit.list, function(ps_fit) ps_fit$coefficients))
  cat("alpha.hat length:", length(alpha_out), "\n")

  se_alpha_out <- unlist(lapply(ebmr$ps_fit.list, function(ps_fit) ps_fit$se))
  cat("se_alpha.hat length:", length(se_alpha_out), "\n")

  cat("nu.hat length:", length(result$nu.hat), "\n")
  cat("w.hat length:", length(result$w.hat), "\n")

}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  traceback()
})
}

cat("\nDone!\n")
