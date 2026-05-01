setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
source("Simulation.r")
library(EBMRalgorithmFast4)

ps_spec <- get_ps_spec("9-alt1")
ps_spec_13 <- list(
  formula.list = ps_spec[["formula.list"]][c(1, 3)],
  h_alpha.list = ps_spec[["h_alpha.list"]][c(1, 3)],
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]]
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
n_val <- 2000

for (rep_i in c(1, 2, 3, 10)) {
  set.seed(12345 + rep_i)
  dat <- setting4.B1(n_val)

  ebmr <- EBMRAlgorithmFast4[["new"]]("y", ps_spec_13, dat, W_func)

  cat(sprintf("\n=== Rep %d ===\n", rep_i))
  for (j in 1:2) {
    ps_fit <- ebmr[["ps_fit.list"]][[j]]
    alpha <- ps_fit[["coefficients"]]
    pi_hat <- ps_fit[["fitted.values"]]
    cat(sprintf("  Model %d:\n", j))
    cat(sprintf("    alpha: %s\n", paste(round(alpha, 4), collapse=", ")))
    cat(sprintf("    ||alpha||=%.4f\n", sqrt(sum(alpha^2))))
    cat(sprintf("    pi range: [%.6f, %.6f]\n", min(pi_hat), max(pi_hat)))
    cat(sprintf("    pi quantiles: %.4f %.4f %.4f %.4f %.4f\n",
        quantile(pi_hat, 0.01), quantile(pi_hat, 0.05),
        quantile(pi_hat, 0.50), quantile(pi_hat, 0.95), quantile(pi_hat, 0.99)))
    cat(sprintf("    pi < 0.02: %d, pi > 0.98: %d\n", sum(pi_hat < 0.02), sum(pi_hat > 0.98)))
  }
}
