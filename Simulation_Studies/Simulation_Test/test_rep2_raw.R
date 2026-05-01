setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
library(EBMRalgorithmFast4)

ps_spec <- get_ps_spec("9-alt1")
ps_spec_13 <- list(
  formula.list = ps_spec[["formula.list"]][c(1, 3)],
  h_alpha.list = ps_spec[["h_alpha.list"]][c(1, 3)],
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]]
)
W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

set.seed(12345 + 2)
dat <- setting4.B1(2000)
ebmr <- EBMRAlgorithmFast4[["new"]]("y", ps_spec_13, dat, W_func)

for (j in 1:2) {
  ps <- ebmr[["ps_fit.list"]][[j]]
  pi_hat <- ps[["fitted.values"]]
  cat(sprintf("\nModel %d:\n", j))
  cat("Alpha:", round(ps[["coefficients"]], 4), "\n")
  cat("||alpha||:", round(sqrt(sum(ps[["coefficients"]]^2)), 3), "\n")
  cat("Pi quantiles:", round(quantile(pi_hat, c(0, 0.01, 0.05, 0.10, 0.50, 0.90, 0.95, 0.99, 1)), 6), "\n")
  cat("Pi < 0.01:", sum(pi_hat < 0.01), "\n")
  cat("Pi > 0.99:", sum(pi_hat > 0.99), "\n")
  cat("Objective:", ps[["gmm_fit"]][["opt"]][["objective"]], "\n")
}
