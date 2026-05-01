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

for (rep_i in c(1, 2, 3, 5, 10, 15, 20)) {
  set.seed(12345 + rep_i)
  dat <- setting4.B1(2000)
  ebmr <- EBMRAlgorithmFast4[["new"]]("y", ps_spec_13, dat, W_func)

  cat(sprintf("\n=== Rep %d ===\n", rep_i))
  for (j in 1:2) {
    ps_fit <- ebmr[["ps_fit.list"]][[j]]
    opt <- ps_fit[["gmm_fit"]][["opt"]]
    cat(sprintf("  Model %d: ||alpha||=%.2f, obj=%.6f, hessian_pd=%s, iters=%d, grad=%.2e, converged=%s\n",
        j, sqrt(sum(ps_fit[["coefficients"]]^2)),
        opt[["objective"]],
        as.character(opt[["hessian_pd"]]),
        opt[["iterations"]],
        opt[["final_grad_norm"]],
        as.character(opt[["converged"]])))
  }
}
