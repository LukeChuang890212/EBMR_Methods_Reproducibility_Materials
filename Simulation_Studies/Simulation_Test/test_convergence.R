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
n_val <- 2000

cat(sprintf("%5s %6s %10s %12s %10s %6s %12s %8s\n",
    "rep", "model", "||alpha||", "Q(alpha)", "hessian_pd", "iters", "grad_norm", "convgd"))

for (rep_i in c(1, 2, 3, 5, 10, 15, 20, 25, 30)) {
  set.seed(12345 + rep_i)
  dat <- setting4.B1(n_val)

  ebmr <- EBMRAlgorithmFast4[["new"]]("y", ps_spec_13, dat, W_func)

  for (j in 1:2) {
    ps_fit <- ebmr[["ps_fit.list"]][[j]]
    alpha <- ps_fit[["coefficients"]]
    gmm_fit <- ps_fit[["gmm_fit"]]
    opt <- gmm_fit[["opt"]]

    cat(sprintf("%5d %6d %10.3f %12.8f %10s %6d %12.2e %8s\n",
        rep_i, j, sqrt(sum(alpha^2)),
        opt[["objective"]],
        as.character(opt[["hessian_pd"]]),
        opt[["iterations"]],
        opt[["final_grad_norm"]],
        as.character(opt[["converged"]])))
  }
}
