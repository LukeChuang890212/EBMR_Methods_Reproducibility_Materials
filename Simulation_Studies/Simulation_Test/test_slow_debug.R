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

for (rep_i in c(1, 2, 3, 5, 10)) {
  set.seed(12345 + rep_i)
  dat <- setting4.B1(n_val)

  cat(sprintf("\n=== Rep %d ===\n", rep_i))

  # Time PS estimation (step 1)
  t1 <- system.time({
    ebmr <- EBMRAlgorithmFast4[["new"]]("y", ps_spec_13, dat, W_func)
  })
  cat(sprintf("PS estimation: %.2f sec\n", t1[3]))

  # Check PS fit details
  for (j in 1:2) {
    ps_fit <- ebmr[["ps_fit.list"]][[j]]
    gmm_fit <- ps_fit[["gmm_fit"]]
    cat(sprintf("  Model %d: ||alpha||=%.4f, obj=%.6f, iter=%d, converged=%s, grad_norm=%.2e\n",
        j, sqrt(sum(ps_fit[["coefficients"]]^2)),
        gmm_fit[["opt"]][["objective"]],
        gmm_fit[["opt"]][["iterations"]],
        gmm_fit[["opt"]][["converged"]],
        gmm_fit[["opt"]][["final_grad_norm"]]))
  }

  # Time ensemble step
  t2 <- system.time({
    res <- ebmr[["EBMR_IPW"]](
      h_nu = function(dat) cbind(u1 = dat[["u1"]], u2 = dat[["u2"]], z1 = dat[["z1"]], z2 = dat[["z2"]], u1_u2 = dat[["u1"]]*dat[["u2"]]),
      se.fit = TRUE, type = "HT"
    )
  })
  cat(sprintf("EBMR_IPW (ensemble+SE): %.2f sec\n", t2[3]))
  cat(sprintf("  mu_ipw=%.4f, se_ipw=%.4f\n", res[["mu_ipw"]], res[["se_ipw"]]))
  cat(sprintf("  nu.hat: %s\n", paste(round(res[["nu.hat"]], 4), collapse=", ")))
  cat(sprintf("  w.hat:  %s\n", paste(round(res[["w.hat"]], 4), collapse=", ")))

  # Check ensemble GMM details
  # Re-run ensemble to get timing breakdown
  ps.matrix <- do.call(cbind, lapply(ebmr[["ps_fit.list"]], function(x) x[["fitted.values"]]))
  cat(sprintf("  PS range model1: [%.4f, %.4f]\n", min(ps.matrix[,1]), max(ps.matrix[,1])))
  cat(sprintf("  PS range model2: [%.4f, %.4f]\n", min(ps.matrix[,2]), max(ps.matrix[,2])))
}
