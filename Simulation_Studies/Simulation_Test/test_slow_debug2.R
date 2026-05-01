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

for (rep_i in c(1, 2, 3)) {
  set.seed(12345 + rep_i)
  dat <- setting4.B1(n_val)
  cat(sprintf("\n=== Rep %d ===\n", rep_i))

  ebmr <- EBMRAlgorithmFast4[["new"]]("y", ps_spec_13, dat, W_func)

  # Time ensemble without SE
  t1 <- system.time({
    res_nose <- ebmr[["EBMR_IPW"]](
      h_nu = function(dat) cbind(u1 = dat[["u1"]], u2 = dat[["u2"]], z1 = dat[["z1"]], z2 = dat[["z2"]], u1_u2 = dat[["u1"]]*dat[["u2"]]),
      se.fit = FALSE, type = "HT"
    )
  })
  cat(sprintf("EBMR_IPW (no SE):  %.2f sec, mu=%.4f\n", t1[3], res_nose[["mu_ipw"]]))

  # Rebuild for SE computation
  ebmr2 <- EBMRAlgorithmFast4[["new"]]("y", ps_spec_13, dat, W_func)
  t2 <- system.time({
    res_se <- ebmr2[["EBMR_IPW"]](
      h_nu = function(dat) cbind(u1 = dat[["u1"]], u2 = dat[["u2"]], z1 = dat[["z1"]], z2 = dat[["z2"]], u1_u2 = dat[["u1"]]*dat[["u2"]]),
      se.fit = TRUE, type = "HT"
    )
  })
  cat(sprintf("EBMR_IPW (with SE): %.2f sec, mu=%.4f, se=%.4f\n", t2[3], res_se[["mu_ipw"]], res_se[["se_ipw"]]))
  cat(sprintf("SE computation overhead: %.2f sec\n", t2[3] - t1[3]))
}
