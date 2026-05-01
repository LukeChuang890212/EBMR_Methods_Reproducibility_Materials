setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
library(EBMRalgorithmFast4)
source("Simulation.r")

ps_spec <- get_ps_spec("9-alt1")
ps_spec_13 <- list(
  formula.list = ps_spec[["formula.list"]][c(1, 3)],
  h_alpha.list = ps_spec[["h_alpha.list"]][c(1, 3)],
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]]
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
n_val <- 2000
mu_true <- get_mu_true("setting4")

for (rep_i in c(1, 2, 3, 5, 10, 15, 20)) {
  set.seed(12345 + rep_i)
  dat <- setting4.B1(n_val)

  ebmr <- EBMRAlgorithmFast4[["new"]]("y", ps_spec_13, dat, W_func)

  cat(sprintf("\n=== Rep %d ===\n", rep_i))
  for (j in 1:2) {
    ps_fit <- ebmr[["ps_fit.list"]][[j]]
    alpha <- ps_fit[["coefficients"]]
    pi_hat <- ps_fit[["fitted.values"]]
    gmm_fit <- ps_fit[["gmm_fit"]]
    hpd <- gmm_fit[["opt"]][["hessian_pd"]]
    cat(sprintf("  Model %d: ||alpha||=%.2f, obj=%.6f, hessian_pd=%s, pi_bnd=%d\n",
        j, sqrt(sum(alpha^2)), gmm_fit[["opt"]][["objective"]],
        ifelse(is.null(hpd), "NULL", as.character(hpd)),
        sum(pi_hat < 0.06 | pi_hat > 0.94)))
  }

  res <- ebmr[["EBMR_IPW"]](
    h_nu = function(dat) cbind(u1 = dat[["u1"]], u2 = dat[["u2"]], z1 = dat[["z1"]], z2 = dat[["z2"]], u1_u2 = dat[["u1"]]*dat[["u2"]]),
    se.fit = FALSE, type = "HT"
  )
  cat(sprintf("  Ensemble: w=(%s), mu=%.4f (true=%.4f)\n",
      paste(round(res[["w.hat"]], 4), collapse=", "),
      res[["mu_ipw"]], mu_true))
}
