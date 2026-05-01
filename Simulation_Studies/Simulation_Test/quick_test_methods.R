setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
source("Data_Generation.r")
source("config/scenarios.R")
old_wd <- getwd(); setwd("../EBMRalgorithmFast4"); source("R/EBMRAlgorithm.r"); setwd(old_wd)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
inv_link_fn <- function(eta) 1 / (1 + exp(eta))
ps_spec <- list(formula.list = list(r ~ y + u1 + u2),
                h_alpha.list = list(c("u1","u2","z1","z2")),
                inv_link = inv_link_fn, outcome = "y",
                alpha_init.list = list(NULL))

set.seed(42); dat <- setting3.A1(2000)
for (m in c("iterated", "gd")) {
  ebmr <- EBMRAlgorithmFast4$new("y", ps_spec, dat, W_func, bounds = 5, gmm_method = m)
  fit <- ebmr$ps_fit.list[[1]]
  cat(sprintf("%s: alpha=[%s], conv=%s, iter=%d, obj=%.6f\n",
      m, paste(round(fit$coefficients, 4), collapse=", "),
      fit$gmm_fit$opt$converged, fit$gmm_fit$opt$iterations, fit$gmm_fit$opt$objective))
}
cat("OK\n")
