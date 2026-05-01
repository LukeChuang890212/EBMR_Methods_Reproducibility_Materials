setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
source("Simulation.r")
library(EBMRalgorithmFast4)
library(numDeriv)

ps_spec <- get_ps_spec("9-alt1")
ps_spec_2 <- list(
  formula.list = ps_spec[["formula.list"]][2],
  h_alpha.list = ps_spec[["h_alpha.list"]][2],
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]]
)
W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
n <- 2000

all_data <- readRDS("Simulation_Data/Setting3.B1_n2000_replicate1000.RDS")

# Manually replicate GMM for a slow rep to trace iteration-by-iteration
# First find a slow rep
cat("=== Finding slow reps ===\n")
for (i in 1:20) {
  dat <- all_data[((i-1)*n + 1):(i*n), ]
  t0 <- Sys.time()
  ebmr <- EBMRAlgorithmFast4$new("y", ps_spec_2, dat, W_func)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  opt <- ebmr$ps_fit.list[[1]]$gmm_fit$opt
  cat(sprintf("Rep %2d: %.3fs  iter=%d  conv=%s  hpd=%s  obj=%.2e  grad=%.2e\n",
      i, elapsed, opt$iterations, opt$converged, opt$hessian_pd,
      opt$objective, opt$final_grad_norm))
}
