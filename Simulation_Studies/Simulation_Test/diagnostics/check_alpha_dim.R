setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
  devtools::load_all("../EPS")
})
W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
dat <- all_data[1:2000, ]
for (j in 1:3) {
  ps_sub <- list(formula.list=ps_spec[["formula.list"]][j], h_alpha.list=ps_spec[["h_alpha.list"]][j],
                 inv_link=ps_spec[["inv_link"]], outcome=ps_spec[["outcome"]])
  ebmr <- EPS[["new"]]("y", ps_sub, dat, W_func)
  gf <- ebmr[["ps_fit.list"]][[1]][["gmm_fit"]]
  cat(sprintf("Model %d: param_dim=%d, esteq_dim=%d (overid=%d)\n",
      j, length(gf[["estimates"]]), ncol(gf[["g.matrix"]]), ncol(gf[["g.matrix"]])-length(gf[["estimates"]])))
}
