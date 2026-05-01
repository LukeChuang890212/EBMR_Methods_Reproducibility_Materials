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
alpha.true <- misspecified_model_alpha.true.list$setting3$miss50$setting3_1_mild_1000
ps_model.true <- function(dat, alpha.true) {
  1 / (1 + exp(cbind(rep(1, nrow(dat)), dat$y, dat$u1, dat$u2) %*% alpha.true))
}

n_reps <- 20

# Current package (jacobian of conv_grad)
t0 <- Sys.time()
for (i in 1:n_reps) {
  dat <- all_data[((i-1)*n + 1):(i*n), ]
  ebmr <- EBMRAlgorithmFast4$new("y", ps_spec_2, dat, W_func)
  result <- ebmr$EBMR_IPW(
    h_nu = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2, u1_u2 = dat$u1*dat$u2),
    true_ps = ps_model.true(dat, alpha.true),
    type = "HT", se.fit = TRUE
  )
}
t_new <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("jacobian(conv_grad): %.1f sec for %d reps (%.2f sec/rep)\n", t_new, n_reps, t_new/n_reps))
