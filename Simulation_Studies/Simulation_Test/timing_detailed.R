setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
source("Simulation.r")
library(EBMRalgorithmFast4)

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

# Time 10 reps with detailed breakdown
n_reps <- 10
cat(sprintf("Timing %d reps with detailed breakdown\n\n", n_reps))

for (i in 1:n_reps) {
  dat <- all_data[((i-1)*n + 1):(i*n), ]

  t0 <- Sys.time()
  ebmr <- EBMRAlgorithmFast4$new("y", ps_spec_2, dat, W_func)
  t_fit <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  t1 <- Sys.time()
  result <- ebmr$EBMR_IPW(
    h_nu = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2, u1_u2 = dat$u1*dat$u2),
    true_ps = ps_model.true(dat, alpha.true),
    type = "HT", se.fit = TRUE
  )
  t_ipw <- as.numeric(difftime(Sys.time(), t1, units = "secs"))

  ps_fit <- ebmr$ps_fit.list[[1]]
  opt <- ps_fit$gmm_fit$opt

  cat(sprintf("Rep %2d: fit=%.3fs ipw=%.3fs total=%.3fs | iter=%d conv=%s hpd=%s grad=%.2e obj=%.2e\n",
      i, t_fit, t_ipw, t_fit + t_ipw,
      opt$iterations, opt$converged, opt$hessian_pd,
      opt$final_grad_norm, opt$objective))
}
