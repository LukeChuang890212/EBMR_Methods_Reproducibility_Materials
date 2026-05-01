setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
  library(EBMRalgorithmFast4)
})

ps_spec <- get_ps_spec("9-alt1")
ps_spec_1 <- list(
  formula.list = ps_spec[["formula.list"]][1],
  h_alpha.list = ps_spec[["h_alpha.list"]][1],
  inv_link     = ps_spec[["inv_link"]],
  outcome      = ps_spec[["outcome"]]
)
W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
all_data <- readRDS("Simulation_Data/Setting3.B1_n2000_replicate1000.RDS")
n_val  <- 2000
mu_true <- get_mu_true("setting3")
h_nu_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1_u2=dat$u1*dat$u2)

n_reps <- 100
mu_vals <- se_vals <- numeric(n_reps)
t0 <- proc.time()["elapsed"]
for (i in 1:n_reps) {
  set.seed(12345 + i)
  dat <- all_data[((i-1)*n_val+1):(i*n_val), ]
  ebmr <- EBMRAlgorithmFast4$new("y", ps_spec_1, dat, W_func)
  res  <- ebmr$EBMR_IPW(h_nu=h_nu_fn, type="HT", se.fit=TRUE)
  mu_vals[i] <- res$mu_ipw
  se_vals[i] <- res$se_ipw
}
elapsed <- proc.time()["elapsed"] - t0

cat(sprintf("Setting3-miss50-scenario9-3_1, n=2000, %d reps\n", n_reps))
cat(sprintf("mu_true = %.4f\n", mu_true))
cat(sprintf("Bias:    %.4f\n", mean(mu_vals) - mu_true))
cat(sprintf("ESD:     %.4f\n", sd(mu_vals)))
cat(sprintf("ASE:     %.4f\n", mean(se_vals)))
cat(sprintf("ASE/ESD: %.4f\n", mean(se_vals)/sd(mu_vals)))
cat(sprintf("Time:    %.1f sec total, %.3f sec/rep\n", elapsed, elapsed/n_reps))
cat(sprintf("Projected (1000 reps, %d cores): %.1f min\n",
    max(1, parallel::detectCores()-2),
    1000 * (elapsed/n_reps) / max(1, parallel::detectCores()-2) / 60))
