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
n_val <- 2000
h_nu_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1_u2=dat$u1*dat$u2)

# Quick correctness check on rep 1
set.seed(12346)
dat1 <- all_data[1:n_val, ]
ebmr1 <- EBMRAlgorithmFast4$new("y", ps_spec_1, dat1, W_func)
res1  <- ebmr1$EBMR_IPW(h_nu=h_nu_fn, type="HT", se.fit=TRUE)
cat(sprintf("Correctness check: mu=%.4f se=%.4f alpha=[%s] hess_pd=%s\n",
    res1$mu_ipw, res1$se_ipw,
    paste(round(ebmr1$ps_fit.list[[1]]$coefficients, 3), collapse=", "),
    as.character(ebmr1$ps_fit.list[[1]]$gmm_fit$opt$hessian_pd)))

# Timing: 20 reps
times <- numeric(20)
for (i in 1:20) {
  set.seed(12345 + i)
  dat <- all_data[((i-1)*n_val+1):(i*n_val), ]
  t0 <- proc.time()["elapsed"]
  ebmr <- EBMRAlgorithmFast4$new("y", ps_spec_1, dat, W_func)
  ebmr$EBMR_IPW(h_nu=h_nu_fn, type="HT", se.fit=TRUE)
  times[i] <- proc.time()["elapsed"] - t0
}
cat(sprintf("Timing: mean=%.3f  min=%.3f  max=%.3f  total=%.1f sec (n=20 reps)\n",
    mean(times), min(times), max(times), sum(times)))
cat(sprintf("Projected wall time (1000 reps, %d cores): %.1f min\n",
    max(1, parallel::detectCores()-2),
    1000 * mean(times) / max(1, parallel::detectCores()-2) / 60))
