setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
source("Simulation.r")
library(EBMRalgorithmFast4)

# Setting4, scenario 9-3_2 = model 2 only from ps_spec "9-alt1"
ps_spec <- get_ps_spec("9-alt1")
ps_spec_2 <- list(
  formula.list = ps_spec[["formula.list"]][2],
  h_alpha.list = ps_spec[["h_alpha.list"]][2],
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]]
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
n_val  <- 2000
mu_true <- get_mu_true("setting4")
n_reps <- 100

cat(sprintf("Setting4, scenario 9-3_2 (model 2 only), n=%d, n_reps=%d\n", n_val, n_reps))
cat(sprintf("mu_true = %.6f\n\n", mu_true))

mu_vals    <- numeric(n_reps)
se_vals    <- numeric(n_reps)
alpha_norms <- numeric(n_reps)
hess_pd    <- logical(n_reps)
alpha_vals <- matrix(NA, n_reps, 4)

t0 <- Sys.time()
for (i in 1:n_reps) {
  set.seed(12345 + i)
  dat <- setting4.B1(n_val)

  ebmr   <- EBMRAlgorithmFast4$new("y", ps_spec_2, dat, W_func)
  result <- ebmr$EBMR_IPW(
    h_nu = function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1_u2=dat$u1*dat$u2),
    type = "HT", se.fit = TRUE
  )

  mu_vals[i]     <- result$mu_ipw
  se_vals[i]     <- result$se_ipw
  ps_fit         <- ebmr$ps_fit.list[[1]]
  alpha_vals[i,] <- ps_fit$coefficients
  alpha_norms[i] <- sqrt(sum(ps_fit$coefficients^2))
  hess_pd[i]     <- ps_fit$gmm_fit$opt$hessian_pd

  if (i %% 20 == 0) cat(sprintf("Rep %d: mu=%.4f se=%.4f alpha_norm=%.2f hess_pd=%s\n",
      i, mu_vals[i], se_vals[i], alpha_norms[i], hess_pd[i]))
}
elapsed <- as.numeric(difftime(Sys.time(), t0, units="secs"))
cat(sprintf("\nTime: %.1f sec (%.2f sec/rep)\n", elapsed, elapsed/n_reps))

cat("\n=== Results (current package) ===\n")
cat(sprintf("Bias:      %.4f\n", mean(mu_vals) - mu_true))
cat(sprintf("ESD:       %.4f\n", sd(mu_vals)))
cat(sprintf("ASE:       %.4f\n", mean(se_vals)))
cat(sprintf("ASE/ESD:   %.4f\n", mean(se_vals)/sd(mu_vals)))

# Classify extreme alpha
extreme <- abs(alpha_vals[,1]) > 10 | abs(alpha_vals[,2]) > 10
cat(sprintf("\nExtreme alpha reps: %d / %d\n", sum(extreme, na.rm=TRUE), n_reps))
cat(sprintf("hessian_pd=TRUE: %d, FALSE: %d\n", sum(hess_pd, na.rm=TRUE), sum(!hess_pd, na.rm=TRUE)))

cat("\nalpha1 quantiles:", round(quantile(alpha_vals[,1], c(0,.1,.5,.9,1), na.rm=TRUE), 2), "\n")
cat("alpha2 quantiles:", round(quantile(alpha_vals[,2], c(0,.1,.5,.9,1), na.rm=TRUE), 2), "\n")

# Compare old file
old_file <- "Simulation_Results/EBMR_IPW_setting4-miss50-scenario9-3_2_n2000_replicate1000_test59.RDS"
if (file.exists(old_file)) {
  old <- readRDS(old_file)
  cat("\n=== Old results (test59) ===\n")
  cat(sprintf("Bias:    %.4f\n", mean(old["mu_ipw",]) - mu_true))
  cat(sprintf("ESD:     %.4f\n", sd(old["mu_ipw",])))
  cat(sprintf("ASE:     %.4f\n", mean(old["se_ipw",])))
  cat(sprintf("ASE/ESD: %.4f\n", mean(old["se_ipw",])/sd(old["mu_ipw",])))
}
