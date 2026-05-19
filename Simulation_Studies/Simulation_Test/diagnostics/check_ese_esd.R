setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
  library(EBMRalgorithmFast4)
})

W_func  <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_nu_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1_u2=dat$u1*dat$u2)

run_check <- function(data_file, ps_spec_id, model_idx, setting, n_val, n_reps = 500, label = "") {
  if (!file.exists(data_file)) { cat("  MISSING:", data_file, "\n"); return(NULL) }
  all_data <- readRDS(data_file)
  ps_spec  <- get_ps_spec(ps_spec_id)
  ps_sub   <- list(formula.list = ps_spec[["formula.list"]][model_idx],
                   h_alpha.list = ps_spec[["h_alpha.list"]][model_idx],
                   inv_link     = ps_spec[["inv_link"]],
                   outcome      = ps_spec[["outcome"]])
  mu_true  <- get_mu_true(setting)

  mu_vals <- se_vals <- numeric(n_reps)
  hess_pd <- logical(n_reps)
  for (i in 1:n_reps) {
    set.seed(12345 + i)
    dat  <- all_data[((i-1)*n_val + 1):(i*n_val), ]
    ebmr <- EBMRAlgorithmFast4$new("y", ps_sub, dat, W_func)
    res  <- ebmr$EBMR_IPW(h_nu=h_nu_fn, type="HT", se.fit=TRUE)
    mu_vals[i] <- res$mu_ipw
    se_vals[i] <- res$se_ipw
    hess_pd[i] <- isTRUE(ebmr$ps_fit.list[[1]]$gmm_fit$opt$hessian_pd)
  }
  cat(sprintf("%-45s  Bias=%+.4f  ESD=%.4f  ASE=%.4f  ASE/ESD=%.3f  hess_pd=%d/%d\n",
      label,
      mean(mu_vals) - mu_true, sd(mu_vals), mean(se_vals),
      mean(se_vals)/sd(mu_vals),
      sum(hess_pd), n_reps))
}

n_reps <- 500
cat(sprintf("=== ESE/ESD check, %d reps each ===\n\n", n_reps))

# setting3 (misspecified B1 data, n=2000), scenario 9-3
run_check("Simulation_Data/Setting3.B1_n2000_replicate1000.RDS",
          "9-alt1", 1, "setting3", 2000, n_reps, "setting3-miss50-9-3_1 (M1 only)")
run_check("Simulation_Data/Setting3.B1_n2000_replicate1000.RDS",
          "9-alt1", 2, "setting3", 2000, n_reps, "setting3-miss50-9-3_2 (M2 only)")
run_check("Simulation_Data/Setting3.B1_n2000_replicate1000.RDS",
          "9-alt1", 3, "setting3", 2000, n_reps, "setting3-miss50-9-3_3 (M3 only)")

# setting4 (misspecified B1 data, n=2000), scenario 9-3
cat("\n")
run_check("Simulation_Data/Setting4.B1_n2000_replicate1000.RDS",
          "9-alt1", 1, "setting4", 2000, n_reps, "setting4-miss50-9-3_1 (M1 only)")
run_check("Simulation_Data/Setting4.B1_n2000_replicate1000.RDS",
          "9-alt1", 2, "setting4", 2000, n_reps, "setting4-miss50-9-3_2 (M2 only)")
run_check("Simulation_Data/Setting4.B1_n2000_replicate1000.RDS",
          "9-alt1", 3, "setting4", 2000, n_reps, "setting4-miss50-9-3_3 (M3 only)")
