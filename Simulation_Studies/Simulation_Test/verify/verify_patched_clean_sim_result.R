setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({ source("Basic_setup.r"); source("Data_Generation.r"); source("Simulation.r") })
x <- readRDS("Simulation_Results/EBMR_IPW_setting1-miss50-scenario9-3_3_n2000_replicate1000_test65.RDS")
mu_true <- get_mu_true("setting1")
c <- clean_sim_result(x, multiplier = 2, max_pct = 0.01, verbose = FALSE)
sc <- c$result
mu_v <- sc["mu_ipw", ]; se_v <- sc["se_ipw", ]
ci_lo <- mu_v - 1.96*se_v; ci_hi <- mu_v + 1.96*se_v
cat(sprintf("Patched clean_sim_result on test65 file:\n"))
cat(sprintf("  n_kept=%d  n_outliers=%d\n", c$n_successful, c$n_outliers))
cat(sprintf("  Bias=%+.4f  ESD=%.4f  ESE(mean)=%.4f  ESE(median)=%.4f  ESE/ESD=%.3f  CP=%.3f\n",
            mean(mu_v) - mu_true, sd(mu_v), mean(se_v), median(se_v),
            mean(se_v)/sd(mu_v), mean(ci_lo <= mu_true & mu_true <= ci_hi)))
cat("\nExpected:  n=990  Bias=-0.2795  ESD=0.1516  ESE=0.0945  ratio=0.624  CP=0.139\n")
