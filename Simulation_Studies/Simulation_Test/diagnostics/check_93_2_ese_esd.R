setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({source("Basic_setup.r"); source("Data_Generation.r"); source("config/scenarios.R"); source("Simulation.r")})

f <- "Simulation_Results/EBMR_IPW_setting3-miss50-scenario9-3_2_n2000_replicate1000_test59.RDS"
sim <- readRDS(f)
cat("dim:", dim(sim), "\n")
mu_true <- get_mu_true("setting3")
cat("mu_true:", mu_true, "\n")
mu_vals <- sim[1,]
se_vals <- sim[3,]
cat("n_reps:", sum(!is.na(mu_vals)), "\n")
cat("Bias:", round(mean(mu_vals, na.rm=TRUE) - mu_true, 4), "\n")
cat("ESD:", round(sd(mu_vals, na.rm=TRUE), 4), "\n")
cat("ESE:", round(mean(se_vals, na.rm=TRUE), 4), "\n")
cat("ESE/ESD:", round(mean(se_vals, na.rm=TRUE)/sd(mu_vals, na.rm=TRUE), 3), "\n")

cleaned <- clean_sim_result(sim, multiplier=3, verbose=TRUE)
sc <- cleaned$result
mu_c <- sc[1,]; se_c <- sc[3,]
cat("\nAfter cleaning:\n")
cat("n_clean:", sum(!is.na(mu_c)), "\n")
cat("Bias:", round(mean(mu_c, na.rm=TRUE) - mu_true, 4), "\n")
cat("ESD:", round(sd(mu_c, na.rm=TRUE), 4), "\n")
cat("ESE:", round(mean(se_c, na.rm=TRUE), 4), "\n")
cat("ESE/ESD:", round(mean(se_c, na.rm=TRUE)/sd(mu_c, na.rm=TRUE), 3), "\n")
