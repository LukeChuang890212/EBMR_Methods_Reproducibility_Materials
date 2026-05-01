setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
source("Simulation.r")
source("Data_Generation.r")
source("config/scenarios.R")

f <- "Simulation_Results/EBMR_IPW_setting3-miss50-scenario9-3_2_n2000_replicate1000_test59.RDS"
res <- readRDS(f)

cleaned <- clean_sim_result(res, verbose = TRUE)
sim <- cleaned$result

mu_v <- sim[1, ]
se_v <- sim[3, ]

mu.true <- get_mu_true("setting3")

bias <- mean(mu_v) - mu.true
esd <- sd(mu_v)
ese <- mean(se_v)
ci_lower <- mu_v - 1.96 * se_v
ci_upper <- mu_v + 1.96 * se_v
cp <- mean(mu.true >= ci_lower & mu.true <= ci_upper)

cat(sprintf("mu.true = %.4f\n", mu.true))
cat(sprintf("Bias    = %.4f\n", bias))
cat(sprintf("ESD     = %.4f\n", esd))
cat(sprintf("ESE     = %.4f\n", ese))
cat(sprintf("ESE/ESD = %.2f\n", ese / esd))
cat(sprintf("CP      = %.4f\n", cp))
