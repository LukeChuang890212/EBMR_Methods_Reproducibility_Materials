## Diagnose bad ESE/ESD for setting3-miss50-scenario7, models 1&2
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

res <- readRDS("Simulation_Results/EBMR_IPW_setting3-miss50-scenario7-1_2_n2000_replicate1000_test61.RDS")
cat("Dimensions:", dim(res), "\n")
cat("Row names:", rownames(res), "\n\n")

mu_ipw <- res[1,]
se_ipw <- res[3,]
mu_true <- res[2,]

valid <- !is.na(mu_ipw) & !is.na(se_ipw)
cat("Valid reps:", sum(valid), "/", ncol(res), "\n\n")

cat("mu_ipw distribution:\n")
print(summary(mu_ipw[valid]))
cat("\nse_ipw distribution:\n")
print(summary(se_ipw[valid]))
cat("\nmu_ipw.true distribution:\n")
print(summary(mu_true[valid]))

# Check for extreme SE values
cat("\n\nExtreme SE (>1):", sum(se_ipw[valid] > 1), "\n")
cat("Extreme SE (>0.5):", sum(se_ipw[valid] > 0.5), "\n")
cat("SE < 0.1:", sum(se_ipw[valid] < 0.1), "\n")
cat("SE quantiles:\n")
print(quantile(se_ipw[valid], c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1)))

# Check for extreme mu values
cat("\nmu_ipw quantiles:\n")
print(quantile(mu_ipw[valid], c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1)))

# Correlation between mu and se
cat("\ncor(mu_ipw, se_ipw):", cor(mu_ipw[valid], se_ipw[valid]), "\n")

# Compare median vs mean SE
cat("\nMean SE:", mean(se_ipw[valid]), "\n")
cat("Median SE:", median(se_ipw[valid]), "\n")
cat("ESD:", sd(mu_ipw[valid]), "\n")
cat("Median SE / ESD:", median(se_ipw[valid]) / sd(mu_ipw[valid]), "\n")
