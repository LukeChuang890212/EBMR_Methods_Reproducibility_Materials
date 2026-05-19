res <- readRDS("Simulation_Results/EBMR_IPW_setting4-miss50-scenario9-3_13_n2000_replicate1000_test59.RDS")
mu <- res["mu_ipw",]
se <- res["se_ipw",]
valid <- !is.na(mu)
mu_v <- mu[valid]
se_v <- se[valid]

# Need mu.true for setting4
source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
source("Simulation.r")
mu.true <- get_mu_true("setting4")

cat(sprintf("mu.true = %.4f\n", mu.true))
cat(sprintf("\n=== BEFORE OUTLIER REMOVAL (n=%d) ===\n", length(mu_v)))
cat(sprintf("Bias = %.4f, ESD = %.4f, ESE = %.4f, ESE/ESD = %.2f\n",
            mean(mu_v) - mu.true, sd(mu_v), mean(se_v), mean(se_v)/sd(mu_v)))
ci_l <- mu_v - 1.96 * se_v; ci_u <- mu_v + 1.96 * se_v
cat(sprintf("CP = %.4f\n", mean(mu.true >= ci_l & mu.true <= ci_u)))
cat(sprintf("Skewness = %.4f, Excess kurtosis = %.4f\n",
            mean((mu_v - mean(mu_v))^3) / sd(mu_v)^3,
            mean((mu_v - mean(mu_v))^4) / sd(mu_v)^4 - 3))

# Outlier removal
source("Simulation.r")
results_full <- matrix(NA, 4, length(mu))
rownames(results_full) <- c("mu_ipw", "mu_ipw.true", "se_ipw", "se_ipw.true")
results_full[1,] <- mu
results_full[2,] <- ifelse(is.na(mu), NA, 0)
results_full[3,] <- se
results_full[4,] <- ifelse(is.na(se), NA, 0)

for (mult in c(1.5, 3)) {
  cleaned <- clean_sim_result(results_full, multiplier = mult, verbose = TRUE)
  sim <- cleaned$result
  mu_c <- sim[1, ]; se_c <- sim[3, ]
  bias <- mean(mu_c) - mu.true
  esd <- sd(mu_c)
  ese <- mean(se_c)
  ci_lower <- mu_c - 1.96 * se_c
  ci_upper <- mu_c + 1.96 * se_c
  cp <- mean(mu.true >= ci_lower & mu.true <= ci_upper)
  cat(sprintf("\n=== AFTER OUTLIER REMOVAL (mult=%.1f, n=%d) ===\n", mult, length(mu_c)))
  cat(sprintf("Bias = %.4f, ESD = %.4f, ESE = %.4f, ESE/ESD = %.2f, CP = %.4f\n",
              bias, esd, ese, ese/esd, cp))
}
