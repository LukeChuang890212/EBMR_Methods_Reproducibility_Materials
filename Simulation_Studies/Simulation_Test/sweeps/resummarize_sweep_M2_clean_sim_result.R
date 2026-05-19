## Re-summarize the M2 sweep using project's clean_sim_result()
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({ source("Simulation.r"); source("Data_Generation.r") })

dat <- readRDS("Simulation_Test/sweep_cnr_s1_miss50_n2000_M2_results.RDS")
mu_true <- get_mu_true("setting1")

cat(sprintf("\n=== M2 sweep re-summarized via clean_sim_result() (mu_true=%.4f) ===\n", mu_true))
cat(sprintf("%-22s | %4s %4s | %7s %7s %7s | %s\n",
            "tag", "NA", "Out", "Bias", "ESD", "ESE", "ESE/ESD CP"))

for (tag in names(dat)) {
  res <- dat[[tag]]$res  # 500 x 7 matrix
  # Build sim_result-style matrix: rows = mu_ipw, mu_ipw.true, se_ipw, se_ipw.true, alpha.hat1..4
  # We don't have true comparators; duplicate mu_ipw/se_ipw into the .true rows so the
  # NA-mask in clean_sim_result (which checks all 4 first rows) doesn't drop anything.
  sim_mat <- rbind(
    mu_ipw      = res[, "mu_ipw"],
    mu_ipw.true = res[, "mu_ipw"],
    se_ipw      = res[, "se_ipw"],
    se_ipw.true = res[, "se_ipw"],
    alpha.hat1  = res[, "a1"],
    alpha.hat2  = res[, "a2"],
    alpha.hat3  = res[, "a3"],
    alpha.hat4  = res[, "a4"]
  )
  cleaned <- clean_sim_result(sim_mat, multiplier = 2, max_pct = 0.01, verbose = FALSE)
  sc <- cleaned$result
  mu_v <- sc["mu_ipw", ]; se_v <- sc["se_ipw", ]
  bias <- mean(mu_v) - mu_true
  esd  <- sd(mu_v)
  ese  <- mean(se_v)
  ratio <- ese / esd
  ci_lo <- mu_v - 1.96*se_v; ci_hi <- mu_v + 1.96*se_v
  cp <- mean(ci_lo <= mu_true & mu_true <= ci_hi)
  cat(sprintf("%-22s | %4d %4d | %+.4f %.4f %.4f | %.3f  %.3f\n",
              tag, cleaned$n_na, cleaned$n_outliers,
              bias, esd, ese, ratio, cp))
}
cat("\n")
