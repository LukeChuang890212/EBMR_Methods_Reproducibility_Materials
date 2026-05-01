setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

files <- c(
  "Simulation_Results/EBMR_IPW_setting3-miss50-scenario9-3_1_n2000_replicate1000_test57.RDS",
  "Simulation_Results/EBMR_IPW_setting3-miss50-scenario9-3_2_n2000_replicate1000_test57.RDS",
  "Simulation_Results/EBMR_IPW_setting3-miss50-scenario9-3_3_n2000_replicate1000_test57.RDS"
)

for (f in files) {
  cat(sprintf("\n=== %s ===\n", basename(f)))
  res <- readRDS(f)
  cat(sprintf("Dimensions: %d x %d\n", nrow(res), ncol(res)))
  cat("Row names:", rownames(res), "\n")

  mu <- res["mu_ipw", ]
  valid <- !is.na(mu)
  mu_v <- mu[valid]
  cat(sprintf("Valid: %d/%d\n", sum(valid), length(mu)))
  cat(sprintf("mu_ipw: mean=%.4f, sd=%.4f, min=%.4f, max=%.4f\n",
              mean(mu_v), sd(mu_v), min(mu_v), max(mu_v)))

  q1 <- quantile(mu_v, 0.25)
  q3 <- quantile(mu_v, 0.75)
  iqr <- q3 - q1
  outliers <- mu_v < (q1 - 3*iqr) | mu_v > (q3 + 3*iqr)
  cat(sprintf("Outliers (3*IQR): %d/%d (%.1f%%)\n", sum(outliers), length(mu_v), 100*mean(outliers)))

  if (sum(outliers) > 0) {
    cat("Outlier values:", head(round(sort(mu_v[outliers]), 4), 20), "\n")
  }

  cat(sprintf("Quantiles: 1%%=%.4f, 5%%=%.4f, 95%%=%.4f, 99%%=%.4f\n",
      quantile(mu_v, 0.01), quantile(mu_v, 0.05), quantile(mu_v, 0.95), quantile(mu_v, 0.99)))
}
cat("\nDone!\n")
