setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

f <- "Simulation_Results/EBMR_IPW_setting3-miss50-scenario9-3_2_n2000_replicate1000_test57.RDS"
res <- readRDS(f)

mu <- res["mu_ipw", ]
q1 <- quantile(mu, 0.25)
q3 <- quantile(mu, 0.75)
iqr <- q3 - q1
is_outlier <- mu < (q1 - 3*iqr) | mu > (q3 + 3*iqr)

cat(sprintf("Outliers: %d/1000\n\n", sum(is_outlier)))

# Check alpha coefficients for outlier vs non-outlier reps
alpha_rows <- grep("alpha.hat", rownames(res))
cat("Alpha coefficient rows:", rownames(res)[alpha_rows], "\n\n")

alpha_mat <- res[alpha_rows, ]

# Max |alpha| for each replicate
max_alpha <- apply(abs(alpha_mat), 2, max)

cat("--- Max |alpha| stats ---\n")
cat(sprintf("  Outlier reps:     mean=%.2f, median=%.2f, max=%.2f\n",
            mean(max_alpha[is_outlier]), median(max_alpha[is_outlier]), max(max_alpha[is_outlier])))
cat(sprintf("  Non-outlier reps: mean=%.2f, median=%.2f, max=%.2f\n",
            mean(max_alpha[!is_outlier]), median(max_alpha[!is_outlier]), max(max_alpha[!is_outlier])))

cat(sprintf("\n  P(max|alpha|>=4.9 | outlier) = %.1f%%\n", 100*mean(max_alpha[is_outlier] >= 4.9)))
cat(sprintf("  P(max|alpha|>=4.9 | no outlier) = %.1f%%\n", 100*mean(max_alpha[!is_outlier] >= 4.9)))

# Show alpha values for first few outliers
outlier_idx <- which(is_outlier)
cat(sprintf("\nFirst 5 outlier reps (alpha values):\n"))
for (k in head(outlier_idx, 5)) {
  cat(sprintf("  Rep %d: mu=%.4f, alpha=[%s], max|a|=%.2f\n",
              k, mu[k],
              paste(round(alpha_mat[, k], 3), collapse=", "),
              max(abs(alpha_mat[, k]))))
}

# nu and w for outlier reps
cat(sprintf("\nFirst 5 outlier reps (nu, w):\n"))
for (k in head(outlier_idx, 5)) {
  cat(sprintf("  Rep %d: nu=%.4f, w=%.4f\n", k, res["nu.hat", k], res["w.hat", k]))
}

cat("\nDone!\n")
