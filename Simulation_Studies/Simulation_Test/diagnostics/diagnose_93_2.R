setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

f <- "Simulation_Results/EBMR_IPW_setting3-miss50-scenario9-3_2_n2000_replicate1000_test59.RDS"
res <- readRDS(f)

cat("Dimensions:", dim(res), "\n")
cat("Rownames:", paste(rownames(res), collapse=", "), "\n\n")

mu <- res[1, ]; se <- res[3, ]
mu_true <- res[2, ]; se_true <- res[4, ]

valid <- !is.na(mu) & !is.na(se)
cat("Valid replicates:", sum(valid), "/", length(mu), "\n\n")

mu_v <- mu[valid]; se_v <- se[valid]

cat("=== mu_ipw ===\n")
cat(sprintf("  mean=%.4f, sd(ESD)=%.4f, range=[%.4f, %.4f]\n",
            mean(mu_v), sd(mu_v), min(mu_v), max(mu_v)))

cat("=== se_ipw (ESE) ===\n")
cat(sprintf("  mean(ESE)=%.4f, sd=%.4f, range=[%.4f, %.4f]\n",
            mean(se_v), sd(se_v), min(se_v), max(se_v)))

cat(sprintf("\nESD=%.4f, mean(ESE)=%.4f, ratio=%.2f\n",
            sd(mu_v), mean(se_v), mean(se_v)/sd(mu_v)))

# Check for outliers
q1 <- quantile(mu_v, 0.25); q3 <- quantile(mu_v, 0.75); iqr <- q3 - q1
is_out_mu <- mu_v < (q1 - 3*iqr) | mu_v > (q3 + 3*iqr)

q1s <- quantile(se_v, 0.25); q3s <- quantile(se_v, 0.75); iqrs <- q3s - q1s
is_out_se <- se_v < (q1s - 3*iqrs) | se_v > (q3s + 3*iqrs)

is_out <- is_out_mu | is_out_se
cat(sprintf("\nOutliers (3*IQR): %d mu, %d se, %d total\n",
            sum(is_out_mu), sum(is_out_se), sum(is_out)))

# After removing outliers
if (sum(!is_out) > 0) {
  mu_clean <- mu_v[!is_out]; se_clean <- se_v[!is_out]
  cat(sprintf("\nAfter removing outliers:\n"))
  cat(sprintf("  ESD=%.4f, mean(ESE)=%.4f, ratio=%.2f\n",
              sd(mu_clean), mean(se_clean), mean(se_clean)/sd(mu_clean)))
}

# Check weights
if ("w.hat1" %in% rownames(res)) {
  w1 <- res["w.hat1", valid]
  cat(sprintf("\nw.hat1: mean=%.4f, range=[%.4f, %.4f]\n",
              mean(w1), min(w1), max(w1)))
  cat(sprintf("  w1 < 0.5: %d (%.1f%%)\n", sum(w1 < 0.5), 100*mean(w1 < 0.5)))
}

# Check nu
if ("nu.hat1" %in% rownames(res)) {
  nu1 <- res["nu.hat1", valid]
  cat(sprintf("\nnu.hat1: mean=%.4f, sd=%.4f, range=[%.4f, %.4f]\n",
              mean(nu1), sd(nu1), min(nu1), max(nu1)))
}

# Show worst outliers
cat("\n=== Worst mu outliers ===\n")
worst_mu <- order(abs(mu_v - mean(mu_v)), decreasing=TRUE)[1:10]
for (idx in worst_mu) {
  rep_i <- which(valid)[idx]
  cat(sprintf("  rep %d: mu=%.4f, se=%.4f\n", rep_i, mu_v[idx], se_v[idx]))
}

cat("\n=== Worst se outliers ===\n")
worst_se <- order(abs(se_v - mean(se_v)), decreasing=TRUE)[1:10]
for (idx in worst_se) {
  rep_i <- which(valid)[idx]
  cat(sprintf("  rep %d: mu=%.4f, se=%.4f\n", rep_i, mu_v[idx], se_v[idx]))
}

cat("\nDone!\n")
