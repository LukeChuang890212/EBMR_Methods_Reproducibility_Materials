res <- readRDS("Simulation_Results/EBMR_IPW_setting4-miss50-scenario9-3_13_n2000_replicate1000_test59.RDS")
mu <- res["mu_ipw",]
valid <- !is.na(mu)
cat("Valid:", sum(valid), "\n")
cat("mu_ipw summary:\n")
print(summary(mu[valid]))
cat("\nmu_ipw quantiles:\n")
print(quantile(mu[valid], c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1)))

cat("\nw.hat1 summary:\n")
print(summary(as.numeric(res["w.hat1", valid])))
cat("w.hat2 summary:\n")
print(summary(as.numeric(res["w.hat2", valid])))

cat("\nnu.hat1 summary:\n")
print(summary(as.numeric(res["nu.hat1", valid])))
cat("nu.hat2 summary:\n")
print(summary(as.numeric(res["nu.hat2", valid])))

# Identify outlier reps
mu_v <- mu[valid]
mu_med <- median(mu_v)
mu_iqr <- IQR(mu_v)
outlier_idx <- which(abs(mu_v - mu_med) > 3 * mu_iqr)
cat("\nOutlier reps (3*IQR):", length(outlier_idx), "\n")

if (length(outlier_idx) > 0) {
  cat("\nOutlier details:\n")
  for (i in head(outlier_idx, 20)) {
    cat(sprintf("  Rep %d: mu=%.4f, w1=%.4f, w2=%.4f, nu1=%.4f, nu2=%.4f\n",
        which(valid)[i], mu_v[i],
        res["w.hat1", which(valid)[i]], res["w.hat2", which(valid)[i]],
        res["nu.hat1", which(valid)[i]], res["nu.hat2", which(valid)[i]]))
  }
}

# Compare weights for outlier vs normal reps
normal_idx <- setdiff(seq_along(mu_v), outlier_idx)
cat("\n=== NORMAL reps (n=", length(normal_idx), ") ===\n")
cat("Mean w.hat1:", mean(res["w.hat1", which(valid)[normal_idx]]), "\n")
cat("Mean w.hat2:", mean(res["w.hat2", which(valid)[normal_idx]]), "\n")

cat("\n=== OUTLIER reps (n=", length(outlier_idx), ") ===\n")
if (length(outlier_idx) > 0) {
  cat("Mean w.hat1:", mean(res["w.hat1", which(valid)[outlier_idx]]), "\n")
  cat("Mean w.hat2:", mean(res["w.hat2", which(valid)[outlier_idx]]), "\n")
}

# Correlation between w.hat1 and extreme mu
cat("\nCorrelation(w.hat1, mu_ipw):", cor(res["w.hat1", valid], mu_v), "\n")
cat("Correlation(w.hat1, abs(mu - median)):", cor(res["w.hat1", valid], abs(mu_v - median(mu_v))), "\n")
