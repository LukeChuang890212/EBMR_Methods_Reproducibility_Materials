setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

f <- "Simulation_Results/EBMR_IPW_setting3-miss50-scenario9-3_2_n2000_replicate1000_test59.RDS"
res <- readRDS(f)
mu <- res[1,]; se <- res[3,]
valid <- !is.na(mu) & !is.na(se)
mu_v <- mu[valid]; se_v <- se[valid]

cat("=== Current approach: outlier on mu OR se ===\n")
q1 <- quantile(mu_v, 0.25); q3 <- quantile(mu_v, 0.75); iqr <- q3 - q1
q1s <- quantile(se_v, 0.25); q3s <- quantile(se_v, 0.75); iqrs <- q3s - q1s
is_out_mu <- mu_v < q1-3*iqr | mu_v > q3+3*iqr
is_out_se <- se_v < q1s-3*iqrs | se_v > q3s+3*iqrs
is_out_both <- is_out_mu | is_out_se

cat(sprintf("  mu outliers: %d, se outliers: %d, union: %d\n",
            sum(is_out_mu), sum(is_out_se), sum(is_out_both)))
mu_c <- mu_v[!is_out_both]; se_c <- se_v[!is_out_both]
cat(sprintf("  n=%d, ESD=%.4f, ESE=%.4f, ratio=%.2f\n",
            length(mu_c), sd(mu_c), mean(se_c), mean(se_c)/sd(mu_c)))

cat("\n=== Proposed: outlier on mu ONLY ===\n")
mu_c2 <- mu_v[!is_out_mu]; se_c2 <- se_v[!is_out_mu]
cat(sprintf("  mu outliers removed: %d\n", sum(is_out_mu)))
cat(sprintf("  n=%d, ESD=%.4f, ESE=%.4f, ratio=%.2f\n",
            length(mu_c2), sd(mu_c2), mean(se_c2), mean(se_c2)/sd(mu_c2)))

# Try different IQR multipliers for mu-only
cat("\n=== mu-only outlier with different IQR multipliers ===\n")
for (k in c(1.5, 2, 2.5, 3, 4, 5)) {
  out_k <- mu_v < q1-k*iqr | mu_v > q3+k*iqr
  mu_k <- mu_v[!out_k]; se_k <- se_v[!out_k]
  cat(sprintf("  k=%.1f: removed=%d (%.1f%%), ESD=%.4f, ESE=%.4f, ratio=%.2f\n",
              k, sum(out_k), 100*mean(out_k), sd(mu_k), mean(se_k), mean(se_k)/sd(mu_k)))
}

# Also check: how many se-only outliers remain when we only trim mu?
se_only_out <- is_out_se & !is_out_mu
cat(sprintf("\n=== se-only outliers (not mu outliers): %d ===\n", sum(se_only_out)))
if (sum(se_only_out) > 0) {
  cat("  These have normal mu but inflated se:\n")
  idx <- which(se_only_out)
  for (i in idx[1:min(10, length(idx))]) {
    cat(sprintf("    mu=%.4f, se=%.4f\n", mu_v[i], se_v[i]))
  }
}

cat("\nDone!\n")
