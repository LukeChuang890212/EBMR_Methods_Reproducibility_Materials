setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

f <- "Simulation_Results/EBMR_IPW_setting3-miss50-scenario9-3_2_n2000_replicate1000_test59.RDS"
res <- readRDS(f)
mu <- res[1,]; se <- res[3,]
valid <- !is.na(mu) & !is.na(se)
mu_v <- mu[valid]; se_v <- se[valid]

# Remove 3*IQR outliers
q1 <- quantile(mu_v, 0.25); q3 <- quantile(mu_v, 0.75); iqr <- q3 - q1
q1s <- quantile(se_v, 0.25); q3s <- quantile(se_v, 0.75); iqrs <- q3s - q1s
is_out <- (mu_v < q1-3*iqr | mu_v > q3+3*iqr) | (se_v < q1s-3*iqrs | se_v > q3s+3*iqrs)

cat("=== Overall ===\n")
cat(sprintf("  ESD=%.4f, ESE=%.4f, ratio=%.2f\n", sd(mu_v), mean(se_v), mean(se_v)/sd(mu_v)))

cat("\n=== After 3*IQR removal ===\n")
mu_c <- mu_v[!is_out]; se_c <- se_v[!is_out]
cat(sprintf("  n=%d, ESD=%.4f, ESE=%.4f, ratio=%.2f\n",
            length(mu_c), sd(mu_c), mean(se_c), mean(se_c)/sd(mu_c)))

# Key question: is the SE too small for the TYPICAL rep?
# Look at the distribution of se_v for non-outlier reps
cat(sprintf("\n=== Non-outlier SE distribution ===\n"))
cat(sprintf("  mean=%.4f, median=%.4f, sd=%.4f\n", mean(se_c), median(se_c), sd(se_c)))
cat(sprintf("  quantiles: 5%%=%.4f, 25%%=%.4f, 50%%=%.4f, 75%%=%.4f, 95%%=%.4f\n",
            quantile(se_c, 0.05), quantile(se_c, 0.25), quantile(se_c, 0.5),
            quantile(se_c, 0.75), quantile(se_c, 0.95)))

# Check: what is the ESD if we progressively remove more extreme reps?
cat("\n=== ESD as we trim tails ===\n")
for (pct in c(100, 99, 98, 95, 90, 80)) {
  lo <- quantile(mu_v, (100-pct)/200)
  hi <- quantile(mu_v, 1 - (100-pct)/200)
  keep <- mu_v >= lo & mu_v <= hi
  cat(sprintf("  Keep %3d%%: n=%d, ESD=%.4f, ESE=%.4f, ratio=%.2f\n",
              pct, sum(keep), sd(mu_v[keep]), mean(se_v[keep]), mean(se_v[keep])/sd(mu_v[keep])))
}

# Compare with model 1 (correct specification)
cat("\n=== For reference: model 1 (correct) ===\n")
f1 <- "Simulation_Results/EBMR_IPW_setting3-miss50-scenario9-3_1_n2000_replicate1000_test59.RDS"
if (file.exists(f1)) {
  res1 <- readRDS(f1)
  mu1 <- res1[1,]; se1 <- res1[3,]
  v1 <- !is.na(mu1) & !is.na(se1)
  cat(sprintf("  ESD=%.4f, ESE=%.4f, ratio=%.2f\n", sd(mu1[v1]), mean(se1[v1]), mean(se1[v1])/sd(mu1[v1])))
  for (pct in c(100, 99, 98, 95, 90)) {
    lo <- quantile(mu1[v1], (100-pct)/200)
    hi <- quantile(mu1[v1], 1 - (100-pct)/200)
    keep <- mu1[v1] >= lo & mu1[v1] <= hi
    cat(sprintf("  Keep %3d%%: n=%d, ESD=%.4f, ESE=%.4f, ratio=%.2f\n",
                pct, sum(keep), sd(mu1[v1][keep]), mean(se1[v1][keep]), mean(se1[v1][keep])/sd(mu1[v1][keep])))
  }
}

# Check: what fraction of reps have se < ESD_clean?
# i.e., the SE says the variability is X, but the actual ESD is larger
esd_clean <- sd(mu_c)
cat(sprintf("\n=== SE vs true ESD (after outlier removal) ===\n"))
cat(sprintf("  ESD_clean = %.4f\n", esd_clean))
cat(sprintf("  Fraction with se < ESD_clean: %.1f%%\n", 100*mean(se_c < esd_clean)))
cat(sprintf("  Fraction with se < 0.5*ESD_clean: %.1f%%\n", 100*mean(se_c < 0.5*esd_clean)))

cat("\nDone!\n")
