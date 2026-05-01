setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

f <- "Simulation_Results/EBMR_IPW_setting3-miss50-scenario9-3_2_n2000_replicate1000_test59.RDS"
res <- readRDS(f)
mu <- res[1,]; se <- res[3,]
valid <- !is.na(mu) & !is.na(se)
mu_v <- mu[valid]; se_v <- se[valid]

q1 <- quantile(mu_v, 0.25); q3 <- quantile(mu_v, 0.75); iqr <- q3 - q1
q1s <- quantile(se_v, 0.25); q3s <- quantile(se_v, 0.75); iqrs <- q3s - q1s
is_out <- (mu_v < q1 - 3*iqr | mu_v > q3 + 3*iqr) | (se_v < q1s - 3*iqrs | se_v > q3s + 3*iqrs)

# Alpha magnitudes
alpha <- rbind(res["alpha.hat1", valid], res["alpha.hat2", valid],
               res["alpha.hat3", valid], res["alpha.hat4", valid])
alpha_norm <- apply(alpha, 2, function(a) sqrt(sum(a^2)))

cat("=== Alpha norm (L2) ===\n")
cat(sprintf("  Outliers (n=%d):     mean=%.1f, median=%.1f, range=[%.1f, %.1f]\n",
    sum(is_out), mean(alpha_norm[is_out]), median(alpha_norm[is_out]),
    min(alpha_norm[is_out]), max(alpha_norm[is_out])))
cat(sprintf("  Non-outliers (n=%d): mean=%.1f, median=%.1f, range=[%.1f, %.1f]\n",
    sum(!is_out), mean(alpha_norm[!is_out]), median(alpha_norm[!is_out]),
    min(alpha_norm[!is_out]), max(alpha_norm[!is_out])))

cat(sprintf("\n  Outliers with ||alpha|| > 5: %d / %d\n", sum(alpha_norm[is_out] > 5), sum(is_out)))
cat(sprintf("  Outliers with ||alpha|| > 10: %d / %d\n", sum(alpha_norm[is_out] > 10), sum(is_out)))
cat(sprintf("  Non-outliers with ||alpha|| > 5: %d / %d\n", sum(alpha_norm[!is_out] > 5), sum(!is_out)))
cat(sprintf("  Non-outliers with ||alpha|| > 10: %d / %d\n", sum(alpha_norm[!is_out] > 10), sum(!is_out)))

# Show a few outliers with small alpha norm (if any)
small_alpha_out <- which(is_out & alpha_norm <= 5)
if (length(small_alpha_out) > 0) {
  cat(sprintf("\n=== Outliers with ||alpha|| <= 5 (%d total) ===\n", length(small_alpha_out)))
  for (i in small_alpha_out[1:min(10, length(small_alpha_out))]) {
    cat(sprintf("  mu=%.4f, se=%.4f, ||alpha||=%.2f\n", mu_v[i], se_v[i], alpha_norm[i]))
  }
}

cat("\nDone!\n")
