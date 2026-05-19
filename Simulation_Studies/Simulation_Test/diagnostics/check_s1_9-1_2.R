## Check the s1 9-1 _2 (M2 = u1_z2 + full) file: ratio, max|alpha|, CP.
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({ source("Basic_setup.r"); source("Data_Generation.r"); source("Simulation.r") })

f <- "Simulation_Results/EBMR_IPW_setting1-miss50-scenario9-1_2_n2000_replicate1000_test65.RDS"
if (!file.exists(f)) { cat("MISSING:", f, "\n"); quit(status=1) }
x <- readRDS(f)
cat(sprintf("File: %s (%.1f KB)\n", basename(f), file.size(f)/1024))
cat(sprintf("Rows (%d): %s\n\n", nrow(x), paste(rownames(x), collapse=", ")))

# Compute mu_true for setting1 (A1 data, continuous y)
set.seed(12345); d <- setting1.A1(1e7); mu_true <- mean(d$y)
cat(sprintf("mu_true(setting1.A1, n=1e7) = %.4f\n\n", mu_true))

report <- function(label, x) {
  mu_v <- x[1, ]; se_v <- x[3, ]
  valid <- !is.na(mu_v) & !is.na(se_v) & is.finite(mu_v) & is.finite(se_v)
  mu_v <- mu_v[valid]; se_v <- se_v[valid]
  bias <- mean(mu_v) - mu_true; esd <- sd(mu_v)
  ese_mean <- mean(se_v); ese_med <- median(se_v); ratio <- ese_mean / esd
  ci_lo <- mu_v - 1.96*se_v; ci_hi <- mu_v + 1.96*se_v
  cp <- mean(ci_lo <= mu_true & mu_true <= ci_hi)
  cat(sprintf("== %s (n=%d) ==\n", label, length(mu_v)))
  cat(sprintf("Bias=%+.4f  ESD=%.4f  ESE(mean)=%.4f  ESE(med)=%.4f  ratio=%.3f  CP=%.3f\n",
              bias, esd, ese_mean, ese_med, ratio, cp))
}

report("RAW (NA/Inf only)", x)
cl <- clean_sim_result(x, multiplier = 2, max_pct = 0.01, verbose = FALSE)
cat(sprintf("\n(clean_sim_result: n_total=%d, n_na=%d, n_outliers=%d, n_used=%d)\n",
            cl$n_total, cl$n_na, cl$n_outliers, cl$n_successful))
report("CLEANED (mult=2)", cl$result)

# M2 alpha distribution -- 4 coefficients
alpha_rows <- grep("^alpha\\.hat", rownames(x))
if (length(alpha_rows) >= 4) {
  a <- x[alpha_rows[1:4], , drop = FALSE]
  max_a <- apply(a, 2, function(v) max(abs(v), na.rm = TRUE))
  cat(sprintf("\nM2 (u1_z2) max|alpha|: median=%.2f  q90=%.2f  q99=%.2f  max=%.2f\n",
              median(max_a, na.rm=T), quantile(max_a, 0.9, na.rm=T),
              quantile(max_a, 0.99, na.rm=T), max(max_a, na.rm=T)))
  cat(sprintf("Number of reps with max|alpha| > 10: %d / %d\n",
              sum(max_a > 10, na.rm=T), length(max_a)))
  for (k in 1:4) {
    cat(sprintf("  %-30s  mean=%+.3f  sd=%.3f  min=%+.3f  max=%+.3f\n",
                rownames(x)[alpha_rows[k]],
                mean(a[k,], na.rm=T), sd(a[k,], na.rm=T),
                min(a[k,], na.rm=T), max(a[k,], na.rm=T)))
  }
}

# Distribution check for se_ipw
cat("\nse_ipw distribution:\n")
print(quantile(x[3,], c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1.0), na.rm=TRUE))
