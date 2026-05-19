## Check the s2 8-1 _1 (M1 only, the CORRECT model) file: Bias/ESD/ESE/CP, max|alpha|.
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({ source("Basic_setup.r"); source("Data_Generation.r"); source("Simulation.r") })

f <- "Simulation_Results/EBMR_IPW_setting2-miss30-scenario8-1_1_n2000_replicate1000_test65.RDS"
if (!file.exists(f)) { cat("MISSING:", f, "\n"); quit(status=1) }
x <- readRDS(f)
cat(sprintf("File: %s (%.1f KB)\n", basename(f), file.size(f)/1024))
cat(sprintf("Rows (%d): %s\n\n", nrow(x), paste(rownames(x), collapse=", ")))

# Compute true mu_true by generating fresh sample
set.seed(12345)
d <- setting2.A1(1e7)
mu_true <- mean(d$y)
cat(sprintf("mu_true(setting2.A1, n=1e7) = %.4f\n\n", mu_true))

report <- function(label, x) {
  mu_v <- x["mu_ipw", ]; se_v <- x["se_ipw", ]
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

report("RAW (NA/Inf filtered)", x)

cl <- clean_sim_result(x, multiplier = 2, max_pct = 0.01, verbose = FALSE)
cat(sprintf("\n(clean_sim_result: n_total=%d, n_na=%d, n_outliers=%d, n_used=%d)\n",
            cl$n_total, cl$n_na, cl$n_outliers, cl$n_successful))
report("CLEANED (mult=2)", cl$result)

# Same with true-PS row (rows 2/4) for reference
mu_t <- x["mu_ipw.true", ]; se_t <- x["se_ipw.true", ]
v <- !is.na(mu_t) & !is.na(se_t) & is.finite(mu_t) & is.finite(se_t)
mu_t <- mu_t[v]; se_t <- se_t[v]
ci_lo_t <- mu_t - 1.96*se_t; ci_hi_t <- mu_t + 1.96*se_t
cat(sprintf("\n[reference] true-PS IPW: Bias=%+.4f  ESD=%.4f  ESE=%.4f  ratio=%.3f  CP=%.3f\n",
            mean(mu_t)-mu_true, sd(mu_t), mean(se_t), mean(se_t)/sd(mu_t),
            mean(ci_lo_t <= mu_true & mu_true <= ci_hi_t)))

# M1 alpha distribution
alpha_rows <- grep("^alpha\\.hat", rownames(x))
if (length(alpha_rows) >= 4) {
  a <- x[alpha_rows[1:4], , drop = FALSE]
  max_a <- apply(a, 2, function(v) max(abs(v), na.rm = TRUE))
  cat(sprintf("\nM1 max|alpha|: median=%.2f  q90=%.2f  q99=%.2f  max=%.2f\n",
              median(max_a, na.rm=T), quantile(max_a, 0.9, na.rm=T),
              quantile(max_a, 0.99, na.rm=T), max(max_a, na.rm=T)))
  # alpha names
  cat("alpha coef summary (mean ± sd):\n")
  for (k in 1:4) {
    cat(sprintf("  %-30s  mean=%+.3f  sd=%.3f  min=%+.3f  max=%+.3f\n",
                rownames(x)[alpha_rows[k]],
                mean(a[k,], na.rm=T), sd(a[k,], na.rm=T),
                min(a[k,], na.rm=T), max(a[k,], na.rm=T)))
  }
}
