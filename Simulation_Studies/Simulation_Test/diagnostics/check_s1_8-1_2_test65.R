## Check the s1 8-1 _2 (M2 only) file: ESE/ESD, max|alpha| distribution.
## Verify whether the M2 cnr/500 override was applied at generation time.
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({ source("Basic_setup.r"); source("Simulation.r") })

f <- "Simulation_Results/EBMR_IPW_setting1-miss50-scenario8-1_2_n2000_replicate1000_test65.RDS"
x <- readRDS(f)
cat(sprintf("Rows (%d): %s\n\n", nrow(x), paste(rownames(x), collapse=", ")))

mu_true <- get_mu_true("setting1")
cat(sprintf("mu_true(setting1) = %.4f\n\n", mu_true))

# RAW
mu_v <- x["mu_ipw", ]; se_v <- x["se_ipw", ]
valid <- !is.na(mu_v) & !is.na(se_v) & is.finite(mu_v) & is.finite(se_v)
cat(sprintf("== RAW (n_valid=%d / %d) ==\n", sum(valid), ncol(x)))
mu_v <- mu_v[valid]; se_v <- se_v[valid]
cat(sprintf("Bias=%+.4f  ESD=%.4f  ESE(mean)=%.4f  ESE(median)=%.4f  ratio=%.3f\n",
            mean(mu_v)-mu_true, sd(mu_v), mean(se_v), median(se_v), mean(se_v)/sd(mu_v)))

# CLEANED via clean_sim_result (mult=2, max_pct=0.01)
cl <- clean_sim_result(x, multiplier = 2, max_pct = 0.01, verbose = FALSE)
xc <- cl$result
mu_v <- xc["mu_ipw", ]; se_v <- xc["se_ipw", ]
cat(sprintf("\n== CLEANED (n=%d, NA=%d, Outliers=%d) ==\n",
            cl$n_successful, cl$n_na, cl$n_outliers))
cat(sprintf("Bias=%+.4f  ESD=%.4f  ESE(mean)=%.4f  ESE(median)=%.4f  ratio=%.3f\n",
            mean(mu_v)-mu_true, sd(mu_v), mean(se_v), median(se_v), mean(se_v)/sd(mu_v)))

# M2 alpha distribution -- 4 coefficients for M2 = u1_z1 + full
alpha_rows <- grep("^alpha\\.hat", rownames(x))
cat(sprintf("\nalpha rows: %s\n", paste(rownames(x)[alpha_rows], collapse=", ")))
if (length(alpha_rows) >= 4) {
  a <- x[alpha_rows[1:4], , drop = FALSE]
  max_a <- apply(a, 2, function(v) max(abs(v), na.rm = TRUE))
  cat(sprintf("max|alpha| across reps: median=%.2f  q90=%.2f  q99=%.2f  max=%.2f\n",
              median(max_a), quantile(max_a, 0.9), quantile(max_a, 0.99), max(max_a)))
  cat(sprintf("                       (cnr/500 expected: max <= ~5; default L-BFGS-B: max ~30-40)\n"))
}
