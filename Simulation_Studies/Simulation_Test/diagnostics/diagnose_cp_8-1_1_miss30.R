## Diagnose: is the CP shortfall in 8-1_1 miss30 driven by alpha blowups?
## Bin reps by max|alpha| and report CP per bin.
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({ source("Basic_setup.r"); source("Data_Generation.r"); source("Simulation.r") })

f <- "Simulation_Results/EBMR_IPW_setting2-miss30-scenario8-1_1_n2000_replicate1000_test65.RDS"
x <- readRDS(f)
set.seed(12345); d <- setting2.A1(1e7); mu_true <- mean(d$y)
cat(sprintf("mu_true = %.4f\n\n", mu_true))

mu_v <- x["mu_ipw", ]; se_v <- x["se_ipw", ]
alpha_rows <- grep("^alpha\\.hat", rownames(x))
a <- x[alpha_rows[1:4], , drop = FALSE]
max_a <- apply(a, 2, function(v) max(abs(v), na.rm = TRUE))

# CP per max|alpha| bin
ci_lo <- mu_v - 1.96*se_v; ci_hi <- mu_v + 1.96*se_v
covered <- (ci_lo <= mu_true) & (mu_true <= ci_hi)

bins <- list(
  list(label="max|a| <= 2",    mask = max_a <= 2),
  list(label="2 < max|a| <= 5",mask = max_a > 2 & max_a <= 5),
  list(label="5 < max|a| <=10",mask = max_a > 5 & max_a <= 10),
  list(label="max|a| > 10",    mask = max_a > 10)
)
cat(sprintf("%-22s | %5s | %s\n", "bin", "n", "Bias        ESD     ESE     ratio  CP"))
cat(strrep("-", 78), "\n", sep="")
for (b in bins) {
  m <- b$mask & !is.na(mu_v) & !is.na(se_v)
  if (sum(m) == 0) next
  mu_b <- mu_v[m]; se_b <- se_v[m]; cov_b <- covered[m]
  cat(sprintf("%-22s | %5d | %+.4f  %.4f  %.4f  %.3f  %.3f\n",
              b$label, sum(m),
              mean(mu_b) - mu_true, sd(mu_b), mean(se_b),
              mean(se_b)/sd(mu_b), mean(cov_b)))
}

cat("\n========= what if we drop the blowups? =========\n")
for (thr in c(20, 10, 5, 3, 2)) {
  m <- max_a <= thr & !is.na(mu_v) & !is.na(se_v)
  mu_b <- mu_v[m]; se_b <- se_v[m]
  cov_b <- covered[m]
  cat(sprintf("Drop reps with max|a| > %3g:  n=%4d  Bias=%+.4f  ESD=%.4f  ESE=%.4f  ratio=%.3f  CP=%.3f\n",
              thr, sum(m), mean(mu_b)-mu_true, sd(mu_b), mean(se_b),
              mean(se_b)/sd(mu_b), mean(cov_b)))
}

cat("\n========= clean_sim_result (mult=2, max_pct=0.01) =========\n")
cl <- clean_sim_result(x, multiplier = 2, max_pct = 0.01, verbose = FALSE)
cat(sprintf("Removed %d outliers + %d NA. The IQR-outlier removal targets mu_hat;\n",
            cl$n_outliers, cl$n_na))
cat("alpha blowups may give extreme mu_hat (outlier-removed) OR moderate mu_hat (kept).\n")
