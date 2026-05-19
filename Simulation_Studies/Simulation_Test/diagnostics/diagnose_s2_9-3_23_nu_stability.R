## Diagnose ESE/ESD discrepancy for EBMR_IPW_setting2-miss50-scenario9-3_23_n2000_replicate1000_test65.RDS
## Hypothesis: nu estimation is very unstable (rep-to-rep), driving extra variability into mu_ipw
## that the analytic SE (computed at the rep's nu hat) doesn't capture.
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({ source("Basic_setup.r"); source("Data_Generation.r"); source("Simulation.r") })

f <- "Simulation_Results/EBMR_IPW_setting2-miss50-scenario9-3_23_n2000_replicate1000_test65.RDS"
x <- readRDS(f)
mu_true <- get_mu_true("setting2")

cat(sprintf("=== %s ===\n", basename(f)))
cat("dim:", dim(x), "\n")
cat("rownames:", paste(rownames(x), collapse=", "), "\n\n")

# --- ESE/ESD with cleaning ---
cleaned <- clean_sim_result(x, multiplier = 2, max_pct = 0.01, verbose = FALSE)
sc <- cleaned$result
mu_v <- sc["mu_ipw", ]; se_v <- sc["se_ipw", ]
bias <- mean(mu_v) - mu_true; esd <- sd(mu_v); ese <- mean(se_v); ese_med <- median(se_v)
ratio <- ese / esd
ci_lo <- mu_v - 1.96*se_v; ci_hi <- mu_v + 1.96*se_v
cp <- mean(ci_lo <= mu_true & mu_true <= ci_hi)
cat(sprintf("[clean] n=%d Out=%d  Bias=%+.4f  ESD=%.4f  ESE(mean)=%.4f  ESE(median)=%.4f  ratio=%.3f  CP=%.3f\n",
            cleaned$n_successful, cleaned$n_outliers, bias, esd, ese, ese_med, ratio, cp))
cat(sprintf("[raw]   n=%d        Bias=%+.4f  ESD=%.4f  ESE(mean)=%.4f  ratio=%.3f\n",
            ncol(x), mean(x['mu_ipw',]) - mu_true, sd(x['mu_ipw',]), mean(x['se_ipw',]),
            mean(x['se_ipw',])/sd(x['mu_ipw',])))

# --- nu.hat distribution ---
cat("\n--- nu.hat distribution across 1000 reps ---\n")
nu1 <- x["nu.hat1", ]; nu2 <- x["nu.hat2", ]
cat(sprintf("nu.hat1  (M2):  mean=%.3f  sd=%.3f  min=%.3f  q25=%.3f  med=%.3f  q75=%.3f  max=%.3f\n",
            mean(nu1, na.rm=T), sd(nu1, na.rm=T),
            min(nu1, na.rm=T), quantile(nu1, 0.25, na.rm=T), median(nu1, na.rm=T),
            quantile(nu1, 0.75, na.rm=T), max(nu1, na.rm=T)))
cat(sprintf("nu.hat2  (M3):  mean=%.3f  sd=%.3f  min=%.3f  q25=%.3f  med=%.3f  q75=%.3f  max=%.3f\n",
            mean(nu2, na.rm=T), sd(nu2, na.rm=T),
            min(nu2, na.rm=T), quantile(nu2, 0.25, na.rm=T), median(nu2, na.rm=T),
            quantile(nu2, 0.75, na.rm=T), max(nu2, na.rm=T)))

# --- w.hat distribution (normalized weights) ---
cat("\n--- w.hat distribution ---\n")
w1 <- x["w.hat1", ]; w2 <- x["w.hat2", ]
cat(sprintf("w.hat1   (M2):  mean=%.3f  sd=%.3f  min=%.3f  q25=%.3f  med=%.3f  q75=%.3f  max=%.3f\n",
            mean(w1, na.rm=T), sd(w1, na.rm=T),
            min(w1, na.rm=T), quantile(w1, 0.25, na.rm=T), median(w1, na.rm=T),
            quantile(w1, 0.75, na.rm=T), max(w1, na.rm=T)))
cat(sprintf("w.hat2   (M3):  mean=%.3f  sd=%.3f  min=%.3f  q25=%.3f  med=%.3f  q75=%.3f  max=%.3f\n",
            mean(w2, na.rm=T), sd(w2, na.rm=T),
            min(w2, na.rm=T), quantile(w2, 0.25, na.rm=T), median(w2, na.rm=T),
            quantile(w2, 0.75, na.rm=T), max(w2, na.rm=T)))

# --- Tabulate w.hat1 in bins to see how often it "flips" ---
cat("\n--- w.hat1 (M2 weight) histogram, 10 bins ---\n")
h <- hist(w1, breaks = seq(0, 1, by = 0.1), plot = FALSE)
for (i in seq_along(h$counts))
  cat(sprintf("  [%.1f, %.1f):  %d reps\n", h$breaks[i], h$breaks[i+1], h$counts[i]))

# --- Is mu_ipw correlated with w.hat? ---
cor_mu_w <- cor(x['mu_ipw',], w1, use="complete.obs")
cor_se_w <- cor(x['se_ipw',], w1, use="complete.obs")
cat(sprintf("\nCor(mu_ipw, w_M2) = %.3f   (large => mu_ipw shifts when weight flips between M2/M3)\n", cor_mu_w))
cat(sprintf("Cor(se_ipw, w_M2) = %.3f\n", cor_se_w))

# --- Decompose: what would ESD be if we conditioned on w? ---
# Split reps by which model dominates (w_M2 > 0.5 vs <= 0.5)
mu_v_all <- x['mu_ipw', ]; w_all <- w1
m2_dom <- w_all > 0.5
cat(sprintf("\nDecomposition by dominant weight:\n"))
cat(sprintf("  M2 dominant (w>0.5):  n=%d   mean(mu)=%.4f   sd(mu)=%.4f\n",
            sum(m2_dom, na.rm=T), mean(mu_v_all[m2_dom], na.rm=T), sd(mu_v_all[m2_dom], na.rm=T)))
cat(sprintf("  M3 dominant (w<=0.5): n=%d   mean(mu)=%.4f   sd(mu)=%.4f\n",
            sum(!m2_dom, na.rm=T), mean(mu_v_all[!m2_dom], na.rm=T), sd(mu_v_all[!m2_dom], na.rm=T)))
cat(sprintf("  Overall sd(mu) = %.4f\n", sd(mu_v_all, na.rm=T)))
cat(sprintf("  -> If most variance is between-cluster (M2 vs M3), then nu flipping is the culprit.\n"))
