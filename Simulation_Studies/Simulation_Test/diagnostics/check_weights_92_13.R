setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

f <- "Simulation_Results/EBMR_IPW_setting4-miss50-scenario9-2_13_n2000_replicate1000_test57.RDS"
res <- readRDS(f)

cat("Rownames:", paste(rownames(res), collapse=", "), "\n\n")

mu <- res[1, ]; se <- res[3, ]
nu1 <- res["nu.hat1", ]; nu2 <- res["nu.hat2", ]
w1 <- res["w.hat1", ]; w2 <- res["w.hat2", ]

valid <- !is.na(mu) & !is.na(se)
mu_v <- mu[valid]; se_v <- se[valid]
nu1_v <- nu1[valid]; nu2_v <- nu2[valid]
w1_v <- w1[valid]; w2_v <- w2[valid]

# Identify outliers (mu or se)
q1 <- quantile(mu_v, 0.25); q3 <- quantile(mu_v, 0.75); iqr <- q3 - q1
is_out_mu <- mu_v < (q1 - 3*iqr) | mu_v > (q3 + 3*iqr)

q1s <- quantile(se_v, 0.25); q3s <- quantile(se_v, 0.75); iqrs <- q3s - q1s
is_out_se <- se_v < (q1s - 3*iqrs) | se_v > (q3s + 3*iqrs)

is_out <- is_out_mu | is_out_se

cat("=== nu.hat (raw GMM estimates) ===\n")
cat(sprintf("Outliers (%d reps):\n", sum(is_out)))
cat(sprintf("  nu1: mean=%.4f, sd=%.4f, range=[%.4f, %.4f]\n",
            mean(nu1_v[is_out]), sd(nu1_v[is_out]), min(nu1_v[is_out]), max(nu1_v[is_out])))
cat(sprintf("  nu2: mean=%.4f, sd=%.4f, range=[%.4f, %.4f]\n",
            mean(nu2_v[is_out]), sd(nu2_v[is_out]), min(nu2_v[is_out]), max(nu2_v[is_out])))

cat(sprintf("\nNormal (%d reps):\n", sum(!is_out)))
cat(sprintf("  nu1: mean=%.4f, sd=%.4f, range=[%.4f, %.4f]\n",
            mean(nu1_v[!is_out]), sd(nu1_v[!is_out]), min(nu1_v[!is_out]), max(nu1_v[!is_out])))
cat(sprintf("  nu2: mean=%.4f, sd=%.4f, range=[%.4f, %.4f]\n",
            mean(nu2_v[!is_out]), sd(nu2_v[!is_out]), min(nu2_v[!is_out]), max(nu2_v[!is_out])))

cat("\n=== w.hat (normalized weights: w_j = nu_j^2 / sum(nu_k^2)) ===\n")
cat(sprintf("Outliers (%d reps):\n", sum(is_out)))
cat(sprintf("  w1: mean=%.4f, sd=%.4f, range=[%.4f, %.4f]\n",
            mean(w1_v[is_out]), sd(w1_v[is_out]), min(w1_v[is_out]), max(w1_v[is_out])))
cat(sprintf("  w2: mean=%.4f, sd=%.4f, range=[%.4f, %.4f]\n",
            mean(w2_v[is_out]), sd(w2_v[is_out]), min(w2_v[is_out]), max(w2_v[is_out])))

cat(sprintf("\nNormal (%d reps):\n", sum(!is_out)))
cat(sprintf("  w1: mean=%.4f, sd=%.4f, range=[%.4f, %.4f]\n",
            mean(w1_v[!is_out]), sd(w1_v[!is_out]), min(w1_v[!is_out]), max(w1_v[!is_out])))
cat(sprintf("  w2: mean=%.4f, sd=%.4f, range=[%.4f, %.4f]\n",
            mean(w2_v[!is_out]), sd(w2_v[!is_out]), min(w2_v[!is_out]), max(w2_v[!is_out])))

# Correlation between w1 and mu/se
cat("\n=== Correlation ===\n")
cat(sprintf("cor(w1, mu_ipw) = %.4f\n", cor(w1_v, mu_v)))
cat(sprintf("cor(w1, se_ipw) = %.4f\n", cor(w1_v, se_v)))

# Which model is "correct"? Check scenario config
cat("\n=== Scenario info ===\n")
cat("File: scenario9-2, model_set=13 (models 1 and 3)\n")
cat("scenario9-2 is misspecified. Which model is correct?\n")

# Show distribution of w1
cat("\n=== w1 distribution (weight on model 1) ===\n")
cat(sprintf("Quantiles: 1%%=%.4f, 5%%=%.4f, 25%%=%.4f, 50%%=%.4f, 75%%=%.4f, 95%%=%.4f, 99%%=%.4f\n",
            quantile(w1_v, 0.01), quantile(w1_v, 0.05), quantile(w1_v, 0.25),
            quantile(w1_v, 0.50), quantile(w1_v, 0.75), quantile(w1_v, 0.95),
            quantile(w1_v, 0.99)))
cat(sprintf("w1 > 0.99: %d (%.1f%%)\n", sum(w1_v > 0.99), 100*mean(w1_v > 0.99)))
cat(sprintf("w1 > 0.95: %d (%.1f%%)\n", sum(w1_v > 0.95), 100*mean(w1_v > 0.95)))
cat(sprintf("w1 < 0.50: %d (%.1f%%)\n", sum(w1_v < 0.50), 100*mean(w1_v < 0.50)))

cat("\nDone!\n")
