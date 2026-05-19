setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
source("Simulation.r")
library(EPS)

data_file <- misspecified_model_all_data_file.list$setting3$miss50[[1]]
all_data <- readRDS(data_file)
n_val <- 2000

ps_spec <- get_ps_spec("9-alt1")
subset_ps_spec <- list(
  formula.list = ps_spec$formula.list[2],
  h_alpha.list = ps_spec$h_alpha.list[2],
  inv_link = ps_spec$inv_link,
  outcome = ps_spec$outcome
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
mu.true <- get_mu_true("setting3")

n_reps <- 1000
mu_vals <- rep(NA, n_reps)
se_vals <- rep(NA, n_reps)

cat("Running 1000 reps...\n")
for (rep_i in 1:n_reps) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]
  tryCatch({
    ebmr <- EPS$new("y", subset_ps_spec, dat, W_func)
    res <- ebmr$EBMR_IPW(
      h_nu = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2),
      se.fit = TRUE, type = "HT"
    )
    mu_vals[rep_i] <- res$mu_ipw
    se_vals[rep_i] <- res$se_ipw
  }, error = function(e) {})
  if (rep_i %% 200 == 0) cat(sprintf("  %d / %d done\n", rep_i, n_reps))
}

valid <- !is.na(mu_vals)
mu_raw <- mu_vals[valid]
se_raw <- se_vals[valid]

# After outlier removal (mult=1.5)
results <- matrix(NA, 4, n_reps)
rownames(results) <- c("mu_ipw", "mu_ipw.true", "se_ipw", "se_ipw.true")
results[1, ] <- mu_vals
results[2, ] <- ifelse(is.na(mu_vals), NA, 0)
results[3, ] <- se_vals
results[4, ] <- ifelse(is.na(se_vals), NA, 0)
cleaned <- clean_sim_result(results, multiplier = 1.5, verbose = TRUE)
sim <- cleaned$result
mu_clean <- sim[1, ]
se_clean <- sim[3, ]

cat(sprintf("\n=== RAW (n=%d) ===\n", length(mu_raw)))
cat(sprintf("mean = %.4f, sd = %.4f\n", mean(mu_raw), sd(mu_raw)))
cat(sprintf("skewness = %.4f\n", mean((mu_raw - mean(mu_raw))^3) / sd(mu_raw)^3))
cat(sprintf("kurtosis = %.4f (excess = %.4f)\n", mean((mu_raw - mean(mu_raw))^4) / sd(mu_raw)^4, mean((mu_raw - mean(mu_raw))^4) / sd(mu_raw)^4 - 3))
sw <- shapiro.test(mu_raw[1:min(5000, length(mu_raw))])
cat(sprintf("Shapiro-Wilk p-value = %.4e\n", sw$p.value))

cat(sprintf("\n=== AFTER OUTLIER REMOVAL (n=%d) ===\n", length(mu_clean)))
cat(sprintf("mean = %.4f, sd = %.4f\n", mean(mu_clean), sd(mu_clean)))
cat(sprintf("skewness = %.4f\n", mean((mu_clean - mean(mu_clean))^3) / sd(mu_clean)^3))
cat(sprintf("kurtosis = %.4f (excess = %.4f)\n", mean((mu_clean - mean(mu_clean))^4) / sd(mu_clean)^4, mean((mu_clean - mean(mu_clean))^4) / sd(mu_clean)^4 - 3))
sw2 <- shapiro.test(mu_clean[1:min(5000, length(mu_clean))])
cat(sprintf("Shapiro-Wilk p-value = %.4e\n", sw2$p.value))

# Quantile comparison with normal
cat("\n=== QUANTILE COMPARISON (after outlier removal) ===\n")
mu_std <- (mu_clean - mean(mu_clean)) / sd(mu_clean)
probs <- c(0.01, 0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975, 0.99)
q_emp <- quantile(mu_std, probs)
q_norm <- qnorm(probs)
cat(sprintf("%-8s %-10s %-10s\n", "prob", "empirical", "normal"))
for (i in seq_along(probs)) {
  cat(sprintf("%-8.3f %-10.4f %-10.4f\n", probs[i], q_emp[i], q_norm[i]))
}

# Tail analysis
cat("\n=== TAIL ANALYSIS (after outlier removal) ===\n")
mu_centered <- mu_clean - mean(mu_clean)
esd <- sd(mu_clean)
for (k in c(1, 1.5, 2, 2.5, 3)) {
  n_beyond <- sum(abs(mu_centered) > k * esd)
  pct <- n_beyond / length(mu_clean) * 100
  expected_pct <- 2 * pnorm(-k) * 100
  cat(sprintf("  |z| > %.1f: %d (%.1f%%), expected under normal: %.1f%%\n",
              k, n_beyond, pct, expected_pct))
}

# QQ plot
png("Simulation_Results/mu_ipw_qqplot.png", width = 1200, height = 600)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# Raw
qqnorm(mu_raw, main = sprintf("Q-Q plot: mu_ipw raw (n=%d)", length(mu_raw)),
       pch = 16, cex = 0.5)
qqline(mu_raw, col = "red", lwd = 2)

# Cleaned
qqnorm(mu_clean, main = sprintf("Q-Q plot: mu_ipw cleaned (n=%d)", length(mu_clean)),
       pch = 16, cex = 0.5)
qqline(mu_clean, col = "red", lwd = 2)

dev.off()

# Histogram with normal overlay
png("Simulation_Results/mu_ipw_hist_normal.png", width = 1200, height = 600)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

hist(mu_raw, breaks = 50, prob = TRUE,
     main = sprintf("mu_ipw raw (n=%d)", length(mu_raw)),
     xlab = "mu_ipw", col = "lightblue", border = "white")
curve(dnorm(x, mean(mu_raw), sd(mu_raw)), add = TRUE, col = "red", lwd = 2)
abline(v = mu.true, col = "darkgreen", lwd = 2, lty = 2)

hist(mu_clean, breaks = 50, prob = TRUE,
     main = sprintf("mu_ipw cleaned (n=%d)", length(mu_clean)),
     xlab = "mu_ipw", col = "lightgreen", border = "white")
curve(dnorm(x, mean(mu_clean), sd(mu_clean)), add = TRUE, col = "red", lwd = 2)
abline(v = mu.true, col = "darkgreen", lwd = 2, lty = 2)

dev.off()

cat("\nPlots saved to Simulation_Results/mu_ipw_qqplot.png and mu_ipw_hist_normal.png\n")
cat("Done!\n")
