setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
source("Simulation.r")
library(EBMRalgorithmFast4)

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
results <- matrix(NA, 4, n_reps)
rownames(results) <- c("mu_ipw", "mu_ipw.true", "se_ipw", "se_ipw.true")
n_truncated <- rep(NA, n_reps)
alpha_norm <- rep(NA, n_reps)

cat("Running 1000 reps with ps truncation INSIDE GMM (ps_min=0.05)...\n")
n_err <- 0
for (rep_i in 1:n_reps) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]
  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", subset_ps_spec, dat, W_func)
    res <- ebmr$EBMR_IPW(
      h_nu = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2),
      se.fit = TRUE, type = "HT"
    )
    results[1, rep_i] <- res$mu_ipw
    results[3, rep_i] <- res$se_ipw

    # Count truncated obs and alpha norm
    ps_vec <- res$ps.matrix[, 1]
    n_truncated[rep_i] <- sum(ps_vec <= 0.05 + 1e-10)
    fit <- ebmr$ps_fit.list[[1]]
    alpha_norm[rep_i] <- sqrt(sum(fit$coefficients^2))
  }, error = function(e) {
    n_err <<- n_err + 1
    if (n_err <= 5) cat(sprintf("  Error rep %d: %s\n", rep_i, e$message))
  })
  if (rep_i %% 200 == 0) cat(sprintf("  %d / %d done (errors: %d)\n", rep_i, n_reps, n_err))
}

cat(sprintf("\nErrors: %d / %d\n", n_err, n_reps))

results[2, ] <- ifelse(is.na(results[1, ]), NA, 0)
results[4, ] <- ifelse(is.na(results[3, ]), NA, 0)

# Truncation summary
valid <- !is.na(n_truncated)
cat(sprintf("\n=== TRUNCATION SUMMARY (ps at floor 0.05) ===\n"))
cat(sprintf("Reps with any truncated obs: %d / %d\n", sum(n_truncated[valid] > 0), sum(valid)))
cat("Quantiles of # truncated obs per rep:\n")
print(quantile(n_truncated[valid], probs = c(0, 0.25, 0.50, 0.75, 0.90, 0.95, 1)))
cat(sprintf("Mean # truncated: %.1f\n", mean(n_truncated[valid])))

# Alpha norm summary
cat(sprintf("\n=== ALPHA NORM ===\n"))
cat(sprintf("Mean ||alpha||: %.4f, Max ||alpha||: %.4f\n",
            mean(alpha_norm[valid]), max(alpha_norm[valid])))

# Before outlier removal
mu_raw <- results[1, valid]
se_raw <- results[3, valid]

cat(sprintf("\n=== BEFORE OUTLIER REMOVAL (n=%d) ===\n", length(mu_raw)))
cat(sprintf("Bias = %.4f, ESD = %.4f, ESE = %.4f, ESE/ESD = %.2f\n",
            mean(mu_raw) - mu.true, sd(mu_raw), mean(se_raw), mean(se_raw)/sd(mu_raw)))
ci_l <- mu_raw - 1.96 * se_raw; ci_u <- mu_raw + 1.96 * se_raw
cat(sprintf("CP = %.4f\n", mean(mu.true >= ci_l & mu.true <= ci_u)))
cat(sprintf("Skewness = %.4f, Excess kurtosis = %.4f\n",
            mean((mu_raw - mean(mu_raw))^3) / sd(mu_raw)^3,
            mean((mu_raw - mean(mu_raw))^4) / sd(mu_raw)^4 - 3))

# After outlier removal
for (mult in c(1.5, 3)) {
  cleaned <- clean_sim_result(results, multiplier = mult, verbose = TRUE)
  sim <- cleaned$result
  mu_v <- sim[1, ]; se_v <- sim[3, ]

  bias <- mean(mu_v) - mu.true
  esd <- sd(mu_v)
  ese <- mean(se_v)
  ci_lower <- mu_v - 1.96 * se_v
  ci_upper <- mu_v + 1.96 * se_v
  cp <- mean(mu.true >= ci_lower & mu.true <= ci_upper)

  cat(sprintf("\n=== AFTER OUTLIER REMOVAL (mult=%.1f, n=%d) ===\n", mult, length(mu_v)))
  cat(sprintf("Bias = %.4f, ESD = %.4f, ESE = %.4f, ESE/ESD = %.2f, CP = %.4f\n",
              bias, esd, ese, ese/esd, cp))
  cat(sprintf("Skewness = %.4f, Excess kurtosis = %.4f\n",
              mean((mu_v - mean(mu_v))^3) / sd(mu_v)^3,
              mean((mu_v - mean(mu_v))^4) / sd(mu_v)^4 - 3))
}

cat("\nDone!\n")
