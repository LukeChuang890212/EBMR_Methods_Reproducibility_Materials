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
mu_vals <- rep(NA, n_reps)
min_pi <- rep(NA, n_reps)
max_rpi <- rep(NA, n_reps)
alpha_norm <- rep(NA, n_reps)
se_vals <- rep(NA, n_reps)
alpha_list <- vector("list", n_reps)

cat("Running 1000 reps to collect diagnostics...\n")
n_err <- 0
for (rep_i in 1:n_reps) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]
  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", subset_ps_spec, dat, W_func)
    res <- ebmr$EBMR_IPW(
      h_nu = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2),
      se.fit = TRUE, type = "HT"
    )
    mu_vals[rep_i] <- res$mu_ipw
    se_vals[rep_i] <- res$se_ipw

    # Extract PS info from fitted PS
    ps_vec <- res$ps.matrix[, 1]
    r_vec <- as.numeric(dat$r)
    min_pi[rep_i] <- min(ps_vec)
    max_rpi[rep_i] <- max(r_vec / ps_vec)

    fit <- ebmr$ps_fit.list[[1]]
    alpha <- fit$coefficients
    alpha_norm[rep_i] <- sqrt(sum(alpha^2))
    alpha_list[[rep_i]] <- alpha
  }, error = function(e) {
    n_err <<- n_err + 1
    if (n_err <= 5) cat(sprintf("  Error in rep %d: %s\n", rep_i, e$message))
  })
  if (rep_i %% 200 == 0) cat(sprintf("  %d / %d done (errors: %d)\n", rep_i, n_reps, n_err))
}

# Summary
valid <- !is.na(mu_vals)
cat(sprintf("\nValid reps: %d / %d\n", sum(valid), n_reps))
cat(sprintf("mu.true = %.4f\n", mu.true))
cat(sprintf("mean(mu) = %.4f, sd(mu) = %.4f\n", mean(mu_vals[valid]), sd(mu_vals[valid])))

# Sort by |mu - mu.true|
deviation <- abs(mu_vals - mu.true)
top_extreme <- order(deviation, decreasing = TRUE, na.last = TRUE)[1:20]

cat("\n=== TOP 20 EXTREME mu_ipw REPS ===\n")
cat(sprintf("%-6s %-10s %-10s %-12s %-12s %-10s\n",
            "Rep", "mu_ipw", "se_ipw", "min(pi)", "max(r/pi)", "||alpha||"))
for (i in top_extreme) {
  if (is.na(mu_vals[i])) next
  cat(sprintf("%-6d %-10.4f %-10.4f %-12.6f %-12.1f %-10.4f\n",
              i, mu_vals[i], se_vals[i], min_pi[i], max_rpi[i], alpha_norm[i]))
}

# Compare with "normal" reps (close to mu.true)
normal_reps <- order(deviation, na.last = TRUE)[1:10]
cat("\n=== TOP 10 NORMAL mu_ipw REPS ===\n")
cat(sprintf("%-6s %-10s %-10s %-12s %-12s %-10s\n",
            "Rep", "mu_ipw", "se_ipw", "min(pi)", "max(r/pi)", "||alpha||"))
for (i in normal_reps) {
  if (is.na(mu_vals[i])) next
  cat(sprintf("%-6d %-10.4f %-10.4f %-12.6f %-12.1f %-10.4f\n",
              i, mu_vals[i], se_vals[i], min_pi[i], max_rpi[i], alpha_norm[i]))
}

# Correlation analysis
cat("\n=== CORRELATIONS ===\n")
cat(sprintf("cor(|mu-mu.true|, min_pi):   %.3f\n", cor(deviation[valid], min_pi[valid], use = "complete")))
cat(sprintf("cor(|mu-mu.true|, max_rpi):  %.3f\n", cor(deviation[valid], max_rpi[valid], use = "complete")))
cat(sprintf("cor(|mu-mu.true|, ||alpha||): %.3f\n", cor(deviation[valid], alpha_norm[valid], use = "complete")))

# Quantiles of min_pi
cat("\n=== QUANTILES of min(pi) ===\n")
print(quantile(min_pi[valid], probs = c(0, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 1)))

cat("\n=== QUANTILES of max(r/pi) ===\n")
print(quantile(max_rpi[valid], probs = c(0.50, 0.75, 0.90, 0.95, 0.99, 1)))

# How many reps have min(pi) < thresholds
for (th in c(0.001, 0.01, 0.05, 0.1)) {
  cat(sprintf("Reps with min(pi) < %.3f: %d\n", th, sum(min_pi[valid] < th)))
}

# Print alpha for top 5 extreme reps
cat("\n=== ALPHA COEFFICIENTS for top 5 extreme reps ===\n")
cnt <- 0
for (i in top_extreme) {
  if (is.na(mu_vals[i])) next
  cnt <- cnt + 1
  if (cnt > 5) break
  cat(sprintf("Rep %d (mu=%.4f): alpha = [%s]\n",
              i, mu_vals[i], paste(round(alpha_list[[i]], 4), collapse=", ")))
}

cat("\nDone!\n")
