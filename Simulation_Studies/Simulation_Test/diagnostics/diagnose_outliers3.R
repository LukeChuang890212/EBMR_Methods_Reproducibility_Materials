setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

f <- "Simulation_Results/EBMR_IPW_setting3-miss50-scenario9-3_1_n2000_replicate1000_test57.RDS"
res <- readRDS(f)

cat("File:", basename(f), "\n")
cat("Dim:", dim(res), "\n\n")

# Scenario 9-3: misspecified, ps_spec_id = "9-alt1"
# PS spec 9-alt1: formulas = (full, u1_z2, u2_z2), h_alpha = (full_sq1, full, full)
# full: r ~ y + u1 + u2 (4 params)
# u1_z2: r ~ y + u1 + z2 (4 params)
# u2_z2: r ~ y + u2 + z2 (4 params)
# This is model 1 only, so only the "full" formula is used
# full_sq1 h_alpha = (u1, u2, z1, z2, u2^2) => 5 h_alpha variables
# esteq_dim = max(5, 4) = 5 (but model 1 has 4 params, 5 h_alpha => overidentified)

mu <- res[1, ]   # mu_ipw
se <- res[3, ]   # se_ipw

valid <- !is.na(mu) & !is.na(se)
cat(sprintf("Valid: %d/%d\n", sum(valid), length(mu)))

mu_v <- mu[valid]; se_v <- se[valid]
q1 <- quantile(mu_v, 0.25); q3 <- quantile(mu_v, 0.75); iqr <- q3 - q1
is_out_mu <- mu_v < (q1 - 3*iqr) | mu_v > (q3 + 3*iqr)

q1s <- quantile(se_v, 0.25); q3s <- quantile(se_v, 0.75); iqrs <- q3s - q1s
is_out_se <- se_v < (q1s - 3*iqrs) | se_v > (q3s + 3*iqrs)

cat(sprintf("\nmu_ipw outliers (3IQR): %d (%.1f%%)\n", sum(is_out_mu), 100*mean(is_out_mu)))
cat(sprintf("se_ipw outliers (3IQR): %d (%.1f%%)\n", sum(is_out_se), 100*mean(is_out_se)))

# Look at the outlier reps' alpha values
alpha_rows <- grep("alpha\\.hat", rownames(res))
cat(sprintf("\nAlpha rows: %s\n", paste(rownames(res)[alpha_rows], collapse=", ")))

alpha <- res[alpha_rows, valid]
cat(sprintf("Alpha dim: %d x %d\n", nrow(alpha), ncol(alpha)))

# Compare alpha for SE outliers vs non-outliers
cat("\nAlpha stats for SE outliers vs non-outliers:\n")
for (k in 1:nrow(alpha)) {
  cat(sprintf("  %s:\n", rownames(res)[alpha_rows[k]]))
  cat(sprintf("    Outliers:     mean=%7.3f, sd=%6.3f, range=[%7.3f, %7.3f]\n",
              mean(alpha[k, is_out_se]), sd(alpha[k, is_out_se]),
              min(alpha[k, is_out_se]), max(alpha[k, is_out_se])))
  cat(sprintf("    Non-outliers: mean=%7.3f, sd=%6.3f, range=[%7.3f, %7.3f]\n",
              mean(alpha[k, !is_out_se]), sd(alpha[k, !is_out_se]),
              min(alpha[k, !is_out_se]), max(alpha[k, !is_out_se])))
}

# Now reproduce a few outlier cases to check convergence
cat("\n\n========== Reproducing outlier reps to check convergence ==========\n")

source("Data_Generation.r")
source("config/scenarios.R")
library(EPS)
library(numDeriv)

# Load the data
# Scenario 9-3 is misspecified, setting3, miss50
# Need to find the data file
source("Basic_setup.r")

data_file <- misspecified_model_all_data_file.list$setting3$miss50[[1]]
cat("Data file:", data_file, "\n")
all_data <- readRDS(data_file)

n_val <- 2000
ps_spec <- get_ps_spec("9-alt1")

# Subset to model 1 only
subset_ps_spec <- list(
  formula.list = ps_spec$formula.list[1],
  h_alpha.list = ps_spec$h_alpha.list[1],
  inv_link = ps_spec$inv_link,
  outcome = ps_spec$outcome
)

W <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

# Find outlier rep indices
outlier_idx <- which(valid)[is_out_se]
cat(sprintf("\nSE outlier rep indices (first 10): %s\n", paste(head(outlier_idx, 10), collapse=", ")))

# Run a few outlier and non-outlier reps
check_reps <- c(head(outlier_idx, 5), head(which(valid)[!is_out_se], 3))
cat(sprintf("Checking reps: %s\n\n", paste(check_reps, collapse=", ")))

for (rep_i in check_reps) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]

  is_outlier_rep <- rep_i %in% outlier_idx
  label <- if (is_outlier_rep) "OUTLIER" else "NORMAL"

  tryCatch({
    ebmr <- EPS$new("y", subset_ps_spec, dat, W)
    fit <- ebmr$ps_fit.list[[1]]
    gmm_fit <- fit$gmm_fit

    cat(sprintf("Rep %d [%s]: alpha=[%s], iter=%d, conv=%s, grad=%.2e, obj=%.6f, se_ipw_stored=%.4f\n",
                rep_i, label,
                paste(round(fit$coefficients, 3), collapse=", "),
                gmm_fit$opt$iterations,
                gmm_fit$opt$converged,
                gmm_fit$opt$final_grad_norm,
                gmm_fit$opt$objective,
                se[rep_i]))
  }, error = function(e) {
    cat(sprintf("Rep %d [%s]: ERROR - %s\n", rep_i, label, conditionMessage(e)))
  })
}

cat("\nDone!\n")
