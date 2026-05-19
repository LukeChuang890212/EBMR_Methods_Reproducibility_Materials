setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

f <- "Simulation_Results/EBMR_IPW_setting4-miss50-scenario9-2_13_n2000_replicate1000_test57.RDS"
res <- readRDS(f)
cat("Dim:", dim(res), "\n")
cat("Rownames:", paste(rownames(res), collapse=", "), "\n\n")

# Summary stats for first 4 rows
for (k in 1:min(4, nrow(res))) {
  v <- res[k, ]
  valid <- !is.na(v)
  q1 <- quantile(v[valid], 0.25); q3 <- quantile(v[valid], 0.75); iqr <- q3 - q1
  out3 <- sum(v[valid] < (q1 - 3*iqr) | v[valid] > (q3 + 3*iqr))
  out30 <- sum(v[valid] < (q1 - 30*iqr) | v[valid] > (q3 + 30*iqr))
  cat(sprintf("%-15s: mean=%10.4f, sd=%8.4f, min=%10.4f, max=%10.4f, out(3IQR)=%d, out(30IQR)=%d, NA=%d\n",
              rownames(res)[k], mean(v[valid]), sd(v[valid]), min(v[valid]), max(v[valid]), out3, out30, sum(!valid)))
}

# Identify outlier reps (based on mu_ipw)
mu <- res[1, ]; se <- res[3, ]
valid <- !is.na(mu) & !is.na(se)
mu_v <- mu[valid]; se_v <- se[valid]

q1 <- quantile(mu_v, 0.25); q3 <- quantile(mu_v, 0.75); iqr <- q3 - q1
is_out_mu <- mu_v < (q1 - 3*iqr) | mu_v > (q3 + 3*iqr)

q1s <- quantile(se_v, 0.25); q3s <- quantile(se_v, 0.75); iqrs <- q3s - q1s
is_out_se <- se_v < (q1s - 3*iqrs) | se_v > (q3s + 3*iqrs)

is_out_any <- is_out_mu | is_out_se

cat(sprintf("\nmu_ipw outliers (3IQR): %d (%.1f%%)\n", sum(is_out_mu), 100*mean(is_out_mu)))
cat(sprintf("se_ipw outliers (3IQR): %d (%.1f%%)\n", sum(is_out_se), 100*mean(is_out_se)))
cat(sprintf("Any outlier (3IQR): %d (%.1f%%)\n", sum(is_out_any), 100*mean(is_out_any)))

# Extreme values
cat("\nSmallest 10 mu_ipw:", paste(round(head(sort(mu_v), 10), 4), collapse=", "), "\n")
cat("Largest 10 mu_ipw: ", paste(round(tail(sort(mu_v), 10), 4), collapse=", "), "\n")
cat("Smallest 10 se_ipw:", paste(round(head(sort(se_v), 10), 4), collapse=", "), "\n")
cat("Largest 10 se_ipw: ", paste(round(tail(sort(se_v), 10), 4), collapse=", "), "\n")

# ======== Reproduce outlier reps to check convergence ========
cat("\n\n========== Reproducing outlier reps to check convergence ==========\n")

source("Data_Generation.r")
source("config/scenarios.R")
library(EPS)
library(numDeriv)
source("Basic_setup.r")

# scenario9-2 is misspecified, setting4, miss50
# model_set = "13" means models 1 and 3
data_file <- misspecified_model_all_data_file.list$setting4$miss50[[1]]
cat("Data file:", data_file, "\n")
all_data <- readRDS(data_file)

n_val <- 2000
ps_spec <- get_ps_spec("9-alt1")

# Subset to models 1 and 3
subset_ps_spec <- list(
  formula.list = ps_spec$formula.list[c(1, 3)],
  h_alpha.list = ps_spec$h_alpha.list[c(1, 3)],
  inv_link = ps_spec$inv_link,
  outcome = ps_spec$outcome
)

W <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

# Find outlier rep indices (from mu or se)
outlier_idx <- which(valid)[is_out_any]
normal_idx <- which(valid)[!is_out_any]
cat(sprintf("\nOutlier rep indices (first 10): %s\n", paste(head(outlier_idx, 10), collapse=", ")))

check_reps <- c(head(outlier_idx, 5), head(normal_idx, 3))
cat(sprintf("Checking reps: %s\n\n", paste(check_reps, collapse=", ")))

for (rep_i in check_reps) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]
  is_outlier_rep <- rep_i %in% outlier_idx
  label <- if (is_outlier_rep) "OUTLIER" else "NORMAL"

  tryCatch({
    ebmr <- EPS$new("y", subset_ps_spec, dat, W)

    # Collect convergence info from both PS models
    iters <- sapply(ebmr$ps_fit.list, function(f) f$gmm_fit$opt$iterations)
    grads <- sapply(ebmr$ps_fit.list, function(f) f$gmm_fit$opt$final_grad_norm)
    convs <- sapply(ebmr$ps_fit.list, function(f) f$gmm_fit$opt$converged)
    objs <- sapply(ebmr$ps_fit.list, function(f) f$gmm_fit$opt$objective)

    cat(sprintf("Rep %d [%s]: iter=[%s], conv=[%s], grad=[%s], obj=[%s], mu=%.4f, se=%.4f\n",
                rep_i, label,
                paste(iters, collapse=","),
                paste(convs, collapse=","),
                paste(sprintf("%.2e", grads), collapse=","),
                paste(sprintf("%.4f", objs), collapse=","),
                mu[rep_i], se[rep_i]))
  }, error = function(e) {
    cat(sprintf("Rep %d [%s]: ERROR - %s\n", rep_i, label, conditionMessage(e)))
  })
}

cat("\nDone!\n")
