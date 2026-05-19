setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Data_Generation.r")
source("config/scenarios.R")
library(EBMRalgorithmFast4)
source("Basic_setup.r")

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

# Check convergence info for outlier vs normal reps
extreme_reps <- c(929, 438, 430, 221, 180, 781, 754, 615, 711, 587)
normal_reps <- c(1, 2, 3, 4, 5, 10, 20, 50)

cat("=== Alpha-step convergence ===\n")
cat(sprintf("%-6s %-8s %-6s %-10s %-12s %-10s %-40s\n",
            "Rep", "Type", "Iter", "GradNorm", "Objective", "Converged", "Alpha"))

for (rep_i in c(extreme_reps, normal_reps)) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]
  label <- if (rep_i %in% extreme_reps) "OUT" else "NORM"

  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", subset_ps_spec, dat, W_func)
    fit <- ebmr$ps_fit.list[[1]]
    gmm <- fit$gmm_fit

    cat(sprintf("%-6d %-8s %-6d %-10.2e %-12.6e %-10s [%s]\n",
                rep_i, label,
                gmm$opt$iterations,
                gmm$opt$final_grad_norm,
                gmm$opt$objective,
                gmm$opt$converged,
                paste(round(fit$coefficients, 4), collapse=", ")))

    # Check PS distribution
    ps <- fit$fitted.values
    r_vec <- as.numeric(dat$r)
    ipw_w <- r_vec / ps
    cat(sprintf("       PS: min=%.6f, max=%.4f | max(r/ps)=%.1f | n(ps<0.01)=%d\n",
                min(ps), max(ps), max(ipw_w), sum(ps < 0.01)))
  }, error = function(e) {
    cat(sprintf("%-6d %-8s ERROR: %s\n", rep_i, label, conditionMessage(e)))
  })
}

cat("\nDone!\n")
