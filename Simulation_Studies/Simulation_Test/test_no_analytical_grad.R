setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Data_Generation.r")
source("config/scenarios.R")
library(EBMRalgorithmFast4)
source("Basic_setup.r")

data_file <- misspecified_model_all_data_file.list$setting4$miss50[[1]]
all_data <- readRDS(data_file)
n_val <- 2000

ps_spec <- get_ps_spec("9-alt1")
subset_ps_spec <- list(
  formula.list = ps_spec$formula.list[c(1, 3)],
  h_alpha.list = ps_spec$h_alpha.list[c(1, 3)],
  inv_link = ps_spec$inv_link,
  outcome = ps_spec$outcome
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_nu_func <- function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2, u1_u2 = dat$u1*dat$u2)

# Load stored results for comparison
f <- "Simulation_Results/EBMR_IPW_setting4-miss50-scenario9-2_13_n2000_replicate1000_test57.RDS"
res <- readRDS(f)
mu_stored <- res[1, ]; w1_stored <- res["w.hat1", ]; w2_stored <- res["w.hat2", ]

# Identify outlier reps
valid <- !is.na(mu_stored)
mu_v <- mu_stored[valid]
q1 <- quantile(mu_v, 0.25); q3 <- quantile(mu_v, 0.75); iqr <- q3 - q1
is_out <- mu_v < (q1 - 3*iqr) | mu_v > (q3 + 3*iqr)
outlier_idx <- which(valid)[is_out]
normal_idx <- which(valid)[!is_out]

cat(sprintf("Found %d outlier reps\n", length(outlier_idx)))
cat("Testing outlier reps:", paste(head(outlier_idx, 10), collapse=", "), "\n\n")

test_reps <- c(head(outlier_idx, 10), head(normal_idx, 3))

for (rep_i in test_reps) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]
  is_outlier <- rep_i %in% outlier_idx
  label <- if (is_outlier) "OUTLIER" else "NORMAL"

  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", subset_ps_spec, dat, W_func)
    result <- ebmr$EBMR_IPW(h_nu = h_nu_func, se.fit = FALSE)
    cat(sprintf("Rep %3d [%s]: w=[%.4f,%.4f], mu=%.4f (stored: w=[%.4f,%.4f], mu=%.4f)\n",
                rep_i, label,
                result$w.hat[1], result$w.hat[2], result$mu_ipw,
                w1_stored[rep_i], w2_stored[rep_i], mu_stored[rep_i]))
  }, error = function(e) {
    cat(sprintf("Rep %3d [%s]: ERROR: %s\n", rep_i, label, conditionMessage(e)))
  })
}

cat("\nDone!\n")
