setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Data_Generation.r")
source("config/scenarios.R")
library(EPS)
source("Basic_setup.r")

f <- "Simulation_Results/EBMR_IPW_setting3-miss50-scenario9-3_2_n2000_replicate1000_test59.RDS"
res <- readRDS(f)
mu <- res[1, ]; se <- res[3, ]

# Reproduce a few extreme reps to check PS values
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

# Check extreme and normal reps
extreme_reps <- c(929, 438, 430, 221, 781, 587)  # worst outliers
normal_reps <- c(1, 2, 3, 4, 5)

cat("=== Checking extreme reps ===\n")
for (rep_i in c(extreme_reps, normal_reps)) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]

  tryCatch({
    ebmr <- EPS$new("y", subset_ps_spec, dat, W_func)
    ps_fitted <- ebmr$ps_fit.list[[1]]$fitted.values

    r_vec <- as.numeric(dat$r)
    y_vec <- dat$y
    ipw_weights <- r_vec / ps_fitted

    cat(sprintf("Rep %3d: mu=%.4f, se=%.4f | PS range=[%.4f,%.4f] | max(r/ps)=%.1f | mean(y)=%.3f\n",
                rep_i, mu[rep_i], se[rep_i],
                min(ps_fitted), max(ps_fitted),
                max(ipw_weights),
                mean(y_vec)))
  }, error = function(e) {
    cat(sprintf("Rep %3d: ERROR: %s\n", rep_i, conditionMessage(e)))
  })
}

# Overall: what fraction of reps have extreme IPW weights?
cat("\n=== IPW weight distribution across reps ===\n")
max_weights <- numeric(20)
for (rep_i in 1:20) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]
  ebmr <- EPS$new("y", subset_ps_spec, dat, W_func)
  ps_fitted <- ebmr$ps_fit.list[[1]]$fitted.values
  r_vec <- as.numeric(dat$r)
  max_weights[rep_i] <- max(r_vec / ps_fitted)
}
cat(sprintf("Max IPW weight across 20 reps: mean=%.1f, max=%.1f, min=%.1f\n",
            mean(max_weights), max(max_weights), min(max_weights)))

cat("\nDone!\n")
