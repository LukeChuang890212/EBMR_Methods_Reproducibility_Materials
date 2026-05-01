setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
source("Simulation.r")
library(EBMRalgorithmFast4)

# First, inspect the old results
old_file <- "Simulation_Results/EBMR_IPW_setting3-miss50-scenario9-3_2_n2000_replicate1000_test59.RDS"
old_result <- readRDS(old_file)
cat("=== OLD RESULTS (test59) ===\n")
cat(sprintf("Dimensions: %d x %d\n", nrow(old_result), ncol(old_result)))
cat(sprintf("Row names: %s\n", paste(rownames(old_result), collapse=", ")))

# Extract old stats
mu_old <- old_result[1, ]
se_old <- old_result[3, ]
valid_old <- !is.na(mu_old)
mu_true <- get_mu_true("setting3")

cat(sprintf("mu.true = %.4f\n", mu_true))
cat(sprintf("Valid reps: %d / %d\n", sum(valid_old), length(mu_old)))

bias_old <- mean(mu_old[valid_old]) - mu_true
esd_old <- sd(mu_old[valid_old])
ese_old <- mean(se_old[valid_old])
ci_l <- mu_old[valid_old] - 1.96 * se_old[valid_old]
ci_u <- mu_old[valid_old] + 1.96 * se_old[valid_old]
cp_old <- mean(mu_true >= ci_l & mu_true <= ci_u)

cat(sprintf("Bias=%.4f, ESD=%.4f, ESE=%.4f, ESE/ESD=%.2f, CP=%.4f\n\n",
    bias_old, esd_old, ese_old, ese_old/esd_old, cp_old))

# Now run with current package
# Setting 3, scenario 9-3, model 2 only
data_file <- misspecified_model_all_data_file.list[["setting3"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
n_val <- 2000

ps_spec <- get_ps_spec("9-alt1")
# Model 2 only (index 2 in ps_spec)
single_ps_spec <- list(
  formula.list = ps_spec[["formula.list"]][2],
  h_alpha.list = ps_spec[["h_alpha.list"]][2],
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]]
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

n_reps <- 1000
results_new <- matrix(NA, 4, n_reps)
rownames(results_new) <- c("mu_ipw", "mu_ipw.true", "se_ipw", "se_ipw.true")

cat("=== RUNNING CURRENT PACKAGE (1000 reps) ===\n")
n_err <- 0
for (rep_i in 1:n_reps) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]
  tryCatch({
    ebmr <- EBMRAlgorithmFast4[["new"]]("y", single_ps_spec, dat, W_func)
    res <- ebmr[["EBMR_IPW"]](
      h_nu = function(dat) cbind(u1 = dat[["u1"]], u2 = dat[["u2"]], z1 = dat[["z1"]], z2 = dat[["z2"]]),
      se.fit = TRUE, type = "HT"
    )
    results_new[1, rep_i] <- res[["mu_ipw"]]
    results_new[3, rep_i] <- res[["se_ipw"]]
  }, error = function(e) {
    n_err <<- n_err + 1
    if (n_err <= 5) cat(sprintf("  Error rep %d: %s\n", rep_i, conditionMessage(e)))
  })
  if (rep_i %% 200 == 0) cat(sprintf("  %d / %d done (errors: %d)\n", rep_i, n_reps, n_err))
}

cat(sprintf("\nErrors: %d / %d\n", n_err, n_reps))

valid_new <- !is.na(results_new[1, ])
mu_new <- results_new[1, valid_new]
se_new <- results_new[3, valid_new]

bias_new <- mean(mu_new) - mu_true
esd_new <- sd(mu_new)
ese_new <- mean(se_new)
ci_l_new <- mu_new - 1.96 * se_new
ci_u_new <- mu_new + 1.96 * se_new
cp_new <- mean(mu_true >= ci_l_new & mu_true <= ci_u_new)

cat(sprintf("\n=== COMPARISON ===\n"))
cat(sprintf("%-12s %10s %10s %10s %10s %10s\n", "", "Bias", "ESD", "ESE", "ESE/ESD", "CP"))
cat(sprintf("%-12s %10.4f %10.4f %10.4f %10.2f %10.4f\n", "Old (test59)", bias_old, esd_old, ese_old, ese_old/esd_old, cp_old))
cat(sprintf("%-12s %10.4f %10.4f %10.4f %10.2f %10.4f\n", "New", bias_new, esd_new, ese_new, ese_new/esd_new, cp_new))
