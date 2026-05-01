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

# Run 200 reps comparing HT vs Hajek (offset by 500 to use different data draws)
offset <- 500
n_reps <- 200
results_ht <- matrix(NA, 2, n_reps)   # mu, se
results_hj <- matrix(NA, 2, n_reps)

for (rep_i in 1:n_reps) {
  idx <- rep_i + offset
  dat <- all_data[((idx - 1) * n_val + 1):(idx * n_val), ]

  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", subset_ps_spec, dat, W_func)

    res_ht <- ebmr$EBMR_IPW(
      h_nu = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2),
      se.fit = TRUE, type = "HT"
    )
    results_ht[1, rep_i] <- res_ht$mu_ipw
    results_ht[2, rep_i] <- res_ht$se_ipw

    res_hj <- ebmr$EBMR_IPW(
      h_nu = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2),
      se.fit = TRUE, type = "Hajek"
    )
    results_hj[1, rep_i] <- res_hj$mu_ipw
    results_hj[2, rep_i] <- res_hj$se_ipw
  }, error = function(e) {
    cat(sprintf("Rep %d: ERROR: %s\n", rep_i, conditionMessage(e)))
  })
}

valid <- !is.na(results_ht[1,]) & !is.na(results_hj[1,])

cat("=== HT estimator ===\n")
cat(sprintf("  mean(mu)=%.4f, ESD=%.4f, mean(ESE)=%.4f, ratio=%.2f\n",
            mean(results_ht[1, valid]), sd(results_ht[1, valid]),
            mean(results_ht[2, valid]), mean(results_ht[2, valid])/sd(results_ht[1, valid])))
cat(sprintf("  mu range=[%.4f, %.4f]\n", min(results_ht[1, valid]), max(results_ht[1, valid])))

cat("\n=== Hajek estimator ===\n")
cat(sprintf("  mean(mu)=%.4f, ESD=%.4f, mean(ESE)=%.4f, ratio=%.2f\n",
            mean(results_hj[1, valid]), sd(results_hj[1, valid]),
            mean(results_hj[2, valid]), mean(results_hj[2, valid])/sd(results_hj[1, valid])))
cat(sprintf("  mu range=[%.4f, %.4f]\n", min(results_hj[1, valid]), max(results_hj[1, valid])))

cat("\nDone!\n")
