setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
source("Simulation.r")
library(EBMRalgorithmFast4)

data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
n_val <- 2000

ps_spec <- get_ps_spec("9-alt1")
subset_ps_spec <- list(
  formula.list = ps_spec[["formula.list"]][c(1, 3)],
  h_alpha.list = ps_spec[["h_alpha.list"]][c(1, 3)],
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]]
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

test_reps <- c(1, 2, 3, 8, 28, 83, 92, 133, 155)

cat("=== HYBRID CUE REFINEMENT TEST ===\n\n")
for (rep_i in test_reps) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]
  is_outlier <- rep_i %in% c(8, 28, 83, 92, 133, 155)

  tryCatch({
    ebmr <- EBMRAlgorithmFast4[["new"]]("y", subset_ps_spec, dat, W_func)
    res <- ebmr[["EBMR_IPW"]](
      h_nu = function(dat) cbind(u1 = dat[["u1"]], u2 = dat[["u2"]], z1 = dat[["z1"]], z2 = dat[["z2"]]),
      se.fit = TRUE, type = "HT"
    )

    cat(sprintf("--- Rep %d [%s] ---\n", rep_i, ifelse(is_outlier, "OUTLIER", "normal")))
    for (m in 1:2) {
      fit <- ebmr[["ps_fit.list"]][[m]]
      cat(sprintf("  Model %d: ||a||=%.2f conv=%s iter=%d grad=%.2e obj=%.6f cue_ref=%s\n",
          m, sqrt(sum(fit[["coefficients"]]^2)),
          fit[["opt"]][["converged"]], fit[["opt"]][["iterations"]],
          fit[["opt"]][["final_grad_norm"]], fit[["opt"]][["objective"]],
          fit[["opt"]][["cue_refined"]]))
    }

    ens_fit <- res[["ensemble_fit"]]
    cat(sprintf("  Ensemble: w=(%.4f,%.4f) conv=%s iter=%d grad=%.2e obj=%.6f cue_ref=%s\n",
        res[["w.hat"]][1], res[["w.hat"]][2],
        ens_fit[["gmm_fit"]][["opt"]][["converged"]], ens_fit[["gmm_fit"]][["opt"]][["iterations"]],
        ens_fit[["gmm_fit"]][["opt"]][["final_grad_norm"]], ens_fit[["gmm_fit"]][["opt"]][["objective"]],
        ens_fit[["gmm_fit"]][["opt"]][["cue_refined"]]))
    cat(sprintf("  mu_ipw = %.4f\n\n", res[["mu_ipw"]]))
  }, error = function(e) {
    cat(sprintf("Rep %d: ERROR - %s\n\n", rep_i, conditionMessage(e)))
  })
}
