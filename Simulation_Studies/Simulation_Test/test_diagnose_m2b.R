setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})

# Check: what does the actual package produce for setting4 model2 and setting3 model2?
# Use the package's WangShaoKim2014 function directly.

library(EBMRalgorithmFast4)

n_val <- 2000
n_reps <- 200

run_package <- function(setting, miss, ps_spec_id, model_idx) {
  ps_spec <- get_ps_spec(ps_spec_id)
  data_file <- misspecified_model_all_data_file.list[[setting]][[miss]][[1]]
  all_data <- readRDS(data_file)
  mu_true <- get_mu_true(setting)

  formula_j <- ps_spec[["formula.list"]][[model_idx]]
  h_alpha_vars <- ps_spec[["h_alpha.list"]][[model_idx]]

  mu_vals <- rep(NA_real_, n_reps)
  se1_vals <- rep(NA_real_, n_reps)
  se2_vals <- rep(NA_real_, n_reps)

  for (rep_i in 1:n_reps) {
    dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
    tryCatch({
      result <- WangShaoKim2014(
        formula = formula_j,
        outcome = dat$y,
        h_alpha = h_alpha_vars,
        inv_link = "logistic",
        data = dat,
        se.fit = TRUE
      )
      mu_vals[rep_i] <- result$mu
      se1_vals[rep_i] <- result$se1
      se2_vals[rep_i] <- result$se2
    }, error = function(e) NULL)
  }

  valid <- !is.na(mu_vals) & !is.na(se1_vals)
  esd <- sd(mu_vals[valid])
  ese1 <- mean(se1_vals[valid])
  ese2 <- mean(se2_vals[valid], na.rm=TRUE)

  cat(sprintf("\n=== Package: %s, model %d ===\n", setting, model_idx))
  cat(sprintf("Valid: %d / %d\n", sum(valid), n_reps))
  cat(sprintf("Bias=%.4f, ESD=%.4f\n", mean(mu_vals[valid]) - mu_true, esd))
  cat(sprintf("ESE1=%.4f, ESE1/ESD=%.3f\n", ese1, ese1/esd))
  cat(sprintf("ESE2=%.4f, ESE2/ESD=%.3f\n", ese2, ese2/esd))
}

run_package("setting4", "miss50", "9-alt1", 2)
run_package("setting4", "miss50", "9-alt1", 3)
run_package("setting3", "miss50", "9-alt1", 2)
