#------------------------------------------------------------------------------#
# Simulation Results Summary Script - Cho et al. Scenarios
#
# This script generates summary tables for Cho et al. simulation results.
# The simulation results contain both Cho2025 and EBMRalgorithmOld output.
# We need to slice the results to get EBMRalgorithmOld part, then save them
# to temporary files, and use the existing summary functions.
#
# Slicing formula: [(nrow/2 + 1) : (nrow - 1)]
# After slicing:
#   - Row 1: mu_ipw
#   - Row 2: mu_ipw.true
#   - Row 3: se_ipw
#   - Row 4: se_ipw.true
#------------------------------------------------------------------------------#

# Load required packages
library(dplyr)
library(knitr)
library(kableExtra)

# Load the simulation framework components
source("Basic_setup.r")
source("Data_Generation.r")
source("Simulation.r")
source("config/scenarios.R")

#------------------------------------------------------------------------------#
# Step 1: Slice and save EBMRalgorithmOld results to temporary files
#------------------------------------------------------------------------------#

slice_and_save_cho_results <- function(settings,
                                        scenario = "cho1",
                                        missing_rates = c("miss50", "miss30"),
                                        n.vector = c(300, 2000),
                                        replicate_num = 1000,
                                        version = "test8") {

  model_sets <- c("1", "2", "3", "12", "13", "23", "123")

  for (setting in settings) {
    for (miss in missing_rates) {
      for (n in n.vector) {
        for (model_set in model_sets) {
          # Original file
          file <- paste0("Simulation_Results/EBMR_IPW_", setting, "-", miss,
                         "-scenario", scenario, "_", model_set, "_n", n,
                         "_replicate", replicate_num, "_", version, ".RDS")

          # Sliced file (temporary)
          sliced_file <- paste0("Simulation_Results/EBMR_IPW_", setting, "-", miss,
                                "-scenario", scenario, "_", model_set, "_n", n,
                                "_replicate", replicate_num, "_", version, "_sliced.RDS")

          if (file.exists(file)) {
            sim_result <- readRDS(file)
            n_rows <- nrow(sim_result)
            start_row <- n_rows / 2 + 1
            end_row <- n_rows - 1
            sliced_result <- sim_result[start_row:end_row, , drop = FALSE]
            saveRDS(sliced_result, sliced_file)
            cat("Sliced:", basename(file), "->", basename(sliced_file), "\n")
          }
        }
      }
    }
  }
}

#------------------------------------------------------------------------------#
# Step 2: Summarize using existing functions (with sliced version suffix)
#------------------------------------------------------------------------------#

summarize_cho_setting <- function(setting,
                                   scenario = "cho1",
                                   missing_rates = c("miss50", "miss30"),
                                   n.vector = c(300, 2000),
                                   replicate_num = 1000,
                                   version = "test8_sliced",
                                   is.original = FALSE) {

  mu.true <- get_mu_true(setting)

  summary_tbl <- summarize_all_settings_with_all_missing_rates(
    settings = setting,
    missing_rates = missing_rates,
    scenario = scenario,
    J = 3,
    n.vector = n.vector,
    all_data_file.list = list(),  # Not used for Cho settings
    alpha_true.list = list(),     # Not used for Cho settings
    version = version,
    is.original = is.original
  )

  return(summary_tbl)
}

#------------------------------------------------------------------------------#
# Run the pipeline
#------------------------------------------------------------------------------#

# Settings to process
cho_settings <- c("Cho_RM2q", "Cho_RM2p", "Cho_RM3q", "Cho_RM3p")

# Step 1: Slice and save
cat("\n")
cat(strrep("=", 70), "\n")
cat("Step 1: Slicing Cho simulation results\n")
cat(strrep("=", 70), "\n\n")

slice_and_save_cho_results(
  settings = cho_settings,
  version = "test8"
)

# Step 2: Summarize each setting
cat("\n")
cat(strrep("=", 70), "\n")
cat("Step 2: Summarizing results\n")
cat(strrep("=", 70), "\n")

for (setting in cho_settings) {
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat(setting, "\n")
  cat(strrep("=", 70), "\n")

  summary_tbls <- summarize_cho_setting(
    setting = setting,
    version = "test8_sliced",
    is.original = FALSE
  )

  print(summary_tbls[[1]])
}

#------------------------------------------------------------------------------#
# Step 3: Clean up - delete sliced files
#------------------------------------------------------------------------------#

cat("\n")
cat(strrep("=", 70), "\n")
cat("Step 3: Cleaning up sliced files\n")
cat(strrep("=", 70), "\n\n")

sliced_files <- list.files("Simulation_Results", pattern = "Cho.*sliced\\.RDS$", full.names = TRUE)
if (length(sliced_files) > 0) {
  file.remove(sliced_files)
  cat("Deleted", length(sliced_files), "sliced files\n")
} else {
  cat("No sliced files to delete\n")
}
