#------------------------------------------------------------------------------#
# Simulation Results Summary Script
#
# This script generates summary tables for simulation results.
# Modify the scenarios below to customize which results to display.
#------------------------------------------------------------------------------#

# Load required packages
library(dplyr)
library(knitr)
library(kableExtra)

# Load the simulation framework components (not Simulation_main.r which runs simulations)
source("Basic_setup.r")
source("Data_Generation.r")
source("Simulation.r")
source("config/scenarios.R")

#------------------------------------------------------------------------------#
# Helper Function
#------------------------------------------------------------------------------#

#' Generate summary table for a scenario
#'
#' @param scenario_id Scenario ID (e.g., "9-1", "7-2")
#' @param settings Settings to summarize (default: from scenario definition)
#' @param version Version string (default: from scenario definition)
#' @param is.original Whether to use original data without outlier removal
get_scenario_summary <- function(scenario_id,
                                  settings = NULL,
                                  version = NULL,
                                  is.original = FALSE) {

  # Get scenario configuration
  scenario <- SCENARIOS[[scenario_id]]
  if (is.null(scenario)) {
    stop(paste("Unknown scenario:", scenario_id))
  }

  # Use scenario defaults if not provided
  if (is.null(settings)) settings <- scenario$settings
  if (is.null(version)) version <- scenario$version

  # Get model parameters based on model type
  params <- get_model_params(scenario$model_type)

  # Generate summary table
  summary_tbls <- summarize_all_settings_with_all_missing_rates(
    settings = settings,
    missing_rates = c("miss50", "miss30"),
    scenario = scenario_id,
    J = 3,
    n.vector = unlist(params$n_vector),
    all_data_file.list = params$data_files,
    alpha_true.list = params$alpha_true,
    version = version,
    is.original = is.original
  )

  return(summary_tbls)
}

#------------------------------------------------------------------------------#
# Scenario 9-1: Correct Model (No-Y Models)
#------------------------------------------------------------------------------#

cat("\n")
cat(strrep("=", 70), "\n")
cat("Scenario 9-1: Correct Model (No-Y Models)\n")
cat(strrep("=", 70), "\n")

summary_tbls <- get_scenario_summary("9-1", version = "test5")

# Setting 11
cat("\nSetting 11:\n")
print(summary_tbls[[1]])

# Setting 12
cat("\nSetting 12:\n")
print(summary_tbls[[2]])

#------------------------------------------------------------------------------#
# Scenario 9-2: Misspecified Model (No-Y Models)
#------------------------------------------------------------------------------#

# cat("\n")
# cat(strrep("=", 70), "\n")
# cat("Scenario 9-2: Misspecified Model (No-Y Models)\n")
# cat(strrep("=", 70), "\n")
#
# summary_tbls <- get_scenario_summary("9-2")
#
# # Setting 11
# cat("\nSetting 11:\n")
# print(summary_tbls[[1]])
#
# # Setting 12
# cat("\nSetting 12:\n")
# print(summary_tbls[[2]])

#------------------------------------------------------------------------------#
# Scenario 7-1: Correct Model (Z1 Interactions)
#------------------------------------------------------------------------------#

cat("\n")
cat(strrep("=", 70), "\n")
cat("Scenario 7-1: Correct Model (Z1 Interactions)\n")
cat(strrep("=", 70), "\n")

summary_tbls <- get_scenario_summary("7-1", version = "test5")

# Setting 11
cat("\nSetting 11:\n")
print(summary_tbls[[1]])

# Setting 12
cat("\nSetting 12:\n")
print(summary_tbls[[2]])

#------------------------------------------------------------------------------#
# Scenario 7-2: Misspecified Model (Z1 Interactions)
#------------------------------------------------------------------------------#

# cat("\n")
# cat(strrep("=", 70), "\n")
# cat("Scenario 7-2: Misspecified Model (Z1 Interactions)\n")
# cat(strrep("=", 70), "\n")
#
# summary_tbls <- get_scenario_summary("7-2")
#
# # Setting 11
# cat("\nSetting 11:\n")
# print(summary_tbls[[1]])
#
# # Setting 12
# cat("\nSetting 12:\n")
# print(summary_tbls[[2]])

#------------------------------------------------------------------------------#
# Scenario 8-1: Correct Model (Z2 Interactions)
#------------------------------------------------------------------------------#

# cat("\n")
# cat(strrep("=", 70), "\n")
# cat("Scenario 8-1: Correct Model (Z2 Interactions)\n")
# cat(strrep("=", 70), "\n")
#
# summary_tbls <- get_scenario_summary("8-1")
#
# # Setting 11
# cat("\nSetting 11:\n")
# print(summary_tbls[[1]])
#
# # Setting 12
# cat("\nSetting 12:\n")
# print(summary_tbls[[2]])

#------------------------------------------------------------------------------#
# Scenario 8-2: Misspecified Model (Z2 Interactions)
#------------------------------------------------------------------------------#

# cat("\n")
# cat(strrep("=", 70), "\n")
# cat("Scenario 8-2: Misspecified Model (Z2 Interactions)\n")
# cat(strrep("=", 70), "\n")
#
# summary_tbls <- get_scenario_summary("8-2")
#
# # Setting 11
# cat("\nSetting 11:\n")
# print(summary_tbls[[1]])
#
# # Setting 12
# cat("\nSetting 12:\n")
# print(summary_tbls[[2]])

#------------------------------------------------------------------------------#
# Legacy Scenarios (1-x series)
#------------------------------------------------------------------------------#

# cat("\n")
# cat(strrep("=", 70), "\n")
# cat("Scenario 1-1\n")
# cat(strrep("=", 70), "\n")
#
# summary_tbls <- get_scenario_summary("1-1",
#                                       settings = c("setting11", "setting12"))
# print(summary_tbls[[1]])
# print(summary_tbls[[2]])

# cat("\n")
# cat(strrep("=", 70), "\n")
# cat("Scenario 1-2\n")
# cat(strrep("=", 70), "\n")
#
# summary_tbls <- get_scenario_summary("1-2",
#                                       settings = c("setting11", "setting12"))
# print(summary_tbls[[1]])
# print(summary_tbls[[2]])

# cat("\n")
# cat(strrep("=", 70), "\n")
# cat("Scenario 1-3\n")
# cat(strrep("=", 70), "\n")
#
# summary_tbls <- get_scenario_summary("1-3",
#                                       settings = c("setting11", "setting12"))
# print(summary_tbls[[1]])
# print(summary_tbls[[2]])

#------------------------------------------------------------------------------#
# Cho et al. Scenarios
#------------------------------------------------------------------------------#

# cat("\n")
# cat(strrep("=", 70), "\n")
# cat("Cho et al. Scenarios\n")
# cat(strrep("=", 70), "\n")
#
# summary_tbls <- summarize_all_settings_with_all_missing_rates(
#   settings = c("Cho_RM2q", "Cho_RM2p", "Cho_RM3q", "Cho_RM3p"),
#   missing_rates = c("miss50", "miss30"),
#   scenario = "cho1",
#   J = 3,
#   n.vector = c(1000, 4000),
#   all_data_file.list = misspecified_model_all_data_file.list,
#   alpha_true.list = misspecified_model_alpha.true.list,
#   version = "test7",
#   is.original = FALSE
# )
# print(summary_tbls[[1]])
# print(summary_tbls[[2]])
# print(summary_tbls[[3]])
# print(summary_tbls[[4]])
