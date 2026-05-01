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
get_scenario_summary <- function(scenario_id,
                                  settings = NULL,
                                  version = NULL) {

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
    version = version
  )

  return(summary_tbls)
}

#------------------------------------------------------------------------------#
# Scenario 9-1: Correct Model (No-Y Models)
#------------------------------------------------------------------------------#

# cat("\n")
# cat(strrep("=", 70), "\n")
# cat("Scenario 9-1: Correct Model (No-Y Models)\n")
# cat(strrep("=", 70), "\n")

# summary_tbls <- get_scenario_summary("9-1", version = "test13")

# # Setting 11
# cat("\nSetting 11:\n")
# print(summary_tbls[[1]])

# # Setting 12
# cat("\nSetting 12:\n")
# print(summary_tbls[[2]])

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

# #------------------------------------------------------------------------------#
# # Scenario 7-1: Correct Model (Z1 Interactions)
# #------------------------------------------------------------------------------#

# cat("\n")
# cat(strrep("=", 70), "\n")
# cat("Scenario 7-1: Correct Model (Z1 Interactions)\n")
# cat(strrep("=", 70), "\n")

# summary_tbls <- get_scenario_summary("7-1", version = "test13")

# # Setting 11
# cat("\nSetting 11:\n")
# print(summary_tbls[[1]])

# # Setting 12
# cat("\nSetting 12:\n")
# print(summary_tbls[[2]])

# #------------------------------------------------------------------------------#
# # Scenario 7-2: Misspecified Model (Z1 Interactions)
# #------------------------------------------------------------------------------#

# cat("\n")
# cat(strrep("=", 70), "\n")
# cat("Scenario 7-2: Misspecified Model (Z1 Interactions)\n")
# cat(strrep("=", 70), "\n")

# summary_tbls <- get_scenario_summary("7-2", version = "test13")

# # Setting 11
# cat("\nSetting 11:\n")
# print(summary_tbls[[1]])

# # Setting 12
# cat("\nSetting 12:\n")
# print(summary_tbls[[2]])

# #------------------------------------------------------------------------------#
# # Scenario 7-3: Misspecified Model (Z1 Interactions)
# #------------------------------------------------------------------------------#

# cat("\n")
# cat(strrep("=", 70), "\n")
# cat("Scenario 7-3: Misspecified Model (Z1 Interactions)\n")
# cat(strrep("=", 70), "\n")

# summary_tbls <- get_scenario_summary("7-3", version = "test13")

# # Setting 11
# cat("\nSetting 11:\n")
# print(summary_tbls[[1]])

# # Setting 12
# cat("\nSetting 12:\n")
# print(summary_tbls[[2]])

# #------------------------------------------------------------------------------#
# # Scenario 7-4: Misspecified Model (Z1 Interactions)
# #------------------------------------------------------------------------------#

# cat("\n")
# cat(strrep("=", 70), "\n")
# cat("Scenario 7-4: Misspecified Model (Z1 Interactions)\n")
# cat(strrep("=", 70), "\n")

# summary_tbls <- get_scenario_summary("7-4", version = "test13")

# # Setting 11
# cat("\nSetting 11:\n")
# print(summary_tbls[[1]])

# # Setting 12
# cat("\nSetting 12:\n")
# print(summary_tbls[[2]])

# #------------------------------------------------------------------------------#
# # Scenario 8-1: Correct Model (Z2 Interactions)
# #------------------------------------------------------------------------------#

# cat("\n")
# cat(strrep("=", 70), "\n")
# cat("Scenario 8-1: Correct Model (Z2 Interactions)\n")
# cat(strrep("=", 70), "\n")

# summary_tbls <- get_scenario_summary("8-1", version = "test13")

# # Setting 11
# cat("\nSetting 11:\n")
# print(summary_tbls[[1]])

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
# Setting 13 Results (Scenarios 7-1, 7-2, 7-3, 7-4, 8-1, 9-1)
#------------------------------------------------------------------------------#

# cat("\n")
# cat(strrep("=", 70), "\n")
# cat("Setting 13 Results\n")
# cat(strrep("=", 70), "\n")

# for (scenario_id in c("7-1", "7-2", "7-3", "7-4", "8-1", "9-1")) {
#   cat("\n")
#   cat(strrep("-", 70), "\n")
#   cat("Scenario", scenario_id, "- Setting 13\n")
#   cat(strrep("-", 70), "\n")

#   tryCatch({
#     summary_tbls <- get_scenario_summary(scenario_id,
#                                           settings = "setting13",
#                                           version = "test14")
#     print(summary_tbls[[1]])
#   }, error = function(e) {
#     cat("Error:", e$message, "\n")
#   })
# }

#------------------------------------------------------------------------------#
# Setting 14 Results (Scenarios 7-1, 7-2, 7-3, 7-4, 8-1, 9-1)
#------------------------------------------------------------------------------#

# cat("\n")
# cat(strrep("=", 70), "\n")
# cat("Setting 14 Results\n")
# cat(strrep("=", 70), "\n")

# for (scenario_id in c("7-1", "7-2", "7-3", "7-4", "8-1", "9-1")) {
#   cat("\n")
#   cat(strrep("-", 70), "\n")
#   cat("Scenario", scenario_id, "- Setting 14\n")
#   cat(strrep("-", 70), "\n")

#   tryCatch({
#     summary_tbls <- get_scenario_summary(scenario_id,
#                                           settings = "setting14",
#                                           version = "test14")
#     print(summary_tbls[[1]])
#   }, error = function(e) {
#     cat("Error:", e$message, "\n")
#   })
# }

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

for (scenario in c("7-1", "8-1", "8-2", "9-1", "8-3", "8-4")) {
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat(paste("Scenario", scenario, "\n"))
  cat(strrep("=", 70), "\n")

  cho_settings <- c("setting13", "setting14")

  version <-  ifelse(scenario %in% c("8-2", "8-3", "8-4"), "test22", "test22")

  for (setting in cho_settings) {
    cat("\n")
    cat(strrep("-", 70), "\n")
    cat(setting, "\n")
    cat(strrep("-", 70), "\n")

    summary_tbls <- summarize_all_settings_with_all_missing_rates(
      settings = setting,
      missing_rates = c("miss50", "miss30"),
      scenario = scenario,
      J = 3,
      n.vector = c(500, 2000),
      all_data_file.list = correct_model_all_data_file.list,
      alpha_true.list = correct_model_alpha.true.list,
      version = version
    )

    print(summary_tbls[[1]])
  }
}

