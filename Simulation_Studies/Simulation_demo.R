#------------------------------------------------------------------------------#
# Simulation Demo Script
#
# This script demonstrates how to run simulations using the EBMR framework.
# Modify the parameters below to customize your simulation run.
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Clear mu.true cache to ensure fresh computation
#------------------------------------------------------------------------------#
if (exists(".mu_true_cache", envir = .GlobalEnv)) {
  rm(".mu_true_cache", envir = .GlobalEnv)
  cat("mu.true cache cleared.\n")
}

#------------------------------------------------------------------------------#
# Install required packages if not already installed
#------------------------------------------------------------------------------#

# CRAN packages
required_packages <- c("dplyr", "knitr", "kableExtra", "parallel", "foreach",
                       "doSNOW", "doParallel", "stringr", "Matrix", "numDeriv",
                       "devtools")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# GitHub packages
if (!requireNamespace("EBMRalgorithmFast2", quietly = TRUE)) {
  devtools::install_github(
    "LukeChuang890212/EBMR_Methods_Reproducibility_Materials",
    subdir = "EBMRalgorithmFast2"
  )
}

#------------------------------------------------------------------------------#
# Load the simulation framework components
#------------------------------------------------------------------------------#
source("Basic_setup.r")
source("Data_Generation.r")
source("Simulation.r")
source("config/scenarios.R")

#------------------------------------------------------------------------------#
# Configuration
#------------------------------------------------------------------------------#

# Override default settings
CONFIG$replicate_num <- 1000    # Number of replicates (default: 1000)
CONFIG$n_default <- 300        # Sample size (default: 2000)

#------------------------------------------------------------------------------#
# View Available Options
#------------------------------------------------------------------------------#

# List all available scenarios
show_scenarios()

# Show current configuration
show_config()

#------------------------------------------------------------------------------#
# Step 1: Generate data (required before running simulations)
#------------------------------------------------------------------------------#

# For setting3 (continuous Y) - correct models, need A1 and A2 data at n=2000
generate_data("setting1.A1", n = 2000, replicate_num = 1000)
generate_data("setting1.A2", n = 2000, replicate_num = 1000)

# For setting4 (binary Y) - correct models, need A1 and A2 at n=2000
generate_data("setting2.A1", n = 2000, replicate_num = 1000)
generate_data("setting2.A2", n = 2000, replicate_num = 1000)

# For setting3 - misspecified models need B1 and B2 at n=2000 and n=500
generate_data("setting1.B1", n = 2000, replicate_num = 1000)
generate_data("setting1.B1", n = 500, replicate_num = 1000)
generate_data("setting1.B2", n = 2000, replicate_num = 1000)
generate_data("setting1.B2", n = 500, replicate_num = 1000)

# For setting4 - misspecified models need B1 and B2 at n=2000 and n=500
generate_data("setting2.B1", n = 2000, replicate_num = 1000)
generate_data("setting2.B1", n = 500, replicate_num = 1000)
generate_data("setting2.B2", n = 2000, replicate_num = 1000)
generate_data("setting2.B2", n = 500, replicate_num = 1000)

# # For setting11 - correct models (7-1, 8-1, 9-1), need A1 and A2 data at n=2000
# generate_data("setting13.A1", n = 2000, replicate_num = 1000)
# generate_data("setting13.A2", n = 2000, replicate_num = 1000)

# # For setting12 - correct models need A1 and A2 at n=2000
# generate_data("setting15.A1", n = 2000, replicate_num = 1000)
# generate_data("setting15.A2", n = 2000, replicate_num = 1000)

# # For setting11 - misspecified models need B1 and B2 at n=2000 and n=300
# generate_data("setting13.B1", n = 2000, replicate_num = 1000)
# generate_data("setting13.B1", n = 500, replicate_num = 1000)
# generate_data("setting13.B2", n = 2000, replicate_num = 1000)
# generate_data("setting13.B2", n = 500, replicate_num = 1000)

# # For setting12 - misspecified models need B1 and B2 at n=2000 and n=300
# generate_data("setting15.B1", n = 2000, replicate_num = 1000)
# generate_data("setting15.B1", n = 500, replicate_num = 1000)
# generate_data("setting15.B2", n = 2000, replicate_num = 1000)
# generate_data("setting15.B2", n = 500, replicate_num = 1000)

# # For setting16 - misspecified models (continuous Y, tilt = y+u1+u2)
# generate_data("setting16.B1", n = 2000, replicate_num = 1000)
# generate_data("setting16.B1", n = 500, replicate_num = 1000)
# generate_data("setting16.B2", n = 2000, replicate_num = 1000)
# generate_data("setting16.B2", n = 500, replicate_num = 1000)

# # For setting17 - misspecified models (binary Y, tilt = y+u1+u2)
# generate_data("setting17.B1", n = 2000, replicate_num = 1000)
# generate_data("setting17.B1", n = 500, replicate_num = 1000)
# generate_data("setting17.B2", n = 2000, replicate_num = 1000)
# generate_data("setting17.B2", n = 500, replicate_num = 1000)

# # For setting18 - misspecified models (continuous Y, tilt = u1+u2+z1+z2)
# generate_data("setting18.B1", n = 2000, replicate_num = 1000)
# generate_data("setting18.B1", n = 500, replicate_num = 1000)
# generate_data("setting18.B2", n = 2000, replicate_num = 1000)
# generate_data("setting18.B2", n = 500, replicate_num = 1000)

# # For setting19 - misspecified models (binary Y, tilt = u1+u2+z1+z2)
# generate_data("setting19.B1", n = 2000, replicate_num = 1000)
# generate_data("setting19.B1", n = 500, replicate_num = 1000)
# generate_data("setting19.B2", n = 2000, replicate_num = 1000)
# generate_data("setting19.B2", n = 500, replicate_num = 1000)

# For Cho settings - need A1 and A2 data at n=2000
generate_data("Cho_M1_gamma005.A1", n = 2000, replicate_num = 1000)
generate_data("Cho_M1_gamma005.A2", n = 2000, replicate_num = 1000)
generate_data("Cho_M1_gamma000.A1", n = 2000, replicate_num = 1000)
generate_data("Cho_M1_gamma000.A2", n = 2000, replicate_num = 1000)
generate_data("Cho_M2_gamma005.A1", n = 2000, replicate_num = 1000)
generate_data("Cho_M2_gamma005.A2", n = 2000, replicate_num = 1000)
generate_data("Cho_M2_gamma000.A1", n = 2000, replicate_num = 1000)
generate_data("Cho_M2_gamma000.A2", n = 2000, replicate_num = 1000)

#------------------------------------------------------------------------------#
# Step 2: Run scenarios 7-1, 7-2, 8-1, 8-2, 9-1, 9-2
# Summarize each scenario immediately after completion
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Single model simulation: 8-1, setting12, n=2000 and n=300, miss50, model 1 only
#------------------------------------------------------------------------------#

# # Get scenario and PS specification
# scenario <- get_scenario("8-1")
# ps_spec <- get_ps_spec(scenario$ps_spec_id)

# # Use only model 1
# model_set <- c(1)
# subset_ps_spec <- list(
#   formula.list = ps_spec$formula.list[model_set],
#   h_alpha.list = ps_spec$h_alpha.list[model_set],
#   inv_link = ps_spec$inv_link
# )

# # Get alpha.true for setting12, miss50
# alpha.true <- correct_model_alpha.true.list$setting12$miss50[[1]]
# cat("alpha.true:", alpha.true, "\n")

# # Define true PS model
# ps_model.true <- function(dat, alpha.true) {
#   1 / (1 + exp(cbind(rep(1, nrow(dat)), dat$y, dat$u1, dat$u2) %*% alpha.true))
# }

# # Run simulations for both sample sizes
# for (n in c(2000, 300)) {
#   cat("\n--- Running simulation for n =", n, "---\n")

#   # Load data for setting12, miss50 (A1 data for correct model)
#   all_data <- readRDS(paste0("Simulation_Data/setting12.A1_n2000_replicate1000.RDS"))

#   # Output file
#   save_file <- paste0("Simulation_Results/EBMR_IPW_setting12-miss50-scenario8-1_1_n", n, "_replicate1000_test6.RDS")

#   # Run simulation
#   simulate(
#     all_data = all_data,
#     ps_model.true = ps_model.true,
#     alpha.true = alpha.true,
#     ps_specifications = subset_ps_spec,
#     n = n,
#     replicate_num = 1000,
#     save_file = save_file
#   )

#   # Print summary
#   if (file.exists(save_file)) {
#     sim_result <- readRDS(save_file)
#     cleaned <- clean_sim_result(sim_result, multiplier = 30, verbose = FALSE)
#     cat("Summary: Total=", cleaned$n_total, ", Successful=", cleaned$n_successful,
#         ", NA=", cleaned$n_na, ", Outliers=", cleaned$n_outliers, "\n", sep = "")
#   }
# }

#------------------------------------------------------------------------------#
# Scenario 8-2: setting16 and setting17 (tilt = y+u1+u2)
#------------------------------------------------------------------------------#
# test46/48: include u2^2, z2^2 in h_alpha and h_nu
# test47/49: include only 1st moments in h_alpha and h_nu
# Each entry: list(scenario, settings)
# runs sequentially as listed
run_list <- list(
  list(scenario = "9-3", settings = c("setting2", "setting1")),
  list(scenario = "9-2", settings = c("setting1", "setting2")),
  # list(scenario = "9-1", settings = c("setting14")),
  # list(scenario = "9-1", settings = c("setting1", "setting2")),
  # list(scenario = "8-1", settings = c("setting1", "setting2")),
  # list(scenario = "7-1", settings = c("setting1", "setting2")),
  list(scenario = "cho1", settings = c("Cho_M1_gamma000")),
  list(scenario = "cho2", settings = c("Cho_M2_gamma000")),
  list(scenario = "cho1", settings = c("Cho_M1_gamma005")),
  list(scenario = "cho2", settings = c("Cho_M2_gamma005"))
)

# run_list <- list(
#   list(scenario = "9-1", settings = c("setting14")),
#   list(scenario = "9-2", settings = c("setting14")),
#   list(scenario = "9-1", settings = c("setting14")),
#   list(scenario = "7-1", settings = c("setting14")),
#   list(scenario = "8-1", settings = c("setting14")),
#   list(scenario = "9-3", settings = c("setting14"))
#   # list(scenario = "cho1", settings = c("Cho_M1_gamma000")),
#   # list(scenario = "cho2", settings = c("Cho_M2_gamma000")),
#   # list(scenario = "cho1", settings = c("Cho_M1_gamma005")),
#   # list(scenario = "cho2", settings = c("Cho_M2_gamma005"))
# )

all_latex_tables <- list()

for (entry in run_list) {
  scenario_id <- entry$scenario
  for (setting in entry$settings) {
    # Run the scenario
    run_scenario(scenario_id,
                 setting = setting,
                 version = "test65",
                 type = "HT")

    # Summarize immediately after completion
    cat("\n")
    cat(strrep("=", 70), "\n")
    cat("Summarizing scenario:", scenario_id, "for", setting, "\n")
    cat(strrep("=", 70), "\n")

    scenario <- SCENARIOS[[scenario_id]]
    params <- get_model_params(scenario$model_type)

    summary_table <- summarize_all_settings_with_all_missing_rates(
      settings = setting,
      missing_rates = c("miss50", "miss30"),
      scenario = scenario_id,
      J = 3,
      n.vector = unlist(params$n_vector),
      all_data_file.list = params$data_files,
      alpha_true.list = params$alpha_true,
      version = "test65",
      type = "HT"
    )

    all_latex_tables <- c(all_latex_tables, summary_table)
  }
}

#------------------------------------------------------------------------------#
# Print all LaTeX tables at the end
#------------------------------------------------------------------------------#
cat("\n\n")
cat(strrep("=", 70), "\n")
cat("                    ALL LATEX TABLES                                \n")
cat(strrep("=", 70), "\n\n")

for (tbl in all_latex_tables) {
  cat(tbl, "\n\n")
}

