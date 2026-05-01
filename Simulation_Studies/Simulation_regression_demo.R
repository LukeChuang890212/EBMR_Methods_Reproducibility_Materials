#------------------------------------------------------------------------------#
# Simulation Demo Script for EBMR_IPW_regression
#
# This script demonstrates how to run simulations for regression coefficient
# estimation using the EBMR_IPW_regression function.
#
# Settings:
#   - Setting 13: Continuous y (Gaussian family)
#     y ~ N(0.2 + 0.6*z1 + 0.6*z2 + 0.3*u1 + 0.3*u2, 1)
#     True theta: (0.2, 0.6, 0.6, 0.3, 0.3)
#
#   - Setting 14: Binary y (Binomial family)
#     y ~ Bernoulli(expit(0.2 + 0.8*z1 + 0.8*z2 + 0.4*u1 + 0.4*u2))
#     True theta: (0.2, 0.8, 0.8, 0.4, 0.4)
#
# Both use correctly specified outcome model: y ~ z1 + z2 + u1 + u2
#
# Model Types:
#   - correct: Uses A1/A2 data (correctly specified PS model)
#   - misspecified: Uses B1/B2 data (misspecified PS model)
#------------------------------------------------------------------------------#

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
source("Simulation_regression.r")

#------------------------------------------------------------------------------#
# Configuration
#------------------------------------------------------------------------------#

# Override default settings
CONFIG_REGRESSION$replicate_num <- 1000    # Number of replicates (default: 1000)
CONFIG_REGRESSION$n_default <- 2000        # Sample size (default: 2000)

#------------------------------------------------------------------------------#
# View Available Options
#------------------------------------------------------------------------------#

# List all available regression scenarios
show_scenarios_regression()

# Show current configuration
show_config_regression()

#------------------------------------------------------------------------------#
# Step 1: Generate data (required before running simulations)
#------------------------------------------------------------------------------#

# For setting 13 (continuous y) - correct models need A1 and A2 data at n=2000
generate_data("setting13.A1", n = 2000, replicate_num = 1000)
generate_data("setting13.A2", n = 2000, replicate_num = 1000)

# For setting 14 (binary y) - correct models need A1 and A2 data at n=2000
generate_data("setting14.A1", n = 2000, replicate_num = 1000)
generate_data("setting14.A2", n = 2000, replicate_num = 1000)

# For setting 13 - misspecified models need B1 and B2 at n=2000 and n=500
generate_data("setting13.B1", n = 2000, replicate_num = 1000)
generate_data("setting13.B1", n = 500, replicate_num = 1000)
generate_data("setting13.B2", n = 2000, replicate_num = 1000)
generate_data("setting13.B2", n = 500, replicate_num = 1000)

# For setting 14 - misspecified models need B1 and B2 at n=2000 and n=500
generate_data("setting14.B1", n = 2000, replicate_num = 1000)
generate_data("setting14.B1", n = 500, replicate_num = 1000)
generate_data("setting14.B2", n = 2000, replicate_num = 1000)
generate_data("setting14.B2", n = 500, replicate_num = 1000)

#------------------------------------------------------------------------------#
# Step 2: Run scenarios
# Summarize each scenario immediately after completion
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Run scenarios 7-1, 8-1, 8-2, 9-1 for both settings
#------------------------------------------------------------------------------#

# Scenarios to run:
# 7-1: Correct model with z1 interactions (full, u1+z1, u2+z1)
# 8-1: Correct model with z1/z2 interactions (full, u1+z1, u2+z2)
# 8-2: Misspecified model with z1/z2 interactions (full, u1+z1, u2+z2)
# 9-1: Correct model with z2 interactions (full, u1+z2, u2+z2)

scenarios_to_run <- c("7-1", "8-1", "9-1")
settings_to_run <- c("setting13", "setting14")

for (scenario_id in scenarios_to_run) {
  for (setting in settings_to_run) {
    # Run the scenario
    run_scenario_regression(
      scenario_id,
      setting = setting,
      version = "v1"
    )

    # Summarize immediately after completion using wide table format
    generate_wide_table_regression(
      setting = setting,
      scenario_id = scenario_id,
      replicate_num = CONFIG_REGRESSION$replicate_num,
      version = "v1"
    )
  }
}

#------------------------------------------------------------------------------#
# Example 2: Run misspecified model scenario (reg-2) for setting13 only
#------------------------------------------------------------------------------#

# Uncomment to run:
# run_scenario_regression(
#   "reg-2",
#   setting = "setting13",
#   version = "v1"
# )
#
# summarize_all_regression(
#   settings = "setting13",
#   missing_rates = c("miss50", "miss30"),
#   scenario_id = "reg-2",
#   n.vector = c(500, 2000),
#   replicate_num = CONFIG_REGRESSION$replicate_num,
#   version = "v1"
# )

#------------------------------------------------------------------------------#
# Example 3: Run single model combination only (for debugging)
#------------------------------------------------------------------------------#

# # Load data
# all_data <- readRDS("Simulation_Data/setting14.A1_n2000_replicate1000.RDS")
#
# # Get PS specification
# ps_spec <- PS_SPEC_REGRESSION
#
# # Use only model 1
# model_set <- c(1)
# subset_ps_spec <- list(
#   formula.list = ps_spec$formula.list[model_set],
#   h_alpha.list = ps_spec$h_alpha.list[model_set],
#   inv_link = ps_spec$inv_link
# )
#
# # Get alpha.true for setting14, miss50
# alpha.true <- correct_model_alpha.true.list.regression$setting14$miss50[[1]]
# cat("alpha.true:", alpha.true, "\n")
#
# # Define true PS model
# ps_model.true <- function(dat, alpha.true) {
#   1 / (1 + exp(cbind(rep(1, nrow(dat)), dat$y, dat$u1, dat$u2) %*% alpha.true))
# }
#
# # Run simulation for n=2000
# save_file <- "Simulation_Results/EBMR_IPW_regression_setting14-miss50-debug_1_n2000_replicate100_test.RDS"
#
# simulate_regression(
#   all_data = all_data,
#   ps_model.true = ps_model.true,
#   alpha.true = alpha.true,
#   ps_specifications = subset_ps_spec,
#   n = 2000,
#   replicate_num = 100,
#   setting = "setting14",
#   save_file = save_file
# )
#
# # Print summary
# if (file.exists(save_file)) {
#   sim_result <- readRDS(save_file)
#   print_regression_summary(sim_result, "setting14", n = 2000, version = "test")
# }

#------------------------------------------------------------------------------#
# Example 4: Quick verification with a few replicates
#------------------------------------------------------------------------------#

# # Quick test with 50 replicates
# CONFIG_REGRESSION$replicate_num <- 50
#
# run_scenario_regression(
#   "reg-1",
#   setting = "setting13",
#   n = 2000,  # Only run n=2000
#   version = "quick_test"
# )
#
# # Check results
# sim_result <- readRDS("Simulation_Results/EBMR_IPW_regression_setting13-miss50-scenarioreg-1_123_n2000_replicate50_quick_test.RDS")
# print_regression_summary(sim_result, "setting13", n = 2000, version = "quick_test")

#------------------------------------------------------------------------------#
# Final Summary Table
#------------------------------------------------------------------------------#

cat("\n\n")
cat("================================================================================\n")
cat("                     REGRESSION SIMULATION COMPLETE\n")
cat("================================================================================\n\n")

cat("To view results, use:\n")
cat("  - summarize_all_regression() for all results\n")
cat("  - print_regression_summary() for detailed single file summary\n")
cat("  - readRDS() to load raw results\n\n")

cat("Output files are saved in: Simulation_Results/\n")
cat("File naming: EBMR_IPW_regression_<setting>-<miss>-scenario<id>_<models>_n<n>_replicate<num>_<version>.RDS\n\n")

cat("================================================================================\n")
