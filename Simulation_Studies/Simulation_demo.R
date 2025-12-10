#------------------------------------------------------------------------------#
# Simulation Demo Script
#
# This script demonstrates how to run simulations using the EBMR framework.
# Modify the parameters below to customize your simulation run.
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
if (!requireNamespace("EBMRalgorithmOld", quietly = TRUE)) {
  devtools::install_github(
    "LukeChuang890212/EBMR_Methods_Reproducibility_Materials",
    subdir = "EBMRalgorithmOld"
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

# For setting11 - correct models (7-1, 8-1, 9-1), need A1 and A2 data at n=2000
# generate_data("setting11.A1", n = 2000, replicate_num = 1000)
# generate_data("setting11.A2", n = 2000, replicate_num = 1000)

# For setting12 - correct models need A1 and A2 at n=2000
# generate_data("setting12.A1", n = 2000, replicate_num = 1000)
# generate_data("setting12.A2", n = 2000, replicate_num = 1000)

# For setting11 - misspecified models need B1 and B2 at n=2000 and n=300
# generate_data("setting11.B1", n = 2000, replicate_num = 1000)
# generate_data("setting11.B1", n = 300, replicate_num = 1000)
# generate_data("setting11.B2", n = 2000, replicate_num = 1000)
# generate_data("setting11.B2", n = 300, replicate_num = 1000)

# For setting12 - misspecified models need B1 and B2 at n=2000 and n=300
# generate_data("setting12.B1", n = 2000, replicate_num = 1000)
# generate_data("setting12.B1", n = 300, replicate_num = 1000)
# generate_data("setting12.B2", n = 2000, replicate_num = 1000)
# generate_data("setting12.B2", n = 300, replicate_num = 1000)

#------------------------------------------------------------------------------#
# Step 2: Run scenarios 7-1, 7-2, 8-1, 8-2, 9-1, 9-2
# Summarize each scenario immediately after completion
#------------------------------------------------------------------------------#

scenarios_to_run <- c("7-1", "7-2")
# scenarios_to_run <- c("7-3", "7-4")
settings_to_run <- c("setting11", "setting12")

for (scenario_id in scenarios_to_run) {
  for (setting in settings_to_run) {
    # Run the scenario
    run_scenario(scenario_id,
                 setting = setting,
                 version = "test5")

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
      version = "test5"
    )

    print(summary_table)
  }
}
