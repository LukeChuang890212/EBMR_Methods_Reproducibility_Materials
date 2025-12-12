#------------------------------------------------------------------------------#
# Propensity Score Specifications and Scenario Registry
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Inverse Link Functions
#------------------------------------------------------------------------------#
INV_LINKS <- list(
  logistic_complement = function(eta) 1 / (1 + exp(eta)),
  logistic = function(eta) exp(eta) / (1 + exp(eta))
)

#------------------------------------------------------------------------------#
# Formula Components (building blocks)
#------------------------------------------------------------------------------#
FORMULAS <- list(
  # Full model
  full = r ~ o(y) + u1 + u2,

  # Single covariate models
  u1_only = r ~ o(y) + u1,
  u2_only = r ~ o(y) + u2,
  y_only = r ~ o(y),
  z1_only = r ~ o(y) + z1,
  z2_only = r ~ o(y) + z2,


  # Two covariate models
  u1_z1 = r ~ o(y) + u1 + z1,
  u1_z2 = r ~ o(y) + u1 + z2,
  u2_z1 = r ~ o(y) + u2 + z1,
  u2_z2 = r ~ o(y) + u2 + z2,

  # No y models
  no_y_u1 = r ~ u1 + z1 + z2,
  no_y_u2 = r ~ u2 + z1 + z2,

  # Cho2025 specific
  cho_x1x2 = r ~ x1 + x2,
  cho_y_x1 = r ~ o(y) + x1,
  cho_y_x2 = r ~ o(y) + x2
)

#------------------------------------------------------------------------------#
# H_X Names Sets (auxiliary variables for estimation)
#------------------------------------------------------------------------------#
H_X_NAMES <- list(
  full = c("u1", "u2", "z1", "z2"),
  u1_subset = c("u1", "z1", "z2"),
  u2_subset = c("u2", "z1", "z2"),
  z_only = c("z1", "z2"),
  u1u2_z1 = c("u1", "u2", "z1"),
  u1u2_z2 = c("u1", "u2", "z2"),
  with_v1 = c("u1", "u2", "z1", "z2", "v1"),
  with_v2 = c("u1", "u2", "z1", "z2", "v2"),
  with_v3 = c("u1", "u2", "v3"),
  with_v4 = c("u1", "u2", "v4"),
  cho = c("x1", "x2")
)

#------------------------------------------------------------------------------#
# Factory Function: Create PS Specification
#------------------------------------------------------------------------------#
create_ps_spec <- function(formulas, h_x_names_list, inv_link = "logistic_complement") {
  # Resolve inv_link if it's a string

  inv_link_fn <- if (is.character(inv_link)) INV_LINKS[[inv_link]] else inv_link

  list(
    formula.list = formulas,
    h_x_names.list = h_x_names_list,
    inv_link = inv_link_fn
  )
}

#------------------------------------------------------------------------------#
# Pre-defined PS Specifications
# (These replace the repetitive full_ps_specifications blocks)
#------------------------------------------------------------------------------#
PS_SPECS <- list(
  # Scenario 1-1: full -> u1_only -> y_only
  `1-1` = create_ps_spec(
    formulas = list(FORMULAS$full, FORMULAS$u1_only, FORMULAS$y_only),
    h_x_names_list = list(H_X_NAMES$full, H_X_NAMES$u1_subset, H_X_NAMES$z_only)
  ),

  # Scenario 1-2: full -> u2_only -> y_only
  `1-2` = create_ps_spec(
    formulas = list(FORMULAS$full, FORMULAS$u2_only, FORMULAS$y_only),
    h_x_names_list = list(H_X_NAMES$full, H_X_NAMES$u2_subset, H_X_NAMES$z_only)
  ),

  # Scenario 1-3: full -> u1_only -> u2_only
  `1-3` = create_ps_spec(
    formulas = list(FORMULAS$full, FORMULAS$u1_only, FORMULAS$u2_only),
    h_x_names_list = list(H_X_NAMES$full, H_X_NAMES$u1_subset, H_X_NAMES$u2_subset)
  ),

  # Scenario 1-4: 5 models with z variants
  `1-4` = create_ps_spec(
    formulas = list(FORMULAS$full, FORMULAS$u1_only, FORMULAS$u2_only,
                    FORMULAS$z1_only, FORMULAS$z2_only),
    h_x_names_list = rep(list(H_X_NAMES$full), 5)
  ),

  # Scenario 2-1: misspecified h_x (z1 only in first)
  `2-1` = create_ps_spec(
    formulas = list(FORMULAS$full, FORMULAS$u1_only, FORMULAS$y_only),
    h_x_names_list = list(H_X_NAMES$u1u2_z1, H_X_NAMES$u1_subset, H_X_NAMES$z_only)
  ),

  # Scenario 2-2: misspecified h_x (z2 only in first)
  `2-2` = create_ps_spec(
    formulas = list(FORMULAS$full, FORMULAS$u1_only, FORMULAS$y_only),
    h_x_names_list = list(H_X_NAMES$u1u2_z2, H_X_NAMES$u1_subset, H_X_NAMES$z_only)
  ),

  # Scenario 2-3: same as 2-2
  `2-3` = create_ps_spec(
    formulas = list(FORMULAS$full, FORMULAS$u1_only, FORMULAS$y_only),
    h_x_names_list = list(H_X_NAMES$u1u2_z2, H_X_NAMES$u1_subset, H_X_NAMES$z_only)
  ),

  # Scenario 3/4: u2 variant (same PS spec, different model type)
  `3` = create_ps_spec(
    formulas = list(FORMULAS$full, FORMULAS$u2_only, FORMULAS$y_only),
    h_x_names_list = list(H_X_NAMES$full, H_X_NAMES$u2_subset, H_X_NAMES$z_only)
  ),

  # Scenario 5/6: u1/u2 combination
  `5` = create_ps_spec(
    formulas = list(FORMULAS$full, FORMULAS$u1_only, FORMULAS$u2_only),
    h_x_names_list = list(H_X_NAMES$full, H_X_NAMES$u1_subset, H_X_NAMES$u2_subset)
  ),

  # Scenario 7: with z1 interactions
  `7` = create_ps_spec(
    formulas = list(FORMULAS$full, FORMULAS$u1_z1, FORMULAS$u2_z1),
    h_x_names_list = list(H_X_NAMES$full, H_X_NAMES$u1_subset, H_X_NAMES$u2_subset)
  ),

  # Scenario 7-alt1: z1 interactions (like 7) with h_x u1u2_z1 for first model
  `7-alt1` = create_ps_spec(
    formulas = list(FORMULAS$full, FORMULAS$u1_z1, FORMULAS$u2_z1),
    h_x_names_list = list(H_X_NAMES$u1u2_z1, H_X_NAMES$u1_subset, H_X_NAMES$u2_subset)
  ),

  # Scenario 7-alt2: z1 interactions (like 7) with h_x u1u2_z2 for first model
  `7-alt2` = create_ps_spec(
    formulas = list(FORMULAS$full, FORMULAS$u1_z1, FORMULAS$u2_z1),
    h_x_names_list = list(H_X_NAMES$u1u2_z2, H_X_NAMES$u1_subset, H_X_NAMES$u2_subset)
  ),

  # Scenario 8: with z2 interactions
  `8` = create_ps_spec(
    formulas = list(FORMULAS$full, FORMULAS$u1_z2, FORMULAS$u2_z2),
    h_x_names_list = list(H_X_NAMES$full, H_X_NAMES$u1_subset, H_X_NAMES$u2_subset)
  ),

  # Scenario 9: no-y models
  `9` = create_ps_spec(
    formulas = list(FORMULAS$full, FORMULAS$u1_z1, FORMULAS$u2_z2),
    h_x_names_list = list(H_X_NAMES$full, H_X_NAMES$u1_subset, H_X_NAMES$u2_subset)
  ),

  # Scenario 9-10 (setting5/6): with v1/v2
  `9-v` = create_ps_spec(
    formulas = list(FORMULAS$full, FORMULAS$full, FORMULAS$full),
    h_x_names_list = list(H_X_NAMES$full, H_X_NAMES$with_v1, H_X_NAMES$with_v2)
  ),

  # Scenario 11-12 (setting5/6): with v3/v4
  `11-v` = create_ps_spec(
    formulas = list(FORMULAS$full, FORMULAS$full, FORMULAS$full),
    h_x_names_list = list(H_X_NAMES$full, H_X_NAMES$with_v3, H_X_NAMES$with_v4)
  ),

  # Cho2025 scenarios
  cho1 = create_ps_spec(
    formulas = list(FORMULAS$cho_x1x2, FORMULAS$cho_y_x1, FORMULAS$cho_y_x2),
    h_x_names_list = rep(list(H_X_NAMES$cho), 3),
    inv_link = "logistic"
  )
)

#------------------------------------------------------------------------------#
# Factory Function: Create Scenario
#------------------------------------------------------------------------------#
create_scenario <- function(
  id,
  description,
  ps_spec_id,
  settings,
  missing_rates = c("miss50", "miss30"),
  model_type = "correct",
  n_vector = NULL,
  version = NULL,
  enabled = TRUE
) {
  list(
    id = id,
    description = description,
    ps_spec_id = ps_spec_id,
    settings = settings,
    missing_rates = missing_rates,
    model_type = model_type,
    n_vector = n_vector,
    version = version,
    enabled = enabled
  )
}

#------------------------------------------------------------------------------#
# Scenario Registry
#------------------------------------------------------------------------------#
SCENARIOS <- list(
  # Main paper scenarios (settings 11, 12 - continuous and binary outcomes)

  `1-1` = create_scenario(
    id = "1-1",
    description = "Correct model: full -> u1 -> y_only",
    ps_spec_id = "1-1",
    settings = c("setting11", "setting12"),
    model_type = "correct",
    version = "test8"
  ),

  `1-2` = create_scenario(
    id = "1-2",
    description = "Correct model: full -> u2 -> y_only",
    ps_spec_id = "1-2",
    settings = c("setting11", "setting12"),
    model_type = "correct",
    version = "test8"
  ),

  `1-3` = create_scenario(
    id = "1-3",
    description = "Correct model: full -> u1 -> u2",
    ps_spec_id = "1-3",
    settings = c("setting11", "setting12"),
    model_type = "correct",
    version = "test8"
  ),

  `1-4` = create_scenario(
    id = "1-4",
    description = "Correct model: 5 models with z variants",
    ps_spec_id = "1-4",
    settings = c("setting11", "setting12"),
    model_type = "correct",
    version = "test6"
  ),

  `2-1` = create_scenario(
    id = "2-1",
    description = "Misspecified: h_x uses z1 only",
    ps_spec_id = "2-1",
    settings = c("setting11", "setting12"),
    model_type = "misspecified",
    version = "test7"
  ),

  `2-2` = create_scenario(
    id = "2-2",
    description = "Misspecified: h_x uses z2 only",
    ps_spec_id = "2-2",
    settings = c("setting11", "setting12"),
    model_type = "misspecified",
    version = "test7"
  ),

  `2-3` = create_scenario(
    id = "2-3",
    description = "Misspecified: variant 3",
    ps_spec_id = "2-3",
    settings = c("setting11", "setting12"),
    model_type = "misspecified",
    version = "test7"
  ),

  `3` = create_scenario(
    id = "3",
    description = "Correct model: u2 variant",
    ps_spec_id = "3",
    settings = c("setting11", "setting12"),
    model_type = "correct",
    version = "test5"
  ),

  `4` = create_scenario(
    id = "4",
    description = "Misspecified: u2 variant",
    ps_spec_id = "3",
    settings = c("setting11", "setting12"),
    model_type = "misspecified",
    version = "test5"
  ),

  `5` = create_scenario(
    id = "5",
    description = "Correct model: u1/u2 combination",
    ps_spec_id = "5",
    settings = c("setting11", "setting12"),
    model_type = "correct",
    version = "test5"
  ),

  `6` = create_scenario(
    id = "6",
    description = "Misspecified: u1/u2 combination",
    ps_spec_id = "5",
    settings = c("setting11", "setting12"),
    model_type = "misspecified",
    version = "test5"
  ),

  `7-1` = create_scenario(
    id = "7-1",
    description = "Correct model: z1 interactions",
    ps_spec_id = "7",
    settings = c("setting11", "setting12"),
    model_type = "correct",
    version = "test5"
  ),

  `7-2` = create_scenario(
    id = "7-2",
    description = "Misspecified: z1 interactions",
    ps_spec_id = "7",
    settings = c("setting11", "setting12"),
    model_type = "misspecified",
    version = "test5"
  ),

  `8-1` = create_scenario(
    id = "8-1",
    description = "Correct model: z1/z2 interactions",
    ps_spec_id = "8",
    settings = c("setting11", "setting12"),
    model_type = "correct",
    version = "test5"
  ),

  `8-2` = create_scenario(
    id = "8-2",
    description = "Misspecified: z1/z2 interactions",
    ps_spec_id = "8",
    settings = c("setting11", "setting12"),
    model_type = "misspecified",
    version = "test5"
  ),

  `9-1` = create_scenario(
    id = "9-1",
    description = "Correct model: no-y models",
    ps_spec_id = "9",
    settings = c("setting11", "setting12"),
    model_type = "correct",
    version = "test5"
  ),

  `9-2` = create_scenario(
    id = "9-2",
    description = "Misspecified: no-y models",
    ps_spec_id = "9",
    settings = c("setting11", "setting12"),
    model_type = "misspecified",
    version = "test5"
  ),

  `7-3` = create_scenario(
    id = "7-3",
    description = "Misspecified: z1 interactions with h_x u1u2_z1 for first model",
    ps_spec_id = "7-alt1",
    settings = c("setting11", "setting12"),
    model_type = "misspecified",
    version = "test_v1"
  ),

  `7-4` = create_scenario(
    id = "7-4",
    description = "Misspecified: z1 interactions with h_x u1u2_z2 for first model",
    ps_spec_id = "7-alt2",
    settings = c("setting11", "setting12"),
    model_type = "misspecified",
    version = "test_v1"
  ),

  # Setting 5/6 scenarios (with auxiliary variables v)

  `9` = create_scenario(
    id = "9",
    description = "Setting 5/6: with v1/v2 (correct)",
    ps_spec_id = "9-v",
    settings = c("setting5", "setting6"),
    model_type = "correct",
    version = "test5"
  ),

  `10` = create_scenario(
    id = "10",
    description = "Setting 5/6: with v1/v2 (misspecified)",
    ps_spec_id = "9-v",
    settings = c("setting5", "setting6"),
    model_type = "misspecified",
    version = "test5"
  ),

  `11` = create_scenario(
    id = "11",
    description = "Setting 5/6: with v3/v4 (correct)",
    ps_spec_id = "11-v",
    settings = c("setting5", "setting6"),
    model_type = "correct",
    version = "test5"
  ),

  `12` = create_scenario(
    id = "12",
    description = "Setting 5/6: with v3/v4 (misspecified)",
    ps_spec_id = "11-v",
    settings = c("setting5", "setting6"),
    model_type = "misspecified",
    version = "test5"
  ),

  # Cho2025 scenarios

  cho1 = create_scenario(
    id = "cho1",
    description = "Cho2025: RM2/RM3 variants",
    ps_spec_id = "cho1",
    settings = c("Cho_RM2q", "Cho_RM2p", "Cho_RM3q", "Cho_RM3p"),
    missing_rates = c("miss30", "miss50"),
    model_type = "correct",
    n_vector = c(1000, 4000),
    version = "test7"
  ),

  cho2 = create_scenario(
    id = "cho2",
    description = "Cho2025: RM2/RM3 misspecified",
    ps_spec_id = "cho1",
    settings = c("Cho_RM2", "Cho_RM3"),
    model_type = "misspecified",
    version = "test5"
  )
)

#------------------------------------------------------------------------------#
# Helper Functions
#------------------------------------------------------------------------------#

#' Get PS specification by ID
get_ps_spec <- function(spec_id) {
  if (!spec_id %in% names(PS_SPECS)) {
    stop(paste("Unknown PS specification:", spec_id,
               "\nAvailable:", paste(names(PS_SPECS), collapse = ", ")))
  }
  PS_SPECS[[spec_id]]
}

#' Get scenario by ID
get_scenario <- function(scenario_id) {
  if (!scenario_id %in% names(SCENARIOS)) {
    stop(paste("Unknown scenario:", scenario_id,
               "\nAvailable:", paste(names(SCENARIOS), collapse = ", ")))
  }
  SCENARIOS[[scenario_id]]
}

#' List all scenarios as a data frame
list_scenarios <- function(enabled_only = FALSE) {
  scenarios_to_list <- if (enabled_only) {
    Filter(function(s) s$enabled, SCENARIOS)
  } else {
    SCENARIOS
  }

  df <- do.call(rbind, lapply(names(scenarios_to_list), function(id) {
    s <- scenarios_to_list[[id]]
    data.frame(
      id = s$id,
      description = s$description,
      settings = paste(s$settings, collapse = ", "),
      model_type = s$model_type,
      version = ifelse(is.null(s$version), "default", s$version),
      enabled = s$enabled,
      stringsAsFactors = FALSE
    )
  }))
  df
}

#' Get all enabled scenario IDs
get_enabled_scenarios <- function() {
  names(Filter(function(s) s$enabled, SCENARIOS))
}

#' Enable/disable a scenario
set_scenario_enabled <- function(scenario_id, enabled = TRUE) {
  if (!scenario_id %in% names(SCENARIOS)) {
    stop(paste("Unknown scenario:", scenario_id))
  }
  SCENARIOS[[scenario_id]]$enabled <<- enabled
  invisible(enabled)
}

#' Get model parameters based on model type
#'
#' @param model_type Either "correct" or "misspecified"
#' @return List with n_vector, data_files, and alpha_true
get_model_params <- function(model_type) {
  if (model_type == "correct") {
    list(
      n_vector = n.vector.list$correct_model,
      data_files = correct_model_all_data_file.list,
      alpha_true = correct_model_alpha.true.list
    )
  } else if (model_type == "misspecified") {
    list(
      n_vector = n.vector.list$misspecified_model,
      data_files = misspecified_model_all_data_file.list,
      alpha_true = misspecified_model_alpha.true.list
    )
  } else {
    stop(paste("Unknown model type:", model_type))
  }
}

#------------------------------------------------------------------------------#
# Configuration Object
#------------------------------------------------------------------------------#
CONFIG <- list(
  replicate_num = 1000,
  n_default = 2000
)

#' Show available scenarios (pretty print)
show_scenarios <- function() {
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("Available Scenarios\n")
  cat(strrep("=", 70), "\n\n")

  df <- list_scenarios()
  for (i in 1:nrow(df)) {
    cat(sprintf("%-8s %s\n", df$id[i], df$description[i]))
    cat(sprintf("         Settings: %s | Model: %s | Version: %s\n",
                df$settings[i], df$model_type[i], df$version[i]))
    cat("\n")
  }
}

#' Show current configuration
show_config <- function() {
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("Current Configuration\n")
  cat(strrep("=", 70), "\n\n")
  cat("  Replicate number:", CONFIG$replicate_num, "\n")
  cat("  Default sample size:", CONFIG$n_default, "\n")
  cat("\n")
}

#' Run a scenario
#'
#' @param scenario_id The scenario ID to run (e.g., "7-1", "9-1")
#' @param setting Which setting to run (e.g., "setting11", "setting12")
#' @param n Sample size (optional - if NULL, runs all sizes from n.vector.list)
#' @param replicate_num Number of replicates (default: from CONFIG or 1000)
#' @param version Version string for output files
#' @param missing_rates Which missing rates to run (default: from scenario)
run_scenario <- function(scenario_id,
                         setting = NULL,
                         n = NULL,
                         replicate_num = NULL,
                         version = NULL,
                         missing_rates = NULL) {

  # Get scenario configuration
  scenario <- get_scenario(scenario_id)
  ps_spec <- get_ps_spec(scenario$ps_spec_id)

  # Use defaults if not provided
  if (is.null(setting)) setting <- scenario$settings[1]
  if (is.null(replicate_num)) replicate_num <- if (!is.null(CONFIG$replicate_num)) CONFIG$replicate_num else 1000
  if (is.null(version)) version <- scenario$version
  if (is.null(missing_rates)) missing_rates <- scenario$missing_rates

  # Get model parameters
  params <- get_model_params(scenario$model_type)

  # Get n.vector from params - this contains all sample sizes to run
  n_vector_list <- params$n_vector

  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("Running Scenario:", scenario_id, "\n")
  cat("Description:", scenario$description, "\n")
  cat(strrep("=", 70), "\n")
  cat("  Setting:", setting, "\n")
  cat("  Replicates:", replicate_num, "\n")
  cat("  Version:", version, "\n")
  cat("  Missing rates:", paste(missing_rates, collapse = ", "), "\n")
  cat("\n")

  # Get data file config and alpha.true for this setting
  data_file_config <- params$data_files[[setting]]
  if (is.null(data_file_config)) {
    stop(paste("No data file configuration for setting:", setting))
  }
  alpha.true_config <- params$alpha_true[[setting]]

  # Run for each missing rate
  for (miss_rate in missing_rates) {
    cat("Processing:", miss_rate, "\n")

    # Get data file path(s) - the config stores full paths in a list
    data_files_for_rate <- data_file_config[[miss_rate]]
    if (is.null(data_files_for_rate)) {
      cat("  WARNING: No data file config for", miss_rate, "\n")
      next
    }

    # Get alpha.true for this missing rate
    alpha.true_for_rate <- alpha.true_config[[miss_rate]]

    # Iterate over data files (each may correspond to different n values)
    for (file_idx in seq_along(data_files_for_rate)) {
      data_file <- data_files_for_rate[[file_idx]]
      alpha.true <- alpha.true_for_rate[[file_idx]]

      # Get n.vector for this file index
      n_vector <- n_vector_list[[file_idx]]

      # If user specified n, filter to just that n
      if (!is.null(n)) {
        if (n %in% n_vector) {
          n_vector <- n
        } else {
          next  # Skip this file if it doesn't contain the requested n
        }
      }

      if (!file.exists(data_file)) {
        cat("  WARNING: Data file not found:", data_file, "\n")
        cat("  Use generate_data() to create it first.\n")
        next
      }

      # Load data
      all_data <- readRDS(data_file)
      cat("  Loaded data:", data_file, "\n")

      # Run for each sample size in this n_vector
      for (current_n in n_vector) {
        cat("  Sample size n =", current_n, "\n")

        # Run simulation for each model combination
        J <- length(ps_spec$formula.list)
        for (model_num in 1:J) {
          model_combinations <- combn(J, model_num)
          for (i in 1:ncol(model_combinations)) {
            model_set <- model_combinations[, i]
            model_str <- paste0(model_set, collapse = "")

            cat("    Running models:", model_str, "\n")

            # Create subset PS specification
            subset_ps_spec <- list(
              formula.list = ps_spec$formula.list[model_set],
              h_x_names.list = ps_spec$h_x_names.list[model_set],
              inv_link = ps_spec$inv_link
            )

            # Output file
            save_file <- paste0(
              "Simulation_Results/EBMR_IPW_", setting, "-", miss_rate,
              "-scenario", scenario_id, "_", model_str,
              "_n", current_n, "_replicate", replicate_num, "_", version, ".RDS"
            )

            # Check if already exists
            if (file.exists(save_file)) {
              cat("      Already exists, skipping\n")
              next
            }

            # Define true PS model function based on setting
            ps_model.true <- if (setting %in% c("Cho_RM2", "Cho_RM2p", "Cho_RM2q")) {
              function(dat, alpha.true) {
                eta <- cbind(rep(1, nrow(dat)), dat$x1, dat$y) %*% alpha.true
                exp(eta) / (1 + exp(eta))
              }
            } else if (setting %in% c("Cho_RM3", "Cho_RM3p", "Cho_RM3q")) {
              function(dat, alpha.true) {
                eta <- cbind(rep(1, nrow(dat)), dat$x2, dat$y) %*% alpha.true
                exp(eta) / (1 + exp(eta))
              }
            } else {
              function(dat, alpha.true) {
                1 / (1 + exp(cbind(rep(1, nrow(dat)), dat$y, dat$u1, dat$u2) %*% alpha.true))
              }
            }

            # Run simulation
            simulate(
              all_data = all_data,
              ps_model.true = ps_model.true,
              alpha.true = alpha.true,
              ps_specifications = subset_ps_spec,
              n = current_n,
              replicate_num = replicate_num,
              save_file = save_file
            )

            # Print console summary after simulation completes
            # Uses clean_sim_result() from Simulation.r for consistent logic
            if (file.exists(save_file)) {
              sim_result <- readRDS(save_file)
              cleaned <- clean_sim_result(sim_result, multiplier = 30, verbose = FALSE)
              cat("      Summary: Total=", cleaned$n_total, ", Successful=", cleaned$n_successful,
                  ", NA=", cleaned$n_na, ", Outliers=", cleaned$n_outliers, "\n", sep = "")
            }
          }
        }
      }
    }
  }

  cat("\nScenario", scenario_id, "completed!\n")
}
