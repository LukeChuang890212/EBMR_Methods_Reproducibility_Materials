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
    h_x_names_list = rep(list(H_X_NAMES$full), 3)
  ),

  # Scenario 7-alt1: z1 interactions (like 7) with h_x u1u2_z1 for first model
  `7-alt1` = create_ps_spec(
    formulas = list(FORMULAS$full, FORMULAS$u1_z1, FORMULAS$u2_z1),
    h_x_names_list = list(H_X_NAMES$u1u2_z1, H_X_NAMES$full, H_X_NAMES$full)
  ),

  # Scenario 7-alt2: z1 interactions (like 7) with h_x u1u2_z2 for first model
  `7-alt2` = create_ps_spec(
    formulas = list(FORMULAS$full, FORMULAS$u1_z1, FORMULAS$u2_z1),
    h_x_names_list = list(H_X_NAMES$u1u2_z2, H_X_NAMES$full, H_X_NAMES$full)
  ),

  # Scenario 8: with z2 interactions
  `8` = create_ps_spec(
    formulas = list(FORMULAS$full, FORMULAS$u1_z2, FORMULAS$u2_z2),
    h_x_names_list = rep(list(H_X_NAMES$full), 3)
  ),

  # Scenario 9: no-y models
  `9` = create_ps_spec(
    formulas = list(FORMULAS$full, FORMULAS$no_y_u1, FORMULAS$no_y_u2),
    h_x_names_list = rep(list(H_X_NAMES$full), 3)
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
