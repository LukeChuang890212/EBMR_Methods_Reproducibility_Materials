#------------------------------------------------------------------------------#
# Propensity Score Specifications and Scenario Registry
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Equivalent Scenario Mapping for Model 2 and Model 3
#
# When model combinations exclude model 1 (e.g., "2", "3", "23"), different
# scenarios may produce identical results if they share the same model 2 and 3
# specifications. This mapping defines which scenarios are equivalent for
# non-model-1 combinations, allowing us to copy results instead of re-running.
#
# Key = current scenario, Value = list of equivalent scenarios (priority order)
#------------------------------------------------------------------------------#
# =============================================================================
# IMPORTANT: Equivalence requires SAME DATA (same model_type: correct or misspecified)
# - correct model (1-1, 1-3, 1-4, 2-1, 3-1, 7-1, 8-1, 9-1) uses A data
# - misspecified model (1-2, 4, 7-2, 7-3, 7-4, 8-2, ..., 9-2) uses B data
# =============================================================================

# Model 1 equivalence: scenarios sharing same model 1 (formula + h_alpha) AND same data
# Correct (A data): {7-1, 8-1, 9-1} — full + full h_alpha
#                    {1-1, 2-1, 3-1} — y_u1 + full h_alpha
#                    1-3 (full_sq1), 1-4 (full_sq2) each unique h_alpha
# Misspecified (B data): {7-2, 8-2, 9-2} — full + full h_alpha
#                        {1-2, 4} — y_u1 + full h_alpha
#                        {7-3, 8-3, 9-3} — full + full_sq1 h_alpha
#                        {7-4, 8-4, 9-4} — full + full_sq2 h_alpha
EQUIVALENT_SCENARIOS_MODEL1 <- list(
  `1-1` = c("2-1", "3-1"),
  `2-1` = c("1-1", "3-1"),
  `3-1` = c("1-1", "2-1"),
  `1-2` = c("4"),
  `4` = c("1-2"),
  `7-1` = c("8-1", "9-1"),
  `8-1` = c("7-1", "9-1"),
  `9-1` = c("7-1", "8-1"),
  `7-2` = c("8-2", "9-2"),
  `8-2` = c("7-2", "9-2"),
  `9-2` = c("7-2", "8-2"),
  `7-3` = c("8-3", "9-3"),
  `8-3` = c("7-3", "9-3"),
  `9-3` = c("7-3", "8-3"),
  `7-4` = c("8-4", "9-4"),
  `8-4` = c("7-4", "9-4"),
  `9-4` = c("7-4", "8-4")
)

# Model 2 equivalence: scenarios sharing same model 2 (formula + h_alpha) AND same data
# Correct (A data): {1-1, 1-3, 1-4, 3-1} — y_z1 + full h_alpha
#                    {7-1, 8-1} — u1_z1 + full h_alpha. 9-1 has u1_z2. 2-1 has y_z1_z2.
# Misspecified (B data): {1-2, 4} — y_z1 + full h_alpha
#                        {7-2, 7-3, 7-4, 8-2, 8-3, 8-4} — all u1_z1 + full h_alpha
#                        {9-2, 9-3, 9-4} — all u1_z2 + full h_alpha
EQUIVALENT_SCENARIOS_MODEL2 <- list(
  `1-1` = c("1-3", "1-4", "3-1"),
  `1-3` = c("1-1", "1-4", "3-1"),
  `1-4` = c("1-1", "1-3", "3-1"),
  `3-1` = c("1-1", "1-3", "1-4"),
  `1-2` = c("4"),
  `4` = c("1-2"),
  `7-1` = c("8-1"),
  `8-1` = c("7-1"),
  `7-2` = c("7-3", "7-4", "8-2", "8-3", "8-4"),
  `7-3` = c("7-2", "7-4", "8-2", "8-3", "8-4"),
  `7-4` = c("7-2", "7-3", "8-2", "8-3", "8-4"),
  `8-2` = c("7-2", "7-3", "7-4", "8-3", "8-4"),
  `8-3` = c("7-2", "7-3", "7-4", "8-2", "8-4"),
  `8-4` = c("7-2", "7-3", "7-4", "8-2", "8-3"),
  `9-2` = c("9-3", "9-4"),
  `9-3` = c("9-2", "9-4"),
  `9-4` = c("9-2", "9-3")
)

# Model 3 equivalence: scenarios sharing same model 3 (formula + h_alpha) AND same data
# Correct (A data): {1-1, 1-3, 1-4, 2-1} — y_z2 + full h_alpha
#                    {8-1, 9-1} — u2_z2 + full h_alpha. 7-1 has u2_z1. 3-1 has y_z2_sq.
# Misspecified (B data): 1-2 — y_z2 + full (alone). 4 — y_z2_sq + full (alone).
#                        {7-2, 7-3, 7-4} — u2_z1. {8-2, 8-3, 8-4, 9-2, 9-3, 9-4} — u2_z2.
EQUIVALENT_SCENARIOS_MODEL3 <- list(
  `1-1` = c("1-3", "1-4", "2-1"),
  `1-3` = c("1-1", "1-4", "2-1"),
  `1-4` = c("1-1", "1-3", "2-1"),
  `2-1` = c("1-1", "1-3", "1-4"),
  `8-1` = c("9-1"),
  `9-1` = c("8-1"),
  `7-2` = c("7-3", "7-4"),
  `7-3` = c("7-2", "7-4"),
  `7-4` = c("7-2", "7-3"),
  `8-2` = c("8-3", "8-4", "9-2", "9-3", "9-4"),
  `8-3` = c("8-2", "8-4", "9-2", "9-3", "9-4"),
  `8-4` = c("8-2", "8-3", "9-2", "9-3", "9-4"),
  `9-2` = c("8-2", "8-3", "8-4", "9-3", "9-4"),
  `9-3` = c("8-2", "8-3", "8-4", "9-2", "9-4"),
  `9-4` = c("8-2", "8-3", "8-4", "9-2", "9-3")
)

# Model 1+2 equivalence: same model 1 AND model 2 AND same data
# Correct (A data): {1-1, 3-1} — y_u1+full for M1, y_z1+full for M2
#                    {7-1, 8-1} — full+full for M1, u1_z1+full for M2
#   1-3/1-4 differ in M1 h_alpha; 2-1 differs in M2
# Misspecified (B data): {1-2, 4} — y_u1+full for M1, y_z1+full for M2
#                        {7-2, 8-2} — full+full for M1, u1_z1+full for M2
EQUIVALENT_SCENARIOS_MODEL12 <- list(
  `1-1` = c("3-1"),
  `3-1` = c("1-1"),
  `1-2` = c("4"),
  `4` = c("1-2"),
  `7-1` = c("8-1"),
  `8-1` = c("7-1"),
  `7-2` = c("8-2"),
  `8-2` = c("7-2")
)

# Model 1+3 equivalence: same model 1 AND model 3 AND same data
# Correct (A data): {1-1, 2-1} — y_u1+full for M1, y_z2+full for M3
#                    {8-1, 9-1} — full+full for M1, u2_z2+full for M3
#   3-1 differs in M3; 1-3/1-4 differ in M1 h_alpha
# Misspecified (B data): {8-2, 9-2} — full+full for M1, u2_z2+full for M3
#                        {8-3, 9-3} — full+full_sq1 for M1, u2_z2+full for M3
#                        {8-4, 9-4} — full+full_sq2 for M1, u2_z2+full for M3
#   1-2 alone (4 has different M3)
EQUIVALENT_SCENARIOS_MODEL13 <- list(
  `1-1` = c("2-1"),
  `2-1` = c("1-1"),
  `8-1` = c("9-1"),
  `9-1` = c("8-1"),
  `8-2` = c("9-2"),
  `9-2` = c("8-2"),
  `8-3` = c("9-3"),
  `9-3` = c("8-3"),
  `8-4` = c("9-4"),
  `9-4` = c("8-4")
)

# Model 2+3 equivalence: same model 2 AND model 3 AND same data
# Correct (A data): {1-1, 1-3, 1-4} — y_z1+full for M2, y_z2+full for M3
#   2-1 differs in M2; 3-1 differs in M3; 7-1 has different M3 from 8-1
# Misspecified (B data):
#   {7-2, 7-3, 7-4} — all u1_z1+full for M2, u2_z1+full for M3
#   {8-2, 8-3, 8-4} — all u1_z1+full for M2, u2_z2+full for M3
#   {9-2, 9-3, 9-4} — all u1_z2+full for M2, u2_z2+full for M3
#   1-2 alone (4 differs in M3)
EQUIVALENT_SCENARIOS_MODEL23 <- list(
  `1-1` = c("1-3", "1-4"),
  `1-3` = c("1-1", "1-4"),
  `1-4` = c("1-1", "1-3"),
  `7-2` = c("7-3", "7-4"),
  `7-3` = c("7-2", "7-4"),
  `7-4` = c("7-2", "7-3"),
  `8-2` = c("8-3", "8-4"),
  `8-3` = c("8-2", "8-4"),
  `8-4` = c("8-2", "8-3"),
  `9-2` = c("9-3", "9-4"),
  `9-3` = c("9-2", "9-4"),
  `9-4` = c("9-2", "9-3")
)

# Model 1+2+3 equivalence: ALL three models identical AND same data
# No equivalence exists (7-1/8-1/9-1 all differ in at least one model)
# 7-2 is unique among misspecified (7-3, 7-4 have different model 1)
EQUIVALENT_SCENARIOS_MODEL123 <- list(
  # No entries - no scenarios share all 3 models with same data
)

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
  full = r ~ y + u1 + u2,
  full_z1 = r ~ y + u1 + u2 + z1,
  full_z2 = r ~ y + u1 + u2 + z2,
  full_u2sq = r ~ y + u1 + u2 + u2^2,

  # Single covariate models
  u1_only = r ~ y + u1,
  u2_only = r ~ y + u2,
  y_only = r ~ y,
  z1_only = r ~ y + z1,
  z2_only = r ~ y + z2,


  # Two covariate models
  u1_z1 = r ~ y + u1 + z1,
  u1_z2 = r ~ y + u1 + z2,
  u2_z1 = r ~ y + u2 + z1,
  u2_z2 = r ~ y + u2 + z2,

  # Interaction models
  y_u1 = r ~ y + u1 + z1 + y:z1,
  y_u2 = r ~ y + u2 + z1 + y:z1,
  y_z1 = r ~ y + u2 + z1 + y:z1,
  y_z2 = r ~ y + u2 + z2 + y:z2,
  y_z1_z2 = r ~ y + z1 + z2,
  y_z2_sq = r ~ y + z2 + z2^2,

  # No y models
  no_y_u1 = r ~ u1 + z1 + z2,
  no_y_u2 = r ~ u2 + z1 + z2,

  # Cho2025 specific
  cho_x1x2 = r ~ x1 + x2,
  cho_y_x1 = r ~ y + x1,
  cho_y_x2 = r ~ y + x2

)

#------------------------------------------------------------------------------#
# H_ALPHA Definitions (auxiliary variable functions for estimation)
#
# Each entry can be:
#   - A character vector of column names (backwards compatible), e.g., c("u1", "u2")
#   - A function: function(data) -> data.frame/matrix with named columns
#------------------------------------------------------------------------------#
H_ALPHA <- list(
  full = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2,
                             u1_u2 = dat$u1*dat$u2, z1_z2 = dat$z1*dat$z2,
                             u1_z1 = dat$u1*dat$z1, z1_u2 = dat$z1*dat$u2,
                             u1_z2 = dat$u1*dat$z2, u2_z2 = dat$u2*dat$z2),
  # full = c("u1", "u2", "z1", "z2"),
  # full = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2,
  #                            u2_sq = dat$u2^2, z2_sq = dat$z2^2),
  u1_subset = c("u1", "z1", "z2"),
  u2_subset = c("u2", "z1", "z2"),
  z_only = c("z1", "z2"),
  u1u2_z1 = c("u1", "u2", "z1"),
  u1u2_z2 = c("u1", "u2", "z2"),
  full_sq1 = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2,
                                 u1_u2 = dat$u1*dat$u2, z1_z2 = dat$z1*dat$z2,
                                 u1_z1 = dat$u1*dat$z1, z1_u2 = dat$z1*dat$u2,
                                 u1_z2 = dat$u1*dat$z2, u2_z2 = dat$u2*dat$z2, u2_sq = dat$u2^2),
  full_sq2 = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2,
                                 u1_u2 = dat$u1*dat$u2, z1_z2 = dat$z1*dat$z2,
                                 u1_z1 = dat$u1*dat$z1, z1_u2 = dat$z1*dat$u2,
                                 u1_z2 = dat$u1*dat$z2, u2_z2 = dat$u2*dat$z2, z2_sq = dat$z2^2),
  # full_sq1 = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2, u2_sq = dat$u2^2),
  # full_sq2 = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2, z2_sq = dat$z2^2),
  with_v1 = c("u1", "u2", "z1", "z2", "v1"),
  with_v2 = c("u1", "u2", "z1", "z2", "v2"),
  with_v3 = c("u1", "u2", "v3"),
  with_v4 = c("u1", "u2", "v4"),
  cho = c("x1", "x2", "x3")
)

#------------------------------------------------------------------------------#
# Factory Function: Create PS Specification
#------------------------------------------------------------------------------#
create_ps_spec <- function(formulas, h_alpha_list, inv_link = "logistic_complement",
                           outcome = "y") {
  # Resolve inv_link if it's a string
  inv_link_fn <- if (is.character(inv_link)) INV_LINKS[[inv_link]] else inv_link

  list(
    formula.list = formulas,
    h_alpha.list = h_alpha_list,
    inv_link = inv_link_fn,
    outcome = outcome
  )
}

#------------------------------------------------------------------------------#
# Pre-defined PS Specifications
# (These replace the repetitive full_ps_specifications blocks)
#------------------------------------------------------------------------------#
PS_SPECS <- list(
    # PS 1: y-interaction models (y:u1, y:z1, y:z2)
  `1` = create_ps_spec(
    formulas = list(FORMULAS$y_u1, FORMULAS$y_z1, FORMULAS$y_z2),
    h_alpha_list = list(H_ALPHA$full, H_ALPHA$full, H_ALPHA$full)
  ),

  # PS 1-alt1: y-interaction models, h_alpha full_sq1 (full + u2^2) for model 1
  `1-alt1` = create_ps_spec(
    formulas = list(FORMULAS$y_u1, FORMULAS$y_z1, FORMULAS$y_z2),
    h_alpha_list = list(H_ALPHA$full_sq1, H_ALPHA$full, H_ALPHA$full)
  ),

  # PS 1-alt2: y-interaction models, h_alpha full_sq2 (full + z2^2) for model 1
  `1-alt2` = create_ps_spec(
    formulas = list(FORMULAS$y_u1, FORMULAS$y_z1, FORMULAS$y_z2),
    h_alpha_list = list(H_ALPHA$full_sq2, H_ALPHA$full, H_ALPHA$full)
  ),

  # PS 2: y-interaction models (y:u1, z1:z2, y:z2)
  `2` = create_ps_spec(
    formulas = list(FORMULAS$y_u1, FORMULAS$y_z1_z2, FORMULAS$y_z2),
    h_alpha_list = list(H_ALPHA$full, H_ALPHA$full, H_ALPHA$full)
  ),

  # PS 3: y-interaction models (y:u1, y:z1, z2^2)
  `3` = create_ps_spec(
    formulas = list(FORMULAS$y_u1, FORMULAS$y_z1, FORMULAS$y_z2_sq),
    h_alpha_list = list(H_ALPHA$full, H_ALPHA$full, H_ALPHA$full)
  ),

  # PS 7: z1 subset models (full, u1_z1, u2_z1)
  `7` = create_ps_spec(
    formulas = list(FORMULAS$full, FORMULAS$u1_z1, FORMULAS$u2_z1),
    h_alpha_list = list(H_ALPHA$full, H_ALPHA$full, H_ALPHA$full)
  ),

  # PS 7-alt1: z1 interactions (like 7), h_alpha full_sq1 (full + u2^2) for model 1
  `7-alt1` = create_ps_spec(
    formulas = list(FORMULAS$full, FORMULAS$u1_z1, FORMULAS$u2_z1),
    h_alpha_list = list(H_ALPHA$full_sq1, H_ALPHA$full, H_ALPHA$full)
  ),

  # PS 7-alt2: z1 interactions (like 7), h_alpha full_sq2 (full + z2^2) for model 1
  `7-alt2` = create_ps_spec(
    formulas = list(FORMULAS$full, FORMULAS$u1_z1, FORMULAS$u2_z1),
    h_alpha_list = list(H_ALPHA$full_sq2, H_ALPHA$full, H_ALPHA$full)
  ),

  # PS 8: z1/z2 mixed models (full, u1_z1, u2_z2)
  `8` = create_ps_spec(
    formulas = list(FORMULAS$full, FORMULAS$u1_z1, FORMULAS$u2_z2),
    h_alpha_list = list(H_ALPHA$full, H_ALPHA$full, H_ALPHA$full)
  ),

  # PS 8-alt1: expanded models (full+u2^2, full+z1, full+z2), h_alpha full for model 1
  `8-alt1` = create_ps_spec(
    formulas = list(FORMULAS$full, FORMULAS$u1_z1, FORMULAS$u2_z2),
    h_alpha_list = list(H_ALPHA$full_sq1, H_ALPHA$full, H_ALPHA$full)
  ),

  # PS 8-alt2: expanded models (full+u2^2, full+z1, full+z2), h_alpha full_sq1 (full + u2^2) for model 1
  `8-alt2` = create_ps_spec(
    formulas = list(FORMULAS$full, FORMULAS$u1_z1, FORMULAS$u2_z2),
    h_alpha_list = list(H_ALPHA$full_sq2, H_ALPHA$full, H_ALPHA$full)
  ),

  # PS 9: z2 subset models (full, u1_z2, u2_z2)
  `9` = create_ps_spec(
    formulas = list(FORMULAS$full, FORMULAS$u1_z2, FORMULAS$u2_z2),
    h_alpha_list = list(H_ALPHA$full, H_ALPHA$full, H_ALPHA$full)
  ),
  
  # PS 9-alt1: expanded models (full, u1_z2, u2_z2), h_alpha full for model 1
  `9-alt1` = create_ps_spec(
    formulas = list(FORMULAS$full, FORMULAS$u1_z2, FORMULAS$u2_z2),
    h_alpha_list = list(H_ALPHA$full_sq1, H_ALPHA$full, H_ALPHA$full)
  ),
  
  # PS 9-alt2: expanded models (full, u1_z2, u2_z2), h_alpha full_sq1 (full + u2^2) for model 1
  `9-alt2` = create_ps_spec(
    formulas = list(FORMULAS$full, FORMULAS$u1_z2, FORMULAS$u2_z2),
    h_alpha_list = list(H_ALPHA$full_sq2, H_ALPHA$full, H_ALPHA$full)
  ),
  
  # Scenario 9-10 (setting5/6): with v1/v2
  `9-v` = create_ps_spec(
    formulas = list(FORMULAS$full, FORMULAS$full, FORMULAS$full),
    h_alpha_list = list(H_ALPHA$full, H_ALPHA$with_v1, H_ALPHA$with_v2)
  ),

  # Scenario 11-12 (setting5/6): with v3/v4
  `11-v` = create_ps_spec(
    formulas = list(FORMULAS$full, FORMULAS$full, FORMULAS$full),
    h_alpha_list = list(H_ALPHA$full, H_ALPHA$with_v3, H_ALPHA$with_v4)
  ),

  # Cho2025 scenarios
  cho1 = create_ps_spec(
    formulas = list(FORMULAS$cho_y_x1, FORMULAS$cho_y_x2, FORMULAS$cho_x1x2),
    h_alpha_list = rep(list(H_ALPHA$cho), 3),
    inv_link = "logistic"
  ),
  
  cho2 = create_ps_spec(
    formulas = list(FORMULAS$cho_y_x2, FORMULAS$cho_y_x1, FORMULAS$cho_x1x2),
    h_alpha_list = rep(list(H_ALPHA$cho), 3),
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
    description = "Correct: y-interaction models (y:u1, y:z1, y:z2)",
    ps_spec_id = "1",
    settings = c("setting11", "setting12"),
    model_type = "correct",
    version = "test8"
  ),

  `1-2` = create_scenario(
    id = "1-2",
    description = "Misspecified: y-interaction models (y:u1, y:z1, y:z2)",
    ps_spec_id = "1",
    settings = c("setting11", "setting12"),
    model_type = "misspecified",
    version = "test8"
  ),

  `1-3` = create_scenario(
    id = "1-3",
    description = "Correct: y-interactions, h_alpha full_sq1 for M1",
    ps_spec_id = "1-alt1",
    settings = c("setting11", "setting12"),
    model_type = "correct",
    version = "test8"
  ),

  `1-4` = create_scenario(
    id = "1-4",
    description = "Correct: y-interactions, h_alpha full_sq2 for M1",
    ps_spec_id = "1-alt2",
    settings = c("setting11", "setting12"),
    model_type = "correct",
    version = "test6"
  ),

  `2-1` = create_scenario(
    id = "2-1",
    description = "Correct: y:u1, z1:z2, y:z2 interaction models",
    ps_spec_id = "2",
    settings = c("setting11", "setting12"),
    model_type = "correct",
    version = "test7"
  ),

  `3-1` = create_scenario(
    id = "3-1",
    description = "Correct: y:u1, y:z1, z2^2 models",
    ps_spec_id = "3",
    settings = c("setting11", "setting12"),
    model_type = "correct",
    version = "test5"
  ),

  `4` = create_scenario(
    id = "4",
    description = "Misspecified: y:u1, y:z1, z2^2 models",
    ps_spec_id = "3",
    settings = c("setting11", "setting12"),
    model_type = "misspecified",
    version = "test5"
  ),

  `7-1` = create_scenario(
    id = "7-1",
    description = "Correct: z1 subset models (full, u1_z1, u2_z1)",
    ps_spec_id = "7",
    settings = c("setting11", "setting12", "setting13", "setting14"),
    model_type = "correct",
    version = "test5"
  ),

  `7-2` = create_scenario(
    id = "7-2",
    description = "Misspecified: z1 subset models (full, u1_z1, u2_z1)",
    ps_spec_id = "7",
    settings = c("setting11", "setting12"),
    model_type = "misspecified",
    version = "test5"
  ),

  `7-3` = create_scenario(
    id = "7-3",
    description = "Misspecified: z1 subset, h_alpha full_sq1 for M1",
    ps_spec_id = "7-alt1",
    settings = c("setting11", "setting12"),
    model_type = "misspecified",
    version = "test_v1"
  ),

  `7-4` = create_scenario(
    id = "7-4",
    description = "Misspecified: z1 subset, h_alpha full_sq2 for M1",
    ps_spec_id = "7-alt2",
    settings = c("setting11", "setting12"),
    model_type = "misspecified",
    version = "test_v1"
  ),

  `8-1` = create_scenario(
    id = "8-1",
    description = "Correct: z1/z2 mixed models (full, u1_z1, u2_z2)",
    ps_spec_id = "8",
    settings = c("setting11", "setting12", "setting13", "setting14"),
    model_type = "correct",
    version = "test5"
  ),

  `8-2` = create_scenario(
    id = "8-2",
    description = "Misspecified: z1/z2 mixed models (full, u1_z1, u2_z2)",
    ps_spec_id = "8",
    settings = c("setting11", "setting12", "setting13"),
    model_type = "misspecified",
    version = "test5"
  ),

  `8-3` = create_scenario(
    id = "8-3",
    description = "Misspecified: expanded models (full+u2^2, full+z1, full+z2), h_alpha full",
    ps_spec_id = "8-alt1",
    settings = c("setting11", "setting12", "setting13"),
    model_type = "misspecified",
    version = "test5"
  ),

  `8-4` = create_scenario(
    id = "8-4",
    description = "Misspecified: expanded models, h_alpha full_sq1 for M1",
    ps_spec_id = "8-alt2",
    settings = c("setting11", "setting12", "setting13"),
    model_type = "misspecified",
    version = "test5"
  ),
  

  `9-1` = create_scenario(
    id = "9-1",
    description = "Correct: z2 subset models (full, u1_z2, u2_z2)",
    ps_spec_id = "9",
    settings = c("setting11", "setting12", "setting13", "setting14"),
    model_type = "correct",
    version = "test5"
  ),

  `9-2` = create_scenario(
    id = "9-2",
    description = "Misspecified: z2 subset models (full, u1_z2, u2_z2)",
    ps_spec_id = "9",
    settings = c("setting11", "setting12"),
    model_type = "misspecified",
    version = "test5"
  ),
  
  `9-3` = create_scenario(
    id = "9-3",
    description = "Misspecified: expanded models (full+u2^2, full+z1, full+z2), h_alpha full",
    ps_spec_id = "9-alt1",
    settings = c("setting11", "setting12", "setting13"),
    model_type = "misspecified",
    version = "test5"
  ),
  
  `9-4` = create_scenario(
    id = "9-4",
    description = "Misspecified: expanded models, h_alpha full_sq1 for M1",
    ps_spec_id = "9-alt2",
    settings = c("setting11", "setting12", "setting13"),
    model_type = "misspecified",
    version = "test5"
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
    settings = c("Cho_RM2q", "Cho_RM2p"),
    missing_rates = c("miss30", "miss50"),
    model_type = "correct",
    version = "test7"
  ),

  cho2 = create_scenario(
    id = "cho2",
    description = "Cho2025: RM2/RM3 misspecified",
    ps_spec_id = "cho2",
    settings = c("Cho_RM3q", "Cho_RM3p"),
    model_type = "correct",
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
                         missing_rates = NULL,
                         type = "HT") {

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
  if (is.null(alpha.true_config)) {
    stop(paste("No alpha.true configuration for setting:", setting,
               "in", scenario$model_type, "model list"))
  }

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
    if (is.null(alpha.true_for_rate)) {
      cat("  WARNING: No alpha.true config for", miss_rate, "\n")
      next
    }

    # Iterate over data files (each may correspond to different n values)
    for (file_idx in seq_along(data_files_for_rate)) {
      data_file <- data_files_for_rate[[file_idx]]
      alpha.true <- alpha.true_for_rate[[file_idx]]
      cat("  alpha.true [file_idx=", file_idx, "]: ",
          paste(round(alpha.true, 4), collapse=", "), "\n", sep="")

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
              h_alpha.list = ps_spec$h_alpha.list[model_set],
              inv_link = ps_spec$inv_link,
              outcome = ps_spec$outcome
            )

            # Output file
            type_suffix <- if (type == "Hajek") "_Hajek" else ""
            save_file <- paste0(
              "Simulation_Results/EBMR_IPW_", setting, "-", miss_rate,
              "-scenario", scenario_id, "_", model_str,
              "_n", current_n, "_replicate", replicate_num, "_", version, type_suffix, ".RDS"
            )

            # Check if already exists
            if (file.exists(save_file)) {
              cat("      Already exists, skipping\n")
              next
            }

            # Check for equivalent results from other scenarios (for model combos without model 1)
            # Helper function to check and copy equivalent results
            check_and_copy_equiv <- function(equiv_list) {
              equiv_scenarios <- equiv_list[[scenario_id]]
              if (!is.null(equiv_scenarios)) {
                for (equiv_scen in equiv_scenarios) {
                  equiv_file <- paste0(
                    "Simulation_Results/EBMR_IPW_", setting, "-", miss_rate,
                    "-scenario", equiv_scen, "_", model_str,
                    "_n", current_n, "_replicate", replicate_num, "_", version, type_suffix, ".RDS"
                  )
                  if (file.exists(equiv_file)) {
                    file.copy(equiv_file, save_file)
                    cat("      Copied from equivalent scenario:", equiv_scen, "\n")
                    return(TRUE)
                  }
                }
              }
              return(FALSE)
            }

            # Check equivalence based on model combination
            copied <- FALSE
            if (model_str == "1") {
              copied <- check_and_copy_equiv(EQUIVALENT_SCENARIOS_MODEL1)
            } else if (model_str == "2") {
              copied <- check_and_copy_equiv(EQUIVALENT_SCENARIOS_MODEL2)
            } else if (model_str == "3") {
              copied <- check_and_copy_equiv(EQUIVALENT_SCENARIOS_MODEL3)
            } else if (model_str == "12") {
              copied <- check_and_copy_equiv(EQUIVALENT_SCENARIOS_MODEL12)
            } else if (model_str == "13") {
              copied <- check_and_copy_equiv(EQUIVALENT_SCENARIOS_MODEL13)
            } else if (model_str == "23") {
              copied <- check_and_copy_equiv(EQUIVALENT_SCENARIOS_MODEL23)
            } else if (model_str == "123") {
              copied <- check_and_copy_equiv(EQUIVALENT_SCENARIOS_MODEL123)
            }
            if (copied) next

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
            } else if (setting %in% c("setting3", "setting4")) {
              # setting3/4: PS uses cbind(1, y, u1, u2), NOT y*u1
              function(dat, alpha.true) {
                X <- cbind(rep(1, nrow(dat)), dat$y, dat$u1, dat$u2)
                1 / (1 + exp(X %*% alpha.true))
              }
            } else {
              function(dat, alpha.true) {
                X <- cbind(rep(1, nrow(dat)), dat$y, dat$u1, dat$u2)
                1 / (1 + exp(X %*% alpha.true))
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
              save_file = save_file,
              type = type
            )

            # Print console summary after simulation completes
            # Uses clean_sim_result() from Simulation.r for consistent logic
            if (file.exists(save_file)) {
              sim_result <- readRDS(save_file)
              cleaned <- clean_sim_result(sim_result, multiplier = 3, verbose = FALSE)
              sim_clean <- cleaned$result

              # Compute mu.true for this setting
              mu.true <- get_mu_true(setting)

              # Compute Bias, ESD, ESE, CP for mu_ipw (row 1) and mu_ipw.true (row 2)
              # Row indices: mu_ipw(1), mu_ipw.true(2), se_ipw(3), se_ipw.true(4)
              bias_ipw <- round(mean(sim_clean[1, ], na.rm = TRUE) - mu.true, 3)
              esd_ipw <- round(sd(sim_clean[1, ], na.rm = TRUE), 3)
              ese_ipw <- round(mean(sim_clean[3, ], na.rm = TRUE), 3)
              ci_lower <- sim_clean[1, ] - 1.96 * sim_clean[3, ]
              ci_upper <- sim_clean[1, ] + 1.96 * sim_clean[3, ]
              cp_ipw <- round(mean((ci_lower <= mu.true) & (mu.true <= ci_upper), na.rm = TRUE), 3)

              bias_true <- round(mean(sim_clean[2, ], na.rm = TRUE) - mu.true, 3)
              esd_true <- round(sd(sim_clean[2, ], na.rm = TRUE), 3)
              ese_true <- round(mean(sim_clean[4, ], na.rm = TRUE), 3)
              ci_lower_true <- sim_clean[2, ] - 1.96 * sim_clean[4, ]
              ci_upper_true <- sim_clean[2, ] + 1.96 * sim_clean[4, ]
              cp_true <- round(mean((ci_lower_true <= mu.true) & (mu.true <= ci_upper_true), na.rm = TRUE), 3)

              cat("      Replicates: ", cleaned$n_successful, "/", cleaned$n_total,
                  " (NA:", cleaned$n_na, ", Outliers:", cleaned$n_outliers, ")\n", sep = "")
              cat("      IPW:      Bias=", bias_ipw, " ESD=", esd_ipw, " ESE=", ese_ipw, " CP=", cp_ipw, "\n", sep = "")
              cat("      IPW.true: Bias=", bias_true, " ESD=", esd_true, " ESE=", ese_true, " CP=", cp_true, "\n", sep = "")
            }
          }
        }
      }
    }
  }

  cat("\nScenario", scenario_id, "completed!\n")
}
