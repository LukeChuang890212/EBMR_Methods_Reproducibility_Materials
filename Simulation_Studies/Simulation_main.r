# devtools::install_github("LukeChuang890212/EBMR_Methods_Reproducibility_Materials/EBMRalgorithm")
source("Simulation.r")
source("Basic_setup.r")

library(dplyr)

#------------------------------------------------------------------------------#
#  Run the simulations for Scenario 1, 2 ----
#------------------------------------------------------------------------------#
full_ps_specifications = list(
  formula.list = list(
    r ~ o(y) + u1 + u2,
    r ~ o(y) + u1,
    r ~ o(y)
  ),
  h_x_names.list = list(
    c("u1", "u2", "z1", "z2"),
    c("u1", "z1", "z2"),
    c("z1", "z2")
  ),
  inv_link = function(eta) 1/(1+exp(eta))
)

##  Scenario 1 ----
simulate_all_settings_with_all_missing_rates(scenario = 1,
                                             full_ps_specifications,
                                             n.vector.list$correct_model,
                                             correct_model_all_data_file.list,
                                             correct_model_alpha.true.list,
                                             version = 2)

##  Scenario 2 ----
simulate_all_settings_with_all_missing_rates(scenario = 2,
                                             full_ps_specifications,
                                             n.vector.list$misspecified_model,
                                             misspecified_model_all_data_file.list,
                                             misspecified_model_alpha.true.list,
                                             version = 2)

#------------------------------------------------------------------------------#
# Run the simulations for Scenario 3, 4 ----
#------------------------------------------------------------------------------#
full_ps_specifications = list(
  formula.list = list(
    r ~ o(y) + u1 + u2,
    r ~ o(y) + u2,
    r ~ o(y)
  ),
  h_x_names.list = list(
    c("u1", "u2", "z1", "z2"),
    c("u2", "z1", "z2"),
    c("z1", "z2")
  ),
  inv_link = function(eta) 1/(1+exp(eta))
)

##  Scenario 3 ----
simulate_all_settings_with_all_missing_rates(scenario = 3,
                                             full_ps_specifications,
                                             n.vector.list$correct_model,
                                             correct_model_all_data_file.list,
                                             correct_model_alpha.true.list,
                                             version = 2)

##  Scenario 4 ----
simulate_all_settings_with_all_missing_rates(scenario = 4,
                                             full_ps_specifications,
                                             n.vector.list$misspecified_model,
                                             misspecified_model_all_data_file.list,
                                             misspecified_model_alpha.true.list,
                                             version = 2)

#------------------------------------------------------------------------------#
# Run the simulations for Scenario 5, 6 ----
#------------------------------------------------------------------------------#
full_ps_specifications = list(
  formula.list = list(
    r ~ o(y) + u1 + u2,
    r ~ o(y) + u1,
    r ~ o(y) + u2
  ),
  h_x_names.list = list(
    c("u1", "u2", "z1", "z2"),
    c("u1", "z1", "z2"),
    c("u2", "z1", "z2")
  ),
  inv_link = function(eta) 1/(1+exp(eta))
)

##  Scenario 5 ----
simulate_all_settings_with_all_missing_rates(scenario = 5,
                                             full_ps_specifications,
                                             n.vector.list$correct_model,
                                             correct_model_all_data_file.list,
                                             correct_model_alpha.true.list,
                                             version = 2)

##  Scenario 6 ----
simulate_all_settings_with_all_missing_rates(scenario = 6,
                                             full_ps_specifications,
                                             n.vector.list$misspecified_model,
                                             misspecified_model_all_data_file.list,
                                             misspecified_model_alpha.true.list,
                                             version = 2)





