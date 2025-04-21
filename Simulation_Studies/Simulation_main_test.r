# devtools::install_github("LukeChuang890212/EBMR_Methods_Reproducibility_Materials/EBMRalgorithm")
setwd("C:/Users/stat-pc/Desktop/NTHU_Research/EBMR_Methods_Reproducibility_Materials/Simulation_Studies")
source("Basic_setup.r")
# source("test_methods.r")
source("Simulation_test.r")

library(dplyr)

settings = c("setting1", "setting2", "setting3")
missing_rates = c("miss50", "miss30")
replicate_num = 1000

W = function(g.matrix){
  return(solve(t(g.matrix)%*%g.matrix/n))
}

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
simulate_all_settings_with_all_missing_rates(settings = settings,
                                             missing_rates = missing_rates,
                                             scenario = 1,
                                             full_ps_specifications,
                                             n.vector.list$correct_model,
                                             correct_model_all_data_file.list,
                                             correct_model_alpha.true.list,
                                             version = "test4")

##  Scenario 2 ----
simulate_all_settings_with_all_missing_rates(settings = settings,
                                             missing_rates = missing_rates,
                                             scenario = 2,
                                             full_ps_specifications,
                                             n.vector.list$misspecified_model,
                                             misspecified_model_all_data_file.list,
                                             misspecified_model_alpha.true.list,
                                             version = "test4")

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
simulate_all_settings_with_all_missing_rates(settings = settings,
                                             missing_rates = missing_rates,
                                             scenario = 3,
                                             full_ps_specifications,
                                             n.vector.list$correct_model,
                                             correct_model_all_data_file.list,
                                             correct_model_alpha.true.list,
                                             version = "test4")

##  Scenario 4 ----
simulate_all_settings_with_all_missing_rates(settings = settings,
                                             missing_rates = missing_rates,
                                             scenario = 4,
                                             full_ps_specifications,
                                             n.vector.list$misspecified_model,
                                             misspecified_model_all_data_file.list,
                                             misspecified_model_alpha.true.list,
                                             version = "test4")

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
simulate_all_settings_with_all_missing_rates(settings = settings,
                                             missing_rates = missing_rates,
                                             scenario = 5,
                                             full_ps_specifications,
                                             n.vector.list$correct_model,
                                             correct_model_all_data_file.list,
                                             correct_model_alpha.true.list,
                                             version = "test4")

##  Scenario 6 ----
simulate_all_settings_with_all_missing_rates(settings = settings,
                                             missing_rates = missing_rates,
                                             scenario = 6,
                                             full_ps_specifications,
                                             n.vector.list$misspecified_model,
                                             misspecified_model_all_data_file.list,
                                             misspecified_model_alpha.true.list,
                                             version = "test4")

#------------------------------------------------------------------------------#
# Run the simulations for Scenario 7, 8 ----
#------------------------------------------------------------------------------#
full_ps_specifications = list(
  formula.list = list(
    r ~ o(y) + u1 + u2,
    r ~ o(y) + u1 + z1,
    r ~ o(y) + u2 + z2
  ),
  h_x_names.list = list(
    c("u1", "u2", "z1", "z2"),
    c("u1", "z1", "z2"),
    c("u2", "z1", "z2")
  ),
  inv_link = function(eta) 1/(1+exp(eta))
)

##  Scenario 7 ----
simulate_all_settings_with_all_missing_rates(settings = settings,
                                             missing_rates = missing_rates,
                                             scenario = 7,
                                             full_ps_specifications,
                                             n.vector.list$correct_model,
                                             correct_model_all_data_file.list,
                                             correct_model_alpha.true.list,
                                             version = "test4")

##  Scenario 8 ----
simulate_all_settings_with_all_missing_rates(settings = settings,
                                             missing_rates = missing_rates,
                                             scenario = 8,
                                             full_ps_specifications,
                                             n.vector.list$misspecified_model,
                                             misspecified_model_all_data_file.list,
                                             misspecified_model_alpha.true.list,
                                             version = "test4")

#------------------------------------------------------------------------------#
# Run the simulations for Scenario 9, 10 ----
#------------------------------------------------------------------------------#
full_ps_specifications = list(
  formula.list = list(
    r ~ o(y) + u1 + u2,
    r ~ o(y) + u1 + u2,
    r ~ o(y) + u1 + u2
  ),
  h_x_names.list = list(
    c("u1", "u2", "z1", "z2"),
    c("u1", "u2", "z1"),
    c("u1", "u2", "z2")
  ),
  inv_link = function(eta) 1/(1+exp(eta))
)

##  Scenario 9 ----
simulate_all_settings_with_all_missing_rates(settings = settings,
                                             missing_rates = missing_rates,
                                             scenario = 9,
                                             full_ps_specifications,
                                             n.vector.list$correct_model,
                                             correct_model_all_data_file.list,
                                             correct_model_alpha.true.list,
                                             version = "test4")

##  Scenario 10 ----
simulate_all_settings_with_all_missing_rates(settings = settings,
                                             missing_rates = missing_rates,
                                             scenario = 10,
                                             full_ps_specifications,
                                             n.vector.list$misspecified_model,
                                             misspecified_model_all_data_file.list,
                                             misspecified_model_alpha.true.list,
                                             version = "test4")





