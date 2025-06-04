# devtools::install_github("LukeChuang890212/EBMR_Methods_Reproducibility_Materials/EBMRalgorithm")
setwd("C:/Users/stat-pc/Desktop/NTHU_Research/EBMR_Methods_Reproducibility_Materials/Simulation_Studies")
source("Basic_setup.r")
# source("test_methods.r")
source("Simulation.r")

library(dplyr)

#------------------------------------------------------------------------------#
#  Generate data ----
#------------------------------------------------------------------------------#
set.seed(1234)

cores = detectCores()
cl <- makeCluster(cores-2)
registerDoParallel(cl)

pb = txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)
parallel_packages = c()

setting.list = list(
  setting11.A1 = setting11.A1,
  setting11.A2 = setting11.A2,
  setting11.B1 = setting11.B1,
  setting11.B2 = setting11.B2,
  setting12.A1 = setting12.A1,
  setting12.A2 = setting12.A2,
  setting12.B1 = setting12.B1,
  setting12.B2 = setting12.B2
)

for(j in 1:length(setting.list)){
  print(names(setting.list)[j])
  start = Sys.time()
  setting = setting.list[[j]]
  dat <- foreach(i=1:n.sim,.combine = 'rbind', .options.snow = opts, .packages = parallel_packages)%dopar%{
    setting(n = 1000)
  }
  close(pb)

  print(Sys.time()-start)

  saveRDS(dat, paste0("Simulation_Data/", names(setting.list)[j], "_n1000_replicate1000.RDS"))
}


setting.list = list(
  setting11.B1 = setting11.B1,
  setting11.B2 = setting11.B2,
  setting12.B1 = setting12.B1,
  setting12.B2 = setting12.B2
)

for(j in 1:length(setting.list)){
  print(names(setting.list)[j])
  start = Sys.time()
  setting = setting.list[[j]]
  dat <- foreach(i=1:n.sim,.combine = 'rbind', .options.snow = opts, .packages = parallel_packages)%dopar%{
    setting(n = 300)
  }
  close(pb)

  print(Sys.time()-start)

  saveRDS(dat, paste0("ChuangData_SM/", names(setting.list)[j], "_n300_replicate1000.RDS"))
}

#------------------------------------------------------------------------------#
#  Setup of the simulations for the continuous and binary outcome settings ----
#  Setting 11 corresponds to the continuous outcome setting
#  setting 12 corresponds to the binary outcome setting
#------------------------------------------------------------------------------#

# settings = c("setting1", "setting2", "setting3")
settings = c("setting11", "setting12")
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

##  Scenario 1-1 ----
simulate_all_settings_with_all_missing_rates(settings = settings,
                                             missing_rates = missing_rates,
                                             scenario = "1-1",
                                             full_ps_specifications,
                                             n.vector.list$correct_model,
                                             correct_model_all_data_file.list,
                                             correct_model_alpha.true.list,
                                             version = "test6")

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

##  Scenario 1-2 ----
simulate_all_settings_with_all_missing_rates(settings = settings,
                                             missing_rates = missing_rates,
                                             scenario = "1-2",
                                             full_ps_specifications,
                                             n.vector.list$correct_model,
                                             correct_model_all_data_file.list,
                                             correct_model_alpha.true.list,
                                             version = "test6")

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

##  Scenario 1-3 ----
simulate_all_settings_with_all_missing_rates(settings = settings,
                                             missing_rates = missing_rates,
                                             scenario = "1-3",
                                             full_ps_specifications,
                                             n.vector.list$correct_model,
                                             correct_model_all_data_file.list,
                                             correct_model_alpha.true.list,
                                             version = "test6")

full_ps_specifications = list(
  formula.list = list(
    r ~ o(y) + u1 + u2,
    r ~ o(y) + u1,
    r ~ o(y) + u2,
    r ~ o(y) + z1,
    r ~ o(y) + z2
  ),
  h_x_names.list = list(
    c("u1", "u2", "z1", "z2"),
    c("u1", "u2", "z1", "z2"),
    c("u1", "u2", "z1", "z2"),
    c("u1", "u2", "z1", "z2"),
    c("u1", "u2", "z1", "z2")
  ),
  inv_link = function(eta) 1/(1+exp(eta))
)

##  Scenario 1-4 ----
simulate_all_settings_with_all_missing_rates(settings = settings,
                                             missing_rates = missing_rates,
                                             scenario = "1-4",
                                             full_ps_specifications,
                                             n.vector.list$correct_model,
                                             correct_model_all_data_file.list,
                                             correct_model_alpha.true.list,
                                             version = "test6")

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

##  Scenario 2-1 ----
simulate_all_settings_with_all_missing_rates(settings = settings,
                                             missing_rates = missing_rates,
                                             scenario = "2-1",
                                             full_ps_specifications,
                                             n.vector.list$misspecified_model,
                                             misspecified_model_all_data_file.list,
                                             misspecified_model_alpha.true.list,
                                             version = "test7")

full_ps_specifications = list(
  formula.list = list(
    r ~ o(y) + u1 + u2,
    r ~ o(y) + u1,
    r ~ o(y)
  ),
  h_x_names.list = list(
    c("u1", "u2", "z1"),
    c("u1", "z1", "z2"),
    c("z1", "z2")
  ),
  inv_link = function(eta) 1/(1+exp(eta))
)

##  Scenario 2-2 ----
simulate_all_settings_with_all_missing_rates(settings = settings,
                                             missing_rates = missing_rates,
                                             scenario = "2-2",
                                             full_ps_specifications,
                                             n.vector.list$misspecified_model,
                                             misspecified_model_all_data_file.list,
                                             misspecified_model_alpha.true.list,
                                             version = "test7")

full_ps_specifications = list(
  formula.list = list(
    r ~ o(y) + u1 + u2,
    r ~ o(y) + u1,
    r ~ o(y)
  ),
  h_x_names.list = list(
    c("u1", "u2", "z2"),
    c("u1", "z1", "z2"),
    c("z1", "z2")
  ),
  inv_link = function(eta) 1/(1+exp(eta))
)

##  Scenario 2-3 ----
simulate_all_settings_with_all_missing_rates(settings = settings,
                                             missing_rates = missing_rates,
                                             scenario = "2-3",
                                             full_ps_specifications,
                                             n.vector.list$misspecified_model,
                                             misspecified_model_all_data_file.list,
                                             misspecified_model_alpha.true.list,
                                             version = "test7")

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
                                             version = "test5")

##  Scenario 4 ----
simulate_all_settings_with_all_missing_rates(settings = settings,
                                             missing_rates = missing_rates,
                                             scenario = 4,
                                             full_ps_specifications,
                                             n.vector.list$misspecified_model,
                                             misspecified_model_all_data_file.list,
                                             misspecified_model_alpha.true.list,
                                             version = "test5")

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
                                             version = "test5")

##  Scenario 6 ----
simulate_all_settings_with_all_missing_rates(settings = settings,
                                             missing_rates = missing_rates,
                                             scenario = 6,
                                             full_ps_specifications,
                                             n.vector.list$misspecified_model,
                                             misspecified_model_all_data_file.list,
                                             misspecified_model_alpha.true.list,
                                             version = "test5")

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
                                             version = "test5")

##  Scenario 8 ----
simulate_all_settings_with_all_missing_rates(settings = settings,
                                             missing_rates = missing_rates,
                                             scenario = 8,
                                             full_ps_specifications,
                                             n.vector.list$misspecified_model,
                                             misspecified_model_all_data_file.list,
                                             misspecified_model_alpha.true.list,
                                             version = "test5")


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
    c("u1", "u2", "z1", "z2", "v1"),
    c("u1", "u2", "z1", "z2", "v2")
  ),
  inv_link = function(eta) 1/(1+exp(eta))
)

##  Scenario 9 ----
simulate_all_settings_with_all_missing_rates(settings = c("setting5", "setting6"),
                                             missing_rates = missing_rates,
                                             scenario = 9,
                                             full_ps_specifications,
                                             n.vector.list$correct_model,
                                             correct_model_all_data_file.list,
                                             correct_model_alpha.true.list,
                                             version = "test5")

##  Scenario 10 ----
simulate_all_settings_with_all_missing_rates(settings = c("setting5", "setting6"),
                                             missing_rates = missing_rates,
                                             scenario = 10,
                                             full_ps_specifications,
                                             n.vector.list$misspecified_model,
                                             misspecified_model_all_data_file.list,
                                             misspecified_model_alpha.true.list,
                                             version = "test5")

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
    c("u1", "u2", "v3"),
    c("u1", "u2", "v4")
  ),
  inv_link = function(eta) 1/(1+exp(eta))
)

##  Scenario 11 ----
simulate_all_settings_with_all_missing_rates(settings = c("setting5", "setting6"),
                                             missing_rates = missing_rates,
                                             scenario = 11,
                                             full_ps_specifications,
                                             n.vector.list$correct_model,
                                             correct_model_all_data_file.list,
                                             correct_model_alpha.true.list,
                                             version = "test5")

##  Scenario 12 ----
simulate_all_settings_with_all_missing_rates(settings = c("setting5", "setting6"),
                                             missing_rates = missing_rates,
                                             scenario = 12,
                                             full_ps_specifications,
                                             n.vector.list$misspecified_model,
                                             misspecified_model_all_data_file.list,
                                             misspecified_model_alpha.true.list,
                                             version = "test5")



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

ebmr = EBMRAlgorithm$new("y", full_ps_specifications, setting1.A1(1000), W)
result = ebmr$EBMR_IPW(h_x_names = c("u1", "u2"))
par(mfrow = c(1, 3))
plot(result$ps.matrix[, c(1, 2)], xlab = expression(pi[1]), ylab = expression(pi[2]))
plot(result$ps.matrix[, c(1, 3)], xlab = expression(pi[1]), ylab = expression(pi[3]))
plot(result$ps.matrix[, c(2, 3)], xlab = expression(pi[2]), ylab = expression(pi[3]))

ebmr = EBMRAlgorithm$new("y", full_ps_specifications, setting2.A1(1000), W)
result = ebmr$EBMR_IPW(h_x_names = c("u1", "u2"))
# plot(result$ps.matrix)
par(mfrow = c(1, 3))
plot(result$ps.matrix[, c(1, 2)], xlab = expression(pi[1]), ylab = expression(pi[2]))
plot(result$ps.matrix[, c(1, 3)], xlab = expression(pi[1]), ylab = expression(pi[3]))
plot(result$ps.matrix[, c(2, 3)], xlab = expression(pi[2]), ylab = expression(pi[3]))
