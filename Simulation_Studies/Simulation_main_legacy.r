# Install EBMRalgorithm package only if not already installed
if (!requireNamespace("EBMRalgorithm", quietly = TRUE)) {
  devtools::install_github(
    "LukeChuang890212/EBMR_Methods_Reproducibility_Materials",
    subdir = "EBMRalgorithm"
  )
}

# Set working directory to the directory containing this script
# Works in RStudio; for Rscript, set working directory before running
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}
source("Basic_setup.r")
# source("test_methods.r")
source("Data_Generation.r")
source("Simulation.r")

#------------------------------------------------------------------------------#
#  Generate data ----
#------------------------------------------------------------------------------#
set.seed(2345)
n = 2000
replicate_num = 1000

cores = detectCores()
cl <- makeCluster(cores-2)
registerDoParallel(cl)

pb = txtProgressBar(max = replicate_num, style = 3)
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
  setting12.B2 = setting12.B2,
  # Cho_RM1 = Cho_RM1,
  Cho_RM2.A1 = Cho_RM2.A1,
  Cho_RM2.A2 = Cho_RM2.A2,
  Cho_RM3.A1 = Cho_RM3.A1,
  Cho_RM3.A2 = Cho_RM3.A2,
  Cho_RM2p.A1 = Cho_RM2p.A1,
  Cho_RM2p.A2 = Cho_RM2p.A2,
  Cho_RM3p.A1 = Cho_RM3p.A1,
  Cho_RM3p.A2 = Cho_RM3p.A2,
  Cho_RM2q.A1 = Cho_RM2q.A1,
  Cho_RM2q.A2 = Cho_RM2q.A2,
  Cho_RM3q.A1 = Cho_RM3q.A1,
  Cho_RM3q.A2 = Cho_RM3q.A2
  # Cho_RM4 = Cho_RM4
)

# for(j in 1:length(setting.list)){
for(j in 1:8){
  print(names(setting.list)[j])
  start = Sys.time()
  setting = setting.list[[j]]
  dat <- foreach(i=1:1000,.combine = 'rbind', .options.snow = opts, .packages = parallel_packages)%dopar%{
    setting(n = n)
  }
  close(pb)

  print(Sys.time()-start)

  saveRDS(dat, paste0("Simulation_Data/", names(setting.list)[j], "_n", n, "_replicate", replicate_num, ".RDS"))
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

  saveRDS(dat, paste0("Simulation_Data/", names(setting.list)[j], "_n300_replicate1000.RDS"))
}

#------------------------------------------------------------------------------#
#  Setup of the simulations for the continuous and binary outcome settings ----
#  Setting 11 corresponds to the continuous outcome setting
#  setting 12 corresponds to the binary outcome setting
#------------------------------------------------------------------------------#

# settings = c("setting1", "setting2", "setting3")
settings = c("setting11", "setting12")
# settings = c("Cho_RM2", "Cho_RM3")
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
                                             version = "test8")

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
                                             version = "test8")

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
                                             version = "test8")

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
    r ~ o(y) + u2 + z1
  ),
  h_x_names.list = list(
    c("u1", "u2", "z1", "z2"),
    c("u1", "u2", "z1", "z2"),
    c("u1", "u2", "z1", "z2")
  ),
  inv_link = function(eta) 1/(1+exp(eta))
)

##  Scenario 7 ----
simulate_all_settings_with_all_missing_rates(settings = settings,
                                             missing_rates = missing_rates,
                                             scenario = "7-2",
                                             full_ps_specifications,
                                             n.vector.list$misspecified_model,
                                             misspecified_model_all_data_file.list,
                                             misspecified_model_alpha.true.list,
                                             version = "test5")

simulate_all_settings_with_all_missing_rates(settings = settings,
                                             missing_rates = missing_rates,
                                             scenario = "7-1",
                                             full_ps_specifications,
                                             n.vector.list$correct_model,
                                             correct_model_all_data_file.list,
                                             correct_model_alpha.true.list,
                                             version = "test5")

full_ps_specifications = list(
  formula.list = list(
    r ~ o(y) + u1 + u2,
    r ~ o(y) + u2 + z1,
    r ~ o(y) + u2 + z2
  ),
  h_x_names.list = list(
    c("u1", "u2", "z1", "z2"),
    c("u1", "u2", "z1", "z2"),
    c("u1", "u2", "z1", "z2")
  ),
  inv_link = function(eta) 1/(1+exp(eta))
)

##  Scenario 8 ----
simulate_all_settings_with_all_missing_rates(settings = settings,
                                             missing_rates = missing_rates,
                                             scenario = "8-2",
                                             full_ps_specifications,
                                             n.vector.list$misspecified_model,
                                             misspecified_model_all_data_file.list,
                                             misspecified_model_alpha.true.list,
                                             version = "test5")

simulate_all_settings_with_all_missing_rates(settings = settings,
                                             missing_rates = missing_rates,
                                             scenario = "8-1",
                                             full_ps_specifications,
                                             n.vector.list$correct_model,
                                             correct_model_all_data_file.list,
                                             correct_model_alpha.true.list,
                                             version = "test5")

full_ps_specifications = list(
  formula.list = list(
    r ~ o(y) + u1 + u2,
    r ~ u1 + z1 + z2,
    r ~ u2 + z1 + z2
  ),
  h_x_names.list = list(
    c("u1", "u2", "z1", "z2"),
    c("u1", "u2", "z1", "z2"),
    c("u1", "u2", "z1", "z2")
  ),
  inv_link = function(eta) 1/(1+exp(eta))
)

##  Scenario 9 ----
simulate_all_settings_with_all_missing_rates(settings = settings,
                                             missing_rates = missing_rates,
                                             scenario = "9-2",
                                             full_ps_specifications,
                                             n.vector.list$misspecified_model,
                                             misspecified_model_all_data_file.list,
                                             misspecified_model_alpha.true.list,
                                             version = "test5")

simulate_all_settings_with_all_missing_rates(settings = settings,
                                             missing_rates = missing_rates,
                                             scenario = "9-1",
                                             full_ps_specifications,
                                             n.vector.list$correct_model,
                                             correct_model_all_data_file.list,
                                             correct_model_alpha.true.list,
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
# Run the simulations for Scenario 11, 12 ----
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

#------------------------------------------------------------------------------#
# Run the simulations for Scenario 13, 14 ----
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

#------------------------------------------------------------------------------#
# Run the simulations for Scenario in Cho 2025 ----
#------------------------------------------------------------------------------#
full_ps_specifications = list(
  formula.list = list(
    r ~ x1 + x2,
    r ~ o(y) + x1,
    r ~ o(y) + x2
  ),
  h_x_names.list = list(
    c("x1", "x2"),
    c("x1", "x2"),
    c("x1", "x2")
  ),
  inv_link = function(eta) exp(eta)/(1+exp(eta))
)


##  Scenario cho1 ----
simulate_all_settings_with_all_missing_rates(settings = c("Cho_RM2q", "Cho_RM2p", "Cho_RM3q", "Cho_RM3p"),
                                             missing_rates = c("miss30", "miss50"),
                                             scenario = "cho1",
                                             full_ps_specifications,
                                             c(1000, 4000),
                                             correct_model_all_data_file.list,
                                             correct_model_alpha.true.list,
                                             version = "test7")

all_data_file.list = correct_model_all_data_file.list
alpha.true.list = correct_model_alpha.true.list

##  Scenario cho2 ----
simulate_all_settings_with_all_missing_rates(settings = c("Cho_RM2", "Cho_RM3"),
                                             missing_rates = missing_rates,
                                             scenario = "cho2",
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
