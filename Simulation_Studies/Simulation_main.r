# devtools::install_github("LukeChuang890212/EBMR_Methods_Reproducibility_Materials/EBMRalgorithm")
source("Simulation.r")

#------------------------------------------------------------------------------#
# Basic Setup ----
#------------------------------------------------------------------------------#
setwd("C:/Users/stat-pc/Desktop/NTHU_Research/EBMR_Methods_Reproducibility_Materials/Simulation_Studies")
data_root = "E:/Other computers/我的電腦/MNAR-Simulation/MNAR_2023/ChuangData_SM/"
replicate_num = 1000
ps_model.true = function(y, u1, u2, r, alpha.true) 1/(1+exp(cbind(rep(1, n), y, u1, u2)%*%alpha.true))
# miss.ps_model.true = function(y, u1, u2, r, n, alpha.true) 1/(1+exp(cbind(rep(1, n), y, u1, u2)%*%alpha.true))*exp(n^(-1/2)*y)

n.vector.list = list(
  correct_model = list(c(1000, 300)),
  misspecified_model = list(c(1000), c(300))
)

correct_model_all_data_file.list = list(
  setting1 = list(
    miss50 =list(
      paste0(data_root, "ChuangData_SM1.1_2000.RDS")
    ),
    miss30 =list(
      paste0(data_root, "ChuangData_SM1.2_2000.RDS")
    )
  ),
  setting2 = list(
    miss50 =list(
      paste0(data_root, "ChuangData_SM2.1_2000.RDS")
    ),
    miss30 =list(
      paste0(data_root, "ChuangData_SM2.2_2000.RDS")
    )
  )
)

misspecified_model_all_data_file.list = list(
  setting1 = list(
    miss50 =list(
      paste0(data_root, "ChuangData_SM1.1.mild2_1000.RDS"),
      paste0(data_root, "ChuangData_SM1.1.mild2_300.RDS")
    ),
    miss30 =list(
      paste0(data_root, "ChuangData_SM1.2.mild2_1000.RDS"),
      paste0(data_root, "ChuangData_SM1.2.mild2_300.RDS")
    )
  ),
  setting2 = list(
    miss50 =list(
      paste0(data_root, "ChuangData_SM2.1.mild2_1000.RDS"),
      paste0(data_root, "ChuangData_SM2.1.mild2_300.RDS")
    ),
    miss30 =list(
      paste0(data_root, "ChuangData_SM2.2.mild2_1000.RDS"),
      paste0(data_root, "ChuangData_SM2.2.mild2_300.RDS")
    )
  )
)


correct_model_alpha.true.list = list(
  setting1 = list(
    miss50 =list(
      setting1_1 = c(-0.2, 0.1, 0.1, 0.1)
    ),
    miss30 =list(
      setting1_2 = c(-1.1, 0.1, 0.1, 0.1)
    )
  ),
  setting2 = list(
    miss50 =list(
      setting2_1 = c(-0.1, 0.5, -0.5, -0.1)
    ),
    miss30 =list(
      setting2_2 = c(-0.8, 0.5, -0.5, -0.1)
    )
  )
)

misspecified_model_alpha.true.list = list(
  setting1 = list(
    miss50 =list(
      setting1_1_mild_1000 = c(0.1, -0.1, -0.1, 0.1),
      setting1_1_mild_300 = c(0.1, -0.1, -0.1, 0.1)
    ),
    miss30 =list(
      setting1_2_mild_1000 = c(-1, 0.1, 0.1, 0.1),
      setting1_2_mild_300 = c(-1, 0.1, 0.1, 0.1)
    )
  ),
  setting2 = list(
    miss50 =list(
      setting2_1_mild_1000 = c(0.05, 0.5, -0.5, -0.1),
      setting2_1_mild_300 = c(0.05, 0.5, -0.5, -0.1)
    ),
    miss30 =list(
      setting2_2_mild_1000 = c(-0.2, -0.5, -0.5, -0.1),
      setting2_2_mild_300 = c(-0.2, -0.5, -0.5, -0.1)
    )
  )
)

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
                                             version = 43)

##  Scenario 2 ----
simulate_all_settings_with_all_missing_rates(scenario = 2,
                                             full_ps_specifications,
                                             n.vector.list$misspecified_model,
                                             misspecified_model_all_data_file.list,
                                             misspecified_model_alpha.true.list,
                                             version = 43)

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
                                             version = 38)

##  Scenario 4 ----
simulate_all_settings_with_all_missing_rates(scenario = 4,
                                             full_ps_specifications,
                                             n.vector.list$misspecified_model,
                                             misspecified_model_all_data_file.list,
                                             misspecified_model_alpha.true.list,
                                             version = 38)

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
                                             version = 45)

##  Scenario 6 ----
simulate_all_settings_with_all_missing_rates(scenario = 6,
                                             full_ps_specifications,
                                             n.vector.list$misspecified_model,
                                             misspecified_model_all_data_file.list,
                                             misspecified_model_alpha.true.list,
                                             version = 45)

#------------------------------------------------------------------------------#
# Display the simulation results ----
#------------------------------------------------------------------------------#

rm.extreme = function(v){
  v = na.omit(v)
  lower = quantile(v, 0.25) - 1.5*IQR(v)
  upper = quantile(v, 0.75) + 1.5*IQR(v)
  print(paste(sum(v < lower | v > upper), length(v)))
  v = v[v >= lower & v <= upper]
  return(v)
}

summarize_results = function(sim_result, pe_index, ese_index, mu.true, is.original){
  process_sim_replicates = switch(is.original,
         `TRUE` = function(v) v,
         `FALSE` = rm.extreme)
  pe = rep(NA, length(pe_index))
  for(i in 1:length(pe_index)){
    pe[i] = mean(process_sim_replicates(sim_result[pe_index[i],]))
  }
  esd = rep(NA, length(pe_index))
  for(i in 1:length(pe_index)){
    esd[i] = sd(process_sim_replicates(sim_result[pe_index[i],]))
  }
  
  bias = pe - mu.true
  mse = round(bias^2+esd^2, 4)
  bias = round(bias, 4)
  mse = round(mse, 4)
  esd = round(esd, 4)
  
  ese = sim_result[ese_index,]
  cp = rep(NA, length(ese_index))
  for(k in 1:length(ese_index)){
    ci = cbind(sim_result[pe_index[k],]-1.96*ese, sim_result[pe_index[k],]+1.96*ese)
    coverage = apply(ci, 1,
                     function(interval) ifelse(mu.true >= interval[1] & mu.true <= interval[2], 1, 0))
    cp[k] = mean(na.omit(coverage)) 
    # se.cp[k] = sd(na.omit(coverage)) 
  }
  ese = round(mean(ese), 4)
  cp = round(cp, 4) 
  
  return(c(bias, esd, ese, mse, cp))
}

summarize_all_model_combinations_and_sample_sizes = function(setting, 
                                                             scenario, 
                                                             J,
                                                             missing_rate,
                                                             n.vector, 
                                                             replicate_num, 
                                                             mu.true,
                                                             version,
                                                             is.original){
  
  result = matrix(NA, 8*length(n.vector),  5)
  j = 1
  for(n in n.vector){
    for(model_num in 1:J){
      model_combinations = combn(J, model_num)
      for(i in 1:ncol(model_combinations)){
        # print(c(model_num, i))
        model_set = model_combinations[, i]
        read_file = paste0(c("Simulation_Results/EBMR_IPW_", setting, "-", missing_rate, "-scenario", scenario, "_", model_set, "_n", n, "_replicate", replicate_num, "_", version, ".RDS"), collapse = "")
        sim_result = readRDS(read_file)
        
        is.ipw_result = model_num == 1 & i == 1
        pe_index = ifelse(is.ipw_result, 2, 1)
        ese_index = ifelse(is.ipw_result, 4, 3)
        
        result[j, ] = summarize_results(sim_result, pe_index, ese_index, mu.true, is.original)
        j = j + 1
      }
    }
  }
  return(result)
}

summarize_all_settings_with_all_missing_rates = function(scenario,
                                                         J,
                                                         n.vector,
                                                         all_data_file.list,
                                                         version,
                                                         is.original = TRUE){
  settings = c("setting1", "setting2")
  missing_rates = c("miss50", "miss30")
  for(setting in settings){
    results_with_all_missing_rates = matrix(NA, 8*length(n.vector), 5*length(missing_rates))
    mu.true = mean(readRDS(all_data_file.list[[setting]][[1]][[1]])$y)
    for(i in 1:length((missing_rates))){
      results_with_all_missing_rates[, ((i-1)*5+1):((i-1)*5+5)] = summarize_all_model_combinations_and_sample_sizes(
        setting = setting,
        scenario = scenario,
        J = J,
        missing_rate = missing_rates[i],
        n.vector = n.vector,
        replicate_num = replicate_num,
        mu.true = mu.true,
        is.original,
        version = version
      )
    }
    estimator_names= rep(c("$\\hat{\\mu}_\\text{IPW}$",
                           "$\\hat{\\mu}_{100}$", "$\\hat{\\mu}_{010}$", "$\\hat{\\mu}_{001}$",
                           "$\\hat{\\mu}_{110}$", "$\\hat{\\mu}_{101}$", "$\\hat{\\mu}_{011}$",
                           "$\\hat{\\mu}_{111}$"), length(n.vector))
    
    results_with_all_missing_rates = cbind(estimator_names, as.data.frame(results_with_all_missing_rates)) %>% 
      as.data.frame
    colnames(res) = c("", rep(c("Bias", "ESD", "ESE", "MSE", "CP"), length(missing_rates)))
    
    print(results_with_all_missing_rates)
    # kable(results_with_all_missing_rates, align = "c", booktabs = TRUE, escape = FALSE, linesep = "") %>%
    #   kable_styling(full_width = FALSE, latex_options = c("hold_position")) %>%
    #   add_header_above(c("", "$50\\%$ missing" = 5, "$30\\%$ missing" = 5)) 
  }
}

summarize_all_settings_with_all_missing_rates(scenario = 6,
                                              J = length(full_ps_specifications),
                                              n.vector.list$misspecified_model,
                                              correct_model_all_data_file.list,
                                              version = 45)

# ChunagData_SM1.3 = readRDS("ChuangData_SM/ChuangData_SM1.3_2000.RDS")
# ChunagData_SM1.4 = readRDS("ChuangData_SM/ChuangData_SM1.4_2000.RDS")
# mu.true = c(mean(ChunagData_SM1.3$y), mean(ChunagData_SM1.4$y))
# 
# sim.dats = c("ChuangData1.3-nested", "ChuangData1.4-nested")
# 
# N.v = c(300, 1000)
# vers= rep("51", 8)
# 
# res = matrix(NA, 8*length(N.v), 5*length(sim.dats))
# for(l in 1:length(sim.dats)){
#   data = sim.dats[l]
#   result = matrix(NA, 8*length(N.v),  5)
#   m = 1
#   for(j in 1:length(N.v)){
#     N = N.v[j]
#     for(model_num in 1:3){
#       model_sets = combn(3, model_num)
#       for(i in 1:ncol(model_sets)){
#         # print(c(model_num, i))
#         model_set = model_sets[, i]
#         sim_result = source(paste0(c("ChuangResults_SM/ChuangChao2023_", data, "_", model_set, "_", N,  "_", vers[i] ,".RData"), collapse = ""))[[1]]
#         
#         if(model_num == 1 & i == 1){
#           pe_index = 2
#           ese_index = 4
#           result[m, ] = summary.sim(sim_result, pe_index, ese_index, mu.true[l])
#           m = m + 1
#         }
#         
#         pe_index = 1
#         ese_index = 3
#         result[m, ] = summary.sim(sim_result, pe_index, ese_index, mu.true[l])
#         m = m + 1
#       }
#     }
#   }
#   res[, ((l-1)*5+1):((l-1)*5+5)] = result
# }
# 
# estimator_names= rep(c("$\\hat{\\mu}_\\text{IPW}$", "$\\hat{\\mu}_{100}$", "$\\hat{\\mu}_{010}$", "$\\hat{\\mu}_{001}$",
#                  "$\\hat{\\mu}_{110}$", "$\\hat{\\mu}_{101}$", "$\\hat{\\mu}_{011}$",
#                  "$\\hat{\\mu}_{111}$"), length(N.v))
# res = cbind(estimator_names, as.data.frame(res)) %>% as.data.frame
# colnames(res) = c("", rep(c("Bias", "ESD", "ASE", "MSE", "CP"), length(sim.dats)))
# 
# kable(res, align = "c", booktabs = TRUE, escape = FALSE, linesep = "") %>%
#   kable_styling(full_width = FALSE, latex_options = c("hold_position")) %>%
#   add_header_above(c("", "$50\\%$ missing" = 5, "$30\\%$ missing" = 5)) 








