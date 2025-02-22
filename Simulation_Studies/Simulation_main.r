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
#  Scenario 1, 2 ----
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
sim_all_settings_with_all_missing_mechanisms(scenario = 1,
                                             full_ps_specifications,
                                             n.vector.list$correct_model,
                                             correct_model_all_data_file.list,
                                             correct_model_alpha.true.list,
                                             version = 43)

##  Scenario 2 ----
sim_all_settings_with_all_missing_mechanisms(scenario = 2,
                                             full_ps_specifications,
                                             n.vector.list$misspecified_model,
                                             misspecified_model_all_data_file.list,
                                             misspecified_model_alpha.true.list,
                                             version = 43)

#------------------------------------------------------------------------------#
# Scenario 3, 4 ----
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
sim_all_settings_with_all_missing_mechanisms(scenario = 3,
                                             full_ps_specifications,
                                             n.vector.list$correct_model,
                                             correct_model_all_data_file.list,
                                             correct_model_alpha.true.list,
                                             version = 38)

##  Scenario 4 ----
sim_all_settings_with_all_missing_mechanisms(scenario = 4,
                                             full_ps_specifications,
                                             n.vector.list$misspecified_model,
                                             misspecified_model_all_data_file.list,
                                             misspecified_model_alpha.true.list,
                                             version = 38)

#------------------------------------------------------------------------------#
# Scenario 5, 6 ----
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
sim_all_settings_with_all_missing_mechanisms(scenario = 5,
                                             full_ps_specifications,
                                             n.vector.list$correct_model,
                                             correct_model_all_data_file.list,
                                             correct_model_alpha.true.list,
                                             version = 45)

##  Scenario 6 ----
sim_all_settings_with_all_missing_mechanisms(scenario = 6,
                                             full_ps_specifications,
                                             n.vector.list$misspecified_model,
                                             misspecified_model_all_data_file.list,
                                             misspecified_model_alpha.true.list,
                                             version = 45)








