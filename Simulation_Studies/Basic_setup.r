#------------------------------------------------------------------------------#
# Basic Setup ----
#------------------------------------------------------------------------------#
# data_root = "E:/Other computers/我的電腦/MNAR-Simulation/MNAR_2023/ChuangData_SM/"
# data_root = "/Users/luke/Library/CloudStorage/GoogleDrive-mingjuresearch@gmail.com/其他電腦/我的電腦/MNAR-Simulation/MNAR_2023/ChuangData_SM/"
data_root = "G:/Other computers/我的電腦/MNAR-Simulation/MNAR_2023/ChuangData_SM/"

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
  ),
  setting3 = list(
    miss50 =list(
      paste0(data_root, "Setting3.A1_n1000_replicate1000.RDS")
    ),
    miss30 =list(
      paste0(data_root, "Setting3.A2_n1000_replicate1000.RDS")
    )
  ),
  setting4 = list(
    miss50 =list(
      paste0(data_root, "Setting4.A1_n1000_replicate1000.RDS")
    ),
    miss30 =list(
      paste0(data_root, "Setting4.A2_n1000_replicate1000.RDS")
    )
  ),
  setting5 = list(
    miss50 =list(
      paste0(data_root, "Setting5.A1_n1000_replicate1000.RDS")
    ),
    miss30 =list(
      paste0(data_root, "Setting5.A2_n1000_replicate1000.RDS")
    )
  ),
  setting6 = list(
    miss50 =list(
      paste0(data_root, "Setting6.A1_n1000_replicate1000.RDS")
    ),
    miss30 =list(
      paste0(data_root, "Setting6.A2_n1000_replicate1000.RDS")
    )
  ),
  setting7 = list(
    miss50 =list(
      paste0(data_root, "setting7.A1_n1000_replicate1000.RDS")
    ),
    miss30 =list(
      paste0(data_root, "setting7.A2_n1000_replicate1000.RDS")
    )
  ),
  setting8 = list(
    miss50 =list(
      paste0(data_root, "setting8.A1_n1000_replicate1000.RDS")
    ),
    miss30 =list(
      paste0(data_root, "setting8.A2_n1000_replicate1000.RDS")
    )
  ),
  setting9 = list(
    miss50 =list(
      paste0(data_root, "setting9.A1_n1000_replicate1000.RDS")
    ),
    miss30 =list(
      paste0(data_root, "setting9.A2_n1000_replicate1000.RDS")
    )
  ),
  setting10 = list(
    miss50 =list(
      paste0(data_root, "setting10.A1_n1000_replicate1000.RDS")
    ),
    miss30 =list(
      paste0(data_root, "setting10.A2_n1000_replicate1000.RDS")
    )
  ),
  setting11 = list(
    miss50 =list(
      paste0(data_root, "setting11.A1_n1000_replicate1000.RDS")
    ),
    miss30 =list(
      paste0(data_root, "setting11.A2_n1000_replicate1000.RDS")
    )
  ),
  setting12 = list(
    miss50 =list(
      paste0(data_root, "setting12.A1_n1000_replicate1000.RDS")
    ),
    miss30 =list(
      paste0(data_root, "setting12.A2_n1000_replicate1000.RDS")
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
  ),
  setting3 = list(
    miss50 =list(
      paste0(data_root, "Setting3.B1_n1000_replicate1000.RDS"),
      paste0(data_root, "Setting3.B1_n300_replicate1000.RDS")
    ),
    miss30 =list(
      paste0(data_root, "Setting3.B2_n1000_replicate1000.RDS"),
      paste0(data_root, "Setting3.B2_n300_replicate1000.RDS")
    )
  ),
  setting4 = list(
    miss50 =list(
      paste0(data_root, "Setting4.B1_n1000_replicate1000.RDS"),
      paste0(data_root, "Setting4.B1_n300_replicate1000.RDS")
    ),
    miss30 =list(
      paste0(data_root, "Setting4.B2_n1000_replicate1000.RDS"),
      paste0(data_root, "Setting4.B2_n300_replicate1000.RDS")
    )
  ),
  setting5 = list(
    miss50 =list(
      paste0(data_root, "Setting5.B1_n1000_replicate1000.RDS"),
      paste0(data_root, "Setting5.B1_n300_replicate1000.RDS")
    ),
    miss30 =list(
      paste0(data_root, "Setting5.B2_n1000_replicate1000.RDS"),
      paste0(data_root, "Setting5.B2_n300_replicate1000.RDS")
    )
  ),
  setting6 = list(
    miss50 =list(
      paste0(data_root, "Setting6.B1_n1000_replicate1000.RDS"),
      paste0(data_root, "Setting6.B1_n300_replicate1000.RDS")
    ),
    miss30 =list(
      paste0(data_root, "Setting6.B2_n1000_replicate1000.RDS"),
      paste0(data_root, "Setting6.B2_n300_replicate1000.RDS")
    )
  ),
  setting7 = list(
    miss50 =list(
      paste0(data_root, "setting7.B1_n1000_replicate1000.RDS"),
      paste0(data_root, "setting7.B1_n300_replicate1000.RDS")
    ),
    miss30 =list(
      paste0(data_root, "setting7.B2_n1000_replicate1000.RDS"),
      paste0(data_root, "setting7.B2_n300_replicate1000.RDS")
    )
  ),
  setting8 = list(
    miss50 =list(
      paste0(data_root, "setting8.B1_n1000_replicate1000.RDS"),
      paste0(data_root, "setting8.B1_n300_replicate1000.RDS")
    ),
    miss30 =list(
      paste0(data_root, "setting8.B2_n1000_replicate1000.RDS"),
      paste0(data_root, "Setting8.B2_n300_replicate1000.RDS")
    )
  ),
  setting9 = list(
    miss50 =list(
      paste0(data_root, "setting9.B1_n1000_replicate1000.RDS"),
      paste0(data_root, "setting9.B1_n300_replicate1000.RDS")
    ),
    miss30 =list(
      paste0(data_root, "setting9.B2_n1000_replicate1000.RDS"),
      paste0(data_root, "setting9.B2_n300_replicate1000.RDS")
    )
  ),
  setting10 = list(
    miss50 =list(
      paste0(data_root, "setting10.B1_n1000_replicate1000.RDS"),
      paste0(data_root, "setting10.B1_n300_replicate1000.RDS")
    ),
    miss30 =list(
      paste0(data_root, "setting10.B2_n1000_replicate1000.RDS"),
      paste0(data_root, "Setting10.B2_n300_replicate1000.RDS")
    )
  ),
  setting11 = list(
    miss50 =list(
      paste0(data_root, "setting11.B1_n1000_replicate1000.RDS"),
      paste0(data_root, "setting11.B1_n300_replicate1000.RDS")
    ),
    miss30 =list(
      paste0(data_root, "setting11.B2_n1000_replicate1000.RDS"),
      paste0(data_root, "setting11.B2_n300_replicate1000.RDS")
    )
  ),
  setting12 = list(
    miss50 =list(
      paste0(data_root, "setting12.B1_n1000_replicate1000.RDS"),
      paste0(data_root, "setting12.B1_n300_replicate1000.RDS")
    ),
    miss30 =list(
      paste0(data_root, "setting12.B2_n1000_replicate1000.RDS"),
      paste0(data_root, "Setting12.B2_n300_replicate1000.RDS")
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
  ),
  setting3 = list(
    miss50 =list(
      setting3_1 = c(-0.2, 0.1, 0.1, 0.1)
    ),
    miss30 =list(
      setting3_2 = c(-1.1, 0.1, 0.1, 0.1)
    )
  ),
  setting5 = list(
    miss50 =list(
      setting3_1 = c(-0.2, 0.1, 0.1, 0.1)
    ),
    miss30 =list(
      setting3_2 = c(-1.1, 0.1, 0.1, 0.1)
    )
  ),
  setting6 = list(
    miss50 =list(
      setting2_1 = c(-0.1, 0.5, -0.5, -0.1)
    ),
    miss30 =list(
      setting2_2 = c(-0.8, 0.5, -0.5, -0.1)
    )
  ),
  setting7 = list(
    miss50 =list(
      setting7_1 = c(0.2, 0.1, -0.5, -0.5)
    ),
    miss30 =list(
      setting7_2 = c(-1, 0.1, -0.5, -0.5)
    )
  ),
  setting8 = list(
    miss50 =list(
      setting8_1 = c(-0.04, 0.2, -0.2, -0.2)
    ),
    miss30 =list(
      setting8_2 = c(-0.9, 0.2, -0.2, -0.2)
    )
  ),
  setting9 = list(
    miss50 =list(
      setting9_1 = c(0.6, 0.2, -1, -0.5)
    ),
    miss30 =list(
      setting9_2 = c(-0.5, 0.2, -1, -0.5)
    )
  ),
  setting10 = list(
    miss50 =list(
      setting10_1 = c(0.4, 0.4, -1, -0.5)
    ),
    miss30 =list(
      setting10_2 = c(-0.6, 0.4, -1, -0.5)
    )
  ),
  setting11 = list(
    miss50 =list(
      setting9_1 = c(1.3, 0.2, -2, -0.5)
    ),
    miss30 =list(
      setting9_2 = c(0.2, 0.2, -2, -0.5)
    )
  ),
  setting12 = list(
    miss50 =list(
      setting10_1 = c(1.1, 0.4, -2, -0.5)
    ),
    miss30 =list(
      setting10_2 = c(0.05, 0.4, -2, -0.5)
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
  ),
  setting3 = list(
    miss50 =list(
      setting3_1_mild_1000 = c(-0.15, 0.1, 0.1, 0.1),
      setting3_1_mild_300 = c(-0.15, 0.1, 0.1, 0.1)
    ),
    miss30 =list(
      setting3_2_mild_1000 = c(-0.9, 0.1, 0.1, 0.1),
      setting3_2_mild_300 = c(-0.9, 0.1, 0.1, 0.1)
    )
  ),
  setting5 = list(
    miss50 =list(
      setting5_1_mild_1000 = c(-0.15, 0.1, 0.1, 0.1),
      setting5_1_mild_300 = c(-0.15, 0.1, 0.1, 0.1)
    ),
    miss30 =list(
      setting5_2_mild_1000 = c(-0.9, 0.1, 0.1, 0.1),
      setting5_2_mild_300 = c(-0.9, 0.1, 0.1, 0.1)
    )
  ),
  setting6 = list(
    miss50 =list(
      setting6_1_mild_1000 = c(0.05, 0.5, -0.5, -0.1),
      setting6_1_mild_300 = c(0.05, 0.5, -0.5, -0.1)
    ),
    miss30 =list(
      setting6_2_mild_1000 = c(-0.2, -0.5, -0.5, -0.1),
      setting6_2_mild_300 = c(-0.2, -0.5, -0.5, -0.1)
    )
  ),
  setting7 = list(
    miss50 =list(
      setting7_1_mild_1000 = c(0.2, 0.1, -0.5, -0.5),
      setting7_1_mild_300 = c(0.2, 0.1, -0.5, -0.5)
    ),
    miss30 =list(
      setting7_2_mild_1000 = c(-1, 0.1, -0.5, -0.5),
      setting7_2_mild_300 = c(-1, 0.1, -0.5, -0.5)
    )
  ),
  setting8 = list(
    miss50 =list(
      setting8_1_mild_1000 = c(-0.04, 0.2, -0.2, -0.2),
      setting8_1_mild_300 = c(-0.04, 0.2, -0.2, -0.2)
    ),
    miss30 =list(
      setting8_2_mild_1000 = c(-0.9, 0.2, -0.2, -0.2),
      setting8_2_mild_300 = c(-0.9, 0.2, -0.2, -0.2)
    )
  ),
  setting9 = list(
    miss50 =list(
      setting7_1_mild_1000 =  c(0.6, 0.2, -1, -0.5),
      setting7_1_mild_300 =  c(0.6, 0.2, -1, -0.5)
    ),
    miss30 =list(
      setting7_2_mild_1000 = c(-0.5, 0.2, -1, -0.5),
      setting7_2_mild_300 = c(-0.5, 0.2, -1, -0.5)
    )
  ),
  setting10 = list(
    miss50 =list(
      setting8_1_mild_1000 = c(0.4, 0.4, -1, -0.5),
      setting8_1_mild_300 = c(0.4, 0.4, -1, -0.5)
    ),
    miss30 =list(
      setting8_2_mild_1000 = c(-0.6, 0.4, -1, -0.5),
      setting8_2_mild_300 = c(-0.6, 0.4, -1, -0.5)
    )
  ),
  setting11 = list(
    miss50 =list(
      setting11_1_mild_1000 =  c(1.3, 0.2, -2, -0.5),
      setting11_1_mild_300 =  c(1.3, 0.2, -2, -0.5)
    ),
    miss30 =list(
      setting11_2_mild_1000 = c(0.2, 0.2, -2, -0.5),
      setting11_2_mild_300 = c(0.2, 0.2, -2, -0.5)
    )
  ),
  setting12 = list(
    miss50 =list(
      setting12_1_mild_1000 = c(1.1, 0.4, -2, -0.5),
      setting12_1_mild_300 = c(1.1, 0.4, -2, -0.5)
    ),
    miss30 =list(
      setting12_2_mild_1000 = c(0.05, 0.4, -2, -0.5),
      setting12_2_mild_300 = c(0.05, 0.4, -2, -0.5)
    )
  )
)
