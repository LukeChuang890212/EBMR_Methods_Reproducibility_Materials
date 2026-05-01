setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({ source("Basic_setup.r"); source("Data_Generation.r"); source("config/scenarios.R") })
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
dat <- readRDS(data_file)[1:2000, ]
for (v in c("u1","u2","z1","z2")) {
  cat(sprintf("  %s: %d unique values, range [%.2f, %.2f]\n",
      v, length(unique(dat[[v]])), min(dat[[v]]), max(dat[[v]])))
}
