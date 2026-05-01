setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
source("Simulation.r")
for (ver in c("test48", "test49", "test53", "test57", "test59")) {
  f <- paste0("Simulation_Results/EBMR_IPW_setting3-miss30-scenario9-2_2_n500_replicate1000_", ver, ".RDS")
  if (file.exists(f)) {
    sim <- readRDS(f)
    cleaned <- clean_sim_result(sim, multiplier = 3, verbose = FALSE)
    cat(sprintf("%s: total=%d, NA=%d, outliers=%d, used=%d\n",
        ver, cleaned$n_total, cleaned$n_na, cleaned$n_outliers, cleaned$n_successful))
  }
}
