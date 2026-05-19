setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
files <- c(
  "Simulation_Results/EBMR_IPW_setting2-miss50-scenario9-3_3_n2000_replicate1000_test61.RDS",
  "Simulation_Results/EBMR_IPW_setting2-miss50-scenario9-3_3_n500_replicate1000_test61.RDS",
  "Simulation_Results/EBMR_IPW_setting1-miss50-scenario9-2_3_n2000_replicate1000_test61.RDS"
)
for (f in files) {
  cat("\n===== ", basename(f), " =====\n")
  x <- readRDS(f)
  cat("dim:", paste(dim(x), collapse=" x "), "\n")
  cat("rownames:\n"); print(rownames(x))
  cat("\nFirst rep (column 1):\n")
  print(x[, 1])
  cat("\nLast rep (column 1000):\n")
  print(x[, ncol(x)])
}
