setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
x <- readRDS("Simulation_Results/EBMR_IPW_setting15-miss50-scenario9-1_12_n2000_replicate1000_test63.RDS")
cat("class:", class(x), "\n")
cat("dim:", dim(x), "\n")
cat("nrow:", nrow(x), "\n")
cat("ncol:", ncol(x), "\n")
if (!is.null(dim(x))) {
  cat("rownames:", paste(rownames(x), collapse=", "), "\n")
  cat("first few cols:\n")
  print(x[, 1:min(3, ncol(x))])
}
