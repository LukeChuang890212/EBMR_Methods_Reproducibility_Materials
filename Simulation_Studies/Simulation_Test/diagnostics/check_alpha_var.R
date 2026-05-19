setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
r <- readRDS("Simulation_Test/test_constrained_nr_cond1e2_results.RDS")
for (nm in names(r)) {
  res <- r[[nm]]$res
  alpha_mat <- res[, c("a1","a2","a3","a4")]
  n_zero <- sum(rowSums(abs(alpha_mat)) == 0, na.rm = TRUE)
  cat(sprintf("\n%s\n", nm))
  cat(sprintf("  reps with alpha == (0,0,0,0): %d / %d\n", n_zero, nrow(alpha_mat)))
  cat(sprintf("  sd(a1, a2, a3, a4) across reps: %.4f, %.4f, %.4f, %.4f\n",
              sd(alpha_mat[,1], na.rm=TRUE), sd(alpha_mat[,2], na.rm=TRUE),
              sd(alpha_mat[,3], na.rm=TRUE), sd(alpha_mat[,4], na.rm=TRUE)))
}
