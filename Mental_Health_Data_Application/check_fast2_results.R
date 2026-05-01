setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Mental_Health_Data_Application")
labels <- c("100", "010", "001", "110", "101", "011", "111")
cat(sprintf("%-8s %10s %10s %10s\n", "Label", "mu_ipw", "se_ipw", "w.hat"))
cat(paste(rep("-", 50), collapse=""), "\n")
for (label in labels) {
  f <- sprintf("MHD_results/popmean_EBMR_IPW_%s.RDS", label)
  if (file.exists(f)) {
    res <- readRDS(f)
    w_str <- paste(round(res$w.hat, 4), collapse=", ")
    cat(sprintf("%-8s %10.4f %10.4f     (%s)\n", label, res$mu_ipw, res$se_ipw, w_str))
  }
}
