setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Data_Generation.r")
source("config/scenarios.R")
library(EPS)
source("Basic_setup.r")

# Find result files with good ESE/ESD ratios
result_files <- list.files("Simulation_Results", pattern="EBMR_IPW.*test5[0-9]", full.names=TRUE)

cat("=== ESE/ESD ratios across scenarios ===\n")
cat(sprintf("%-70s %-6s %-8s %-8s %-6s %-6s\n", "File", "Valid", "ESD", "ESE", "Ratio", "Outl%"))

for (f in sort(result_files)) {
  res <- tryCatch(readRDS(f), error=function(e) NULL)
  if (is.null(res)) next
  mu <- res[1,]; se <- res[3,]
  valid <- !is.na(mu) & !is.na(se)
  if (sum(valid) < 100) next

  mu_v <- mu[valid]; se_v <- se[valid]
  esd <- sd(mu_v); ese <- mean(se_v)

  q1 <- quantile(mu_v, 0.25); q3 <- quantile(mu_v, 0.75); iqr <- q3 - q1
  q1s <- quantile(se_v, 0.25); q3s <- quantile(se_v, 0.75); iqrs <- q3s - q1s
  is_out <- (mu_v < q1-3*iqr | mu_v > q3+3*iqr) | (se_v < q1s-3*iqrs | se_v > q3s+3*iqrs)

  cat(sprintf("%-70s %-6d %-8.4f %-8.4f %-6.2f %-6.1f\n",
              basename(f), sum(valid), esd, ese, ese/esd, 100*mean(is_out)))
}

cat("\nDone!\n")
