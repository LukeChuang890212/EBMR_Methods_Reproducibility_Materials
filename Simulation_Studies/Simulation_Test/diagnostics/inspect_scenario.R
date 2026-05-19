setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
source("config/scenarios.R")
ps <- get_ps_spec("9-3")
cat("Models in 9-3:\n")
for (i in seq_along(ps$formula.list)) cat("  Model", i, ":", deparse(ps$formula.list[[i]]), "\n")
cat("h_alpha:\n")
for (i in seq_along(ps$h_alpha.list)) {
  h <- ps$h_alpha.list[[i]]
  cat("  Model", i, ":", if (is.character(h)) paste(h, collapse=", ") else "function", "\n")
}
cat("outcome:", ps$outcome, "\n")

# Check setting4 data files
cat("\nSetting4 data files:\n")
f <- list.files("Simulation_Data", pattern="Setting4", full.names=FALSE)
cat(paste(f, collapse="\n"), "\n")
