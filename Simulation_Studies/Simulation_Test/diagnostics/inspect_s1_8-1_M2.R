setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({ source("Basic_setup.r"); source("Data_Generation.r"); source("Simulation.r") })
mu_true <- get_mu_true("setting1")
for (miss in c("miss50", "miss30")) {
  f <- sprintf("Simulation_Results/EBMR_IPW_setting1-%s-scenario8-1_2_n2000_replicate1000_test65.RDS", miss)
  cat(sprintf("\n=== %s ===\n", basename(f)))
  x <- readRDS(f)
  cleaned <- clean_sim_result(x, multiplier=2, max_pct=0.01, verbose=FALSE)
  sc <- cleaned$result
  mu_v <- sc["mu_ipw",]; se_v <- sc["se_ipw",]
  bias <- mean(mu_v)-mu_true; esd <- sd(mu_v); ese <- mean(se_v); ese_med <- median(se_v)
  cat(sprintf("n=%d Out=%d  Bias=%+.4f ESD=%.4f ESE(mean)=%.4f ESE(median)=%.4f ratio=%.3f\n",
              cleaned$n_successful, cleaned$n_outliers, bias, esd, ese, ese_med, ese/esd))
  cat("top5 se_ipw:", round(sort(x["se_ipw",], decreasing=TRUE)[1:5], 3), "\n")
  a <- x[grep("^alpha.hat", rownames(x)),,drop=FALSE]
  for (i in 1:nrow(a)) cat(rownames(a)[i], ": sd=", round(sd(a[i,], na.rm=T),3), " max|a|=", round(max(abs(a[i,]), na.rm=T),1), "\n", sep="")
}
