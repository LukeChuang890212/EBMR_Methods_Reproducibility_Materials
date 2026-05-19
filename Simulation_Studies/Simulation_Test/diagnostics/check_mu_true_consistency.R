setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r"); source("Simulation.r")
})

# Direct computation
set.seed(12345)
d1 <- setting4.A1(1e7)
cat(sprintf("DIRECT setting4.A1(1e7): mean(y) = %.6f\n", mean(d1$y)))
set.seed(12345)
d2 <- setting4.A2(1e7)
cat(sprintf("DIRECT setting4.A2(1e7): mean(y) = %.6f\n", mean(d2$y)))

# What get_mu_true returns
cat(sprintf("get_mu_true('setting4') = %.6f\n", get_mu_true("setting4")))

# Show z2 line in setting4.A1
cat("\n--- setting4.A1 body (z2 line) ---\n")
b <- deparse(body(setting4.A1))
cat(grep("z2", b, value = TRUE), sep = "\n")
