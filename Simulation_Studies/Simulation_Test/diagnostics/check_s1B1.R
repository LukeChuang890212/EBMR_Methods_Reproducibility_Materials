setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({ source("Data_Generation.r") })

# Average missing rate over many reps at n=2000 and n=500
set.seed(2026)
for (nn in c(2000, 500)) {
  rates <- replicate(200, 1 - mean(setting1.B1(nn)$r))
  cat(sprintf("setting1.B1 n=%d: mean missing = %.4f (sd=%.4f) over 200 reps\n",
              nn, mean(rates), sd(rates)))
}
# Big single sample (uses n>=1000 branch -> alpha0 = -0.0062)
set.seed(1)
big <- setting1.B1(1e6)
cat(sprintf("setting1.B1 n=1e6: missing = %.4f\n", 1 - mean(big$r)))
