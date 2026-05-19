setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({source("Data_Generation.r")})

set.seed(12345)
d <- setting4.A1(1e7)
cat(sprintf("setting4.A1: mu_true = %.6f, missing rate = %.4f\n",
            mean(d$y), 1 - mean(d$r)))
d <- setting4.A2(1e7)
cat(sprintf("setting4.A2: mu_true = %.6f, missing rate = %.4f\n",
            mean(d$y), 1 - mean(d$r)))
d <- setting4.B1(2000)
cat(sprintf("setting4.B1 n=2000: missing rate = %.4f\n", 1 - mean(d$r)))
d <- setting4.B1(500)
cat(sprintf("setting4.B1 n=500:  missing rate = %.4f\n", 1 - mean(d$r)))
d <- setting4.B2(2000)
cat(sprintf("setting4.B2 n=2000: missing rate = %.4f\n", 1 - mean(d$r)))
d <- setting4.B2(500)
cat(sprintf("setting4.B2 n=500:  missing rate = %.4f\n", 1 - mean(d$r)))
