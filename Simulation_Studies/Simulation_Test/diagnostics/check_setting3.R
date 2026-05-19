setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({source("Data_Generation.r")})

set.seed(12345)
for (fn in c("setting3.A1", "setting3.A2")) {
  d <- get(fn)(1e7)
  cat(sprintf("%s: missing rate = %.4f\n", fn, 1 - mean(d$r)))
}
for (fn in c("setting3.B1", "setting3.B2")) {
  for (n in c(2000, 500)) {
    d <- get(fn)(n)
    cat(sprintf("%s n=%d: missing rate = %.4f\n", fn, n, 1 - mean(d$r)))
  }
}
