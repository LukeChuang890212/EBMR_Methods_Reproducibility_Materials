setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({ source("Data_Generation.r") })

set.seed(12345)
for (fn in c("setting1.A1", "setting1.A2", "setting2.A1", "setting2.A2")) {
  d <- get(fn)(1e7)
  cat(sprintf("%s: missing = %.4f, mean(y) = %.4f\n", fn, 1 - mean(d$r), mean(d$y)))
}
for (fn in c("setting1.B1", "setting1.B2", "setting2.B1", "setting2.B2")) {
  for (n in c(2000, 500)) {
    d <- get(fn)(n)
    cat(sprintf("%s n=%d: missing = %.4f\n", fn, n, 1 - mean(d$r)))
  }
}
