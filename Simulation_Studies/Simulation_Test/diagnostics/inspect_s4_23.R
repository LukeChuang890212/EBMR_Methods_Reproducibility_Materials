setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
res <- readRDS("Simulation_Results/EBMR_IPW_setting4-miss50-scenario9-3_23_n2000_replicate1000_test59.RDS")

mu_true <- 0.75
mu  <- res["mu_ipw", ]
se  <- res["se_ipw", ]
nu1 <- res["nu.hat1", ]
nu2 <- res["nu.hat2", ]

cat("=== Ensemble weight distribution ===\n")
cat(sprintf("nu1: mean=%.3f  sd=%.3f  min=%.3f  max=%.3f\n",
    mean(nu1), sd(nu1), min(nu1), max(nu1)))
cat(sprintf("nu2: mean=%.3f  sd=%.3f  min=%.3f  max=%.3f\n\n",
    mean(nu2), sd(nu2), min(nu2), max(nu2)))

# Check switching: how often is nu1 > 0.9 vs nu1 < 0.1
cat(sprintf("nu1 > 0.9 (model2 dominates): %d/%d (%.1f%%)\n",
    sum(nu1 > 0.9), length(nu1), 100*mean(nu1 > 0.9)))
cat(sprintf("nu1 < 0.1 (model3 dominates): %d/%d (%.1f%%)\n",
    sum(nu1 < 0.1), length(nu1), 100*mean(nu1 < 0.1)))
cat(sprintf("0.1 <= nu1 <= 0.9 (mixed):    %d/%d (%.1f%%)\n\n",
    sum(nu1 >= 0.1 & nu1 <= 0.9), length(nu1), 100*mean(nu1 >= 0.1 & nu1 <= 0.9)))

# Compare mu and SE by dominant model
m2_dom <- nu1 > 0.9
m3_dom <- nu1 < 0.1
mixed  <- nu1 >= 0.1 & nu1 <= 0.9

cat(sprintf("%-20s  %6s  %8s  %8s  %8s\n", "Group", "n", "mean(mu)", "ESD", "ESE"))
for (grp in list(list(m2_dom, "model2 dominates"),
                 list(m3_dom, "model3 dominates"),
                 list(mixed,  "mixed"))) {
  idx <- grp[[1]]
  lbl <- grp[[2]]
  if (sum(idx) > 1) {
    cat(sprintf("%-20s  %6d  %8.4f  %8.4f  %8.4f\n",
        lbl, sum(idx), mean(mu[idx]), sd(mu[idx]), mean(se[idx])))
  }
}

# Alpha stats by dominant model
cat("\n=== Alpha1 (intercept m2), Alpha2 (y-coef m2) by dominant model ===\n")
a1 <- res["alpha.hat1", ]
a2 <- res["alpha.hat2", ]
a5 <- res["alpha.hat5", ]
a6 <- res["alpha.hat6", ]
for (grp in list(list(m2_dom, "model2 dom"),
                 list(m3_dom, "model3 dom"),
                 list(mixed,  "mixed"))) {
  idx <- grp[[1]]
  lbl <- grp[[2]]
  if (sum(idx) > 1) {
    cat(sprintf("%-12s  a1=%7.2f  a2=%7.2f  a5=%7.2f  a6=%7.2f\n",
        lbl, mean(a1[idx]), mean(a2[idx]), mean(a5[idx]), mean(a6[idx])))
  }
}
