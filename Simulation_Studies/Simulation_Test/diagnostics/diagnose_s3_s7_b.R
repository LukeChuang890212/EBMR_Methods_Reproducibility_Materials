## Diagnose: check nu/w, alpha, and identify problematic reps
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

res <- readRDS("Simulation_Results/EBMR_IPW_setting3-miss50-scenario7-1_2_n2000_replicate1000_test61.RDS")

mu_ipw <- res[1,]
se_ipw <- res[3,]
nu <- res["nu.hat",]
w <- res["w.hat",]

cat("nu distribution:\n")
print(summary(as.numeric(nu)))
cat("\nw distribution:\n")
print(summary(as.numeric(w)))

# Check if nu is stored as single value (J=1 skipped) or vector
cat("\nClass of nu[1]:", class(nu[[1]]), "\n")
cat("Length of nu[1]:", length(nu[[1]]), "\n")
if (length(nu[[1]]) > 1) {
  nu1 <- sapply(nu, function(x) x[1])
  nu2 <- sapply(nu, function(x) x[2])
  w1 <- sapply(w, function(x) x[1])
  w2 <- sapply(w, function(x) x[2])
  cat("\nnu1:", summary(nu1), "\n")
  cat("nu2:", summary(nu2), "\n")
  cat("w1:", summary(w1), "\n")
  cat("w2:", summary(w2), "\n")
  cat("\nw1 bimodal check: mean=", mean(w1), "median=", median(w1), "sd=", sd(w1), "\n")
} else {
  cat("\nnu values (first 20):", head(as.numeric(nu), 20), "\n")
  cat("w values (first 20):", head(as.numeric(w), 20), "\n")
}

# Look at problematic reps (mu > 5)
extreme <- which(mu_ipw > 5)
cat("\n\nExtreme reps (mu > 5):", length(extreme), "\n")
if (length(extreme) > 0) {
  for (i in head(extreme, 5)) {
    cat(sprintf("  Rep %d: mu=%.2f, se=%.2f, nu=%s, w=%s\n",
        i, mu_ipw[i], se_ipw[i],
        paste(round(unlist(nu[i]), 4), collapse=","),
        paste(round(unlist(w[i]), 4), collapse=",")))
    cat(sprintf("         alpha=(%s)\n",
        paste(round(c(res["alpha.hat1",i], res["alpha.hat2",i],
                      res["alpha.hat3",i], res["alpha.hat4",i]), 4), collapse=",")))
  }
}

# Alpha distribution
cat("\nalpha distributions:\n")
for (k in 1:4) {
  a <- as.numeric(res[paste0("alpha.hat", k),])
  cat(sprintf("  alpha%d: mean=%.4f sd=%.4f min=%.4f max=%.4f\n",
      k, mean(a), sd(a), min(a), max(a)))
}

# Check convergence: are alpha SE reasonable?
cat("\nalpha SE distributions:\n")
for (k in 1:4) {
  se <- as.numeric(res[paste0("se_alpha.hat", k),])
  cat(sprintf("  se_alpha%d: mean=%.4f sd=%.4f min=%.4f max=%.4f\n",
      k, mean(se, na.rm=TRUE), sd(se, na.rm=TRUE), min(se, na.rm=TRUE), max(se, na.rm=TRUE)))
}
