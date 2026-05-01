setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

old <- readRDS("Simulation_Results/EBMR_IPW_setting3-miss50-scenario9-3_2_n2000_replicate1000_test59.RDS")

cat("Rownames:", paste(rownames(old), collapse=", "), "\n\n")

mu_vals <- as.numeric(old["mu_ipw", ])
se_vals <- as.numeric(old["se_ipw", ])
mu_true <- as.numeric(old["mu_ipw.true", ])
se_true <- as.numeric(old["se_ipw.true", ])

cat("=== mu_ipw ===\n")
cat("True value (first):", mu_true[1], "\n")
cat("Mean:", mean(mu_vals, na.rm=TRUE), "\n")
cat("Median:", median(mu_vals, na.rm=TRUE), "\n")
cat("ESD:", sd(mu_vals, na.rm=TRUE), "\n")
cat("ASE:", mean(se_vals, na.rm=TRUE), "\n")
cat("ASE.true:", mean(se_true, na.rm=TRUE), "\n")
cat("NAs in mu:", sum(is.na(mu_vals)), "\n")
cat("NAs in se:", sum(is.na(se_vals)), "\n")
cat("Range mu:", range(mu_vals, na.rm=TRUE), "\n")
cat("Range se:", range(se_vals, na.rm=TRUE), "\n\n")

# Outlier detection
q <- quantile(mu_vals, c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99), na.rm=TRUE)
iqr <- q["75%"] - q["25%"]
lower <- q["25%"] - 3*iqr
upper <- q["75%"] + 3*iqr
outliers <- which(mu_vals < lower | mu_vals > upper)
cat("Outliers (3*IQR):", length(outliers), "out of", length(mu_vals), "\n")
if (length(outliers) > 0) {
  cat("Outlier mu values:", sort(round(mu_vals[outliers], 4)), "\n")
}

# Check alpha values for outlier reps
cat("\n=== Alpha for outlier reps (first 10) ===\n")
alpha_rows <- grep("^alpha", rownames(old))
for (idx in head(outliers, 10)) {
  alphas <- as.numeric(old[alpha_rows, idx])
  cat(sprintf("  Rep %d: mu=%.3f, ||alpha||=%.2f, alpha=(%s)\n",
      idx, mu_vals[idx], sqrt(sum(alphas^2)), paste(round(alphas, 3), collapse=", ")))
}

# ESD/ASE comparison
cat("\n=== ESD vs ASE ===\n")
cat(sprintf("All:     ESD=%.4f, ASE=%.4f, ratio=%.4f\n",
    sd(mu_vals, na.rm=TRUE), mean(se_vals, na.rm=TRUE),
    mean(se_vals, na.rm=TRUE) / sd(mu_vals, na.rm=TRUE)))

# Without outliers
keep <- !(1:length(mu_vals) %in% outliers)
cat(sprintf("No outl: ESD=%.4f, ASE=%.4f, ratio=%.4f\n",
    sd(mu_vals[keep], na.rm=TRUE), mean(se_vals[keep], na.rm=TRUE),
    mean(se_vals[keep], na.rm=TRUE) / sd(mu_vals[keep], na.rm=TRUE)))

# w.hat values
w_vals <- as.numeric(old["w.hat", ])
cat("\n=== w.hat ===\n")
cat("Unique values:", length(unique(round(w_vals, 4))), "\n")
cat("Range:", range(w_vals, na.rm=TRUE), "\n")
cat("Mean:", mean(w_vals, na.rm=TRUE), "\n")
