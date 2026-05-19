setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

f <- "Simulation_Results/EBMR_IPW_setting3-miss50-scenario9-3_2_n500_replicate1000_test57.RDS"
res <- readRDS(f)

cat("Dimensions:", dim(res), "\n")
cat("Rownames:", paste(rownames(res), collapse=", "), "\n\n")

# NA counts per row
cat("=== NA counts per row ===\n")
for (i in 1:nrow(res)) {
  cat(sprintf("  %s: %d / %d NAs\n", rownames(res)[i], sum(is.na(res[i,])), ncol(res)))
}

# Which replicates have NAs?
na_cols <- which(apply(res[1:4, , drop=FALSE], 2, function(x) any(is.na(x))))
cat(sprintf("\nReplicates with NA in rows 1-4: %d / %d\n", length(na_cols), ncol(res)))

if (length(na_cols) > 0) {
  # Check pattern: which rows are NA?
  cat("\n=== NA pattern in first 20 NA reps ===\n")
  check <- na_cols[1:min(20, length(na_cols))]
  for (col_i in check) {
    na_rows <- which(is.na(res[, col_i]))
    cat(sprintf("  Rep %d: NA in rows [%s]\n", col_i, paste(rownames(res)[na_rows], collapse=", ")))
  }

  # Are ALL rows NA or just some?
  all_na <- apply(res[, na_cols, drop=FALSE], 2, function(x) all(is.na(x)))
  partial_na <- !all_na
  cat(sprintf("\n  All rows NA: %d, Partial NA: %d\n", sum(all_na), sum(partial_na)))

  # For partial NAs, which rows are present?
  if (sum(partial_na) > 0) {
    cat("\n=== Partial NA reps (first 10) ===\n")
    partial_cols <- na_cols[partial_na]
    for (col_i in partial_cols[1:min(10, length(partial_cols))]) {
      present <- which(!is.na(res[, col_i]))
      cat(sprintf("  Rep %d: present rows [%s], mu_ipw=%s, se_ipw=%s\n",
                  col_i,
                  paste(rownames(res)[present], collapse=", "),
                  ifelse(is.na(res[1, col_i]), "NA", sprintf("%.4f", res[1, col_i])),
                  ifelse(is.na(res[3, col_i]), "NA", sprintf("%.4f", res[3, col_i]))))
    }
  }
}

cat("\nDone!\n")
