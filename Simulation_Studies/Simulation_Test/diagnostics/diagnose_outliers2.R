setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

f <- "Simulation_Results/EBMR_IPW_setting3-miss50-scenario9-3_1_n2000_replicate1000_test57.RDS"
res <- readRDS(f)

cat("File:", f, "\n")
cat("Class:", class(res), "\n")
cat("Dim:", dim(res), "\n")
cat("Colnames:", paste(head(colnames(res), 20), collapse=", "), "\n")
cat("Rownames:", paste(head(rownames(res), 20), collapse=", "), "\n\n")

# Show first few rows/cols
cat("Head:\n")
print(res[1:min(5, nrow(res)), 1:min(10, ncol(res))])
cat("\n")

# If it's replicate x metrics matrix
if (nrow(res) > ncol(res)) {
  cat("Assuming rows = replicates, cols = metrics\n\n")
  for (k in 1:ncol(res)) {
    v <- res[, k]
    valid <- !is.na(v)
    q1 <- quantile(v[valid], 0.25); q3 <- quantile(v[valid], 0.75); iqr <- q3 - q1
    outlier3 <- sum(v[valid] < (q1 - 3*iqr) | v[valid] > (q3 + 3*iqr))
    nm <- if (!is.null(colnames(res))) colnames(res)[k] else paste0("col", k)
    cat(sprintf("%-25s: mean=%10.4f, sd=%8.4f, min=%10.4f, max=%10.4f, outliers(3IQR)=%d\n",
                nm, mean(v[valid]), sd(v[valid]), min(v[valid]), max(v[valid]), outlier3))
  }
} else {
  cat("Assuming cols = replicates, rows = metrics\n\n")
  for (k in 1:nrow(res)) {
    v <- res[k, ]
    valid <- !is.na(v)
    q1 <- quantile(v[valid], 0.25); q3 <- quantile(v[valid], 0.75); iqr <- q3 - q1
    outlier3 <- sum(v[valid] < (q1 - 3*iqr) | v[valid] > (q3 + 3*iqr))
    nm <- if (!is.null(rownames(res))) rownames(res)[k] else paste0("row", k)
    cat(sprintf("%-25s: mean=%10.4f, sd=%8.4f, min=%10.4f, max=%10.4f, outliers(3IQR)=%d\n",
                nm, mean(v[valid]), sd(v[valid]), min(v[valid]), max(v[valid]), outlier3))
  }

  # Detailed look at mu_ipw row
  mu_row <- grep("mu_ipw|mu\\.ipw|mu_hat|estimate", rownames(res), ignore.case = TRUE)
  if (length(mu_row) > 0) {
    cat(sprintf("\nDetailed mu row (%s):\n", rownames(res)[mu_row[1]]))
    mu <- res[mu_row[1], ]
    valid <- !is.na(mu)
    mu_v <- mu[valid]
    q1 <- quantile(mu_v, 0.25); q3 <- quantile(mu_v, 0.75); iqr <- q3 - q1
    cat(sprintf("  Q1=%.4f, Q3=%.4f, IQR=%.4f\n", q1, q3, iqr))
    cat(sprintf("  Smallest 10: %s\n", paste(round(head(sort(mu_v), 10), 4), collapse=", ")))
    cat(sprintf("  Largest 10:  %s\n", paste(round(tail(sort(mu_v), 10), 4), collapse=", ")))
  }

  # Check iterations row
  iter_row <- grep("iter", rownames(res), ignore.case = TRUE)
  if (length(iter_row) > 0) {
    cat(sprintf("\nIterations row (%s):\n", rownames(res)[iter_row[1]]))
    iters <- res[iter_row[1], ]
    valid_i <- !is.na(iters)
    cat(sprintf("  mean=%.1f, median=%.0f, max=%.0f\n",
                mean(iters[valid_i]), median(iters[valid_i]), max(iters[valid_i])))
    cat(sprintf("  At max (>=200): %d (%.1f%%)\n",
                sum(iters[valid_i] >= 200), 100*mean(iters[valid_i] >= 200)))

    # Cross with mu outliers
    if (length(mu_row) > 0) {
      mu <- res[mu_row[1], ]
      both_valid <- !is.na(mu) & !is.na(iters)
      mu_bv <- mu[both_valid]; iters_bv <- iters[both_valid]
      q1 <- quantile(mu_bv, 0.25); q3 <- quantile(mu_bv, 0.75); iqr <- q3 - q1
      is_out <- mu_bv < (q1 - 3*iqr) | mu_bv > (q3 + 3*iqr)
      cat(sprintf("\n  Outlier reps: mean_iter=%.1f, at_max=%d/%d\n",
                  mean(iters_bv[is_out]), sum(iters_bv[is_out] >= 200), sum(is_out)))
      cat(sprintf("  Normal reps:  mean_iter=%.1f, at_max=%d/%d\n",
                  mean(iters_bv[!is_out]), sum(iters_bv[!is_out] >= 200), sum(!is_out)))
    }
  }

  # Check grad norm
  gn_row <- grep("grad", rownames(res), ignore.case = TRUE)
  if (length(gn_row) > 0) {
    cat(sprintf("\nGrad norm row (%s):\n", rownames(res)[gn_row[1]]))
    gn <- res[gn_row[1], ]
    valid_g <- !is.na(gn)
    cat(sprintf("  mean=%.6f, median=%.6f, max=%.6f\n",
                mean(gn[valid_g]), median(gn[valid_g]), max(gn[valid_g])))
    cat(sprintf("  >1e-4: %d, >1e-2: %d, >1: %d\n",
                sum(gn[valid_g] > 1e-4), sum(gn[valid_g] > 1e-2), sum(gn[valid_g] > 1)))
  }
}

cat("\nDone!\n")
