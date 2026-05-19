# Diagnose outliers in scenario9-3 test57 results
f <- "Simulation_Results/EBMR_IPW_setting3-miss50-scenario9-3_1_n2000_replicate1000_test57.RDS"
res <- readRDS(f)

cat("File:", f, "\n")
cat("Class:", class(res), "\n")
cat("Names:", paste(names(res), collapse=", "), "\n\n")

# Check structure
if (is.list(res)) {
  for (nm in names(res)) {
    val <- res[[nm]]
    if (is.numeric(val) && length(val) > 1) {
      valid <- !is.na(val)
      cat(sprintf("%-20s: n=%d, valid=%d, mean=%.4f, sd=%.4f, min=%.4f, max=%.4f\n",
                  nm, length(val), sum(valid), mean(val[valid]), sd(val[valid]),
                  min(val[valid]), max(val[valid])))
    } else if (is.numeric(val)) {
      cat(sprintf("%-20s: %.6f\n", nm, val))
    } else {
      cat(sprintf("%-20s: %s\n", nm, paste(head(val), collapse=", ")))
    }
  }
}

# Focus on mu_ipw
if ("mu_ipw" %in% names(res)) {
  mu <- res$mu_ipw
  valid <- !is.na(mu)
  mu_v <- mu[valid]
  cat(sprintf("\nmu_ipw: %d valid out of %d\n", sum(valid), length(mu)))

  q1 <- quantile(mu_v, 0.25); q3 <- quantile(mu_v, 0.75); iqr <- q3 - q1
  outlier3 <- mu_v < (q1 - 3*iqr) | mu_v > (q3 + 3*iqr)
  outlier2 <- mu_v < (q1 - 2*iqr) | mu_v > (q3 + 2*iqr)
  cat(sprintf("Outliers (3*IQR): %d (%.1f%%)\n", sum(outlier3), 100*mean(outlier3)))
  cat(sprintf("Outliers (2*IQR): %d (%.1f%%)\n", sum(outlier2), 100*mean(outlier2)))

  cat(sprintf("\nQuantiles: 1%%=%.4f, 5%%=%.4f, 25%%=%.4f, 50%%=%.4f, 75%%=%.4f, 95%%=%.4f, 99%%=%.4f\n",
              quantile(mu_v, 0.01), quantile(mu_v, 0.05), q1, median(mu_v), q3,
              quantile(mu_v, 0.95), quantile(mu_v, 0.99)))

  # Show extreme values
  sorted <- sort(mu_v)
  cat("\nSmallest 10:", paste(round(head(sorted, 10), 4), collapse=", "), "\n")
  cat("Largest 10: ", paste(round(tail(sorted, 10), 4), collapse=", "), "\n")
}

# Check convergence info if available
if ("converged" %in% names(res)) {
  conv <- res$converged
  cat(sprintf("\nConvergence: %d/%d (%.1f%%)\n", sum(conv, na.rm=TRUE), sum(!is.na(conv)),
              100*mean(conv, na.rm=TRUE)))
}
if ("iterations" %in% names(res)) {
  iters <- res$iterations
  valid_i <- !is.na(iters)
  cat(sprintf("Iterations: mean=%.1f, median=%.0f, max=%d\n",
              mean(iters[valid_i]), median(iters[valid_i]), max(iters[valid_i])))

  # How many hit the max (200)?
  at_max <- sum(iters[valid_i] >= 200)
  cat(sprintf("At max iterations (>=200): %d (%.1f%%)\n", at_max, 100*at_max/sum(valid_i)))

  # Cross-tabulate: are outliers the ones that hit max iterations?
  if ("mu_ipw" %in% names(res)) {
    mu <- res$mu_ipw
    valid_both <- !is.na(mu) & !is.na(iters)
    mu_vb <- mu[valid_both]; iters_vb <- iters[valid_both]
    q1 <- quantile(mu_vb, 0.25); q3 <- quantile(mu_vb, 0.75); iqr <- q3 - q1
    is_outlier <- mu_vb < (q1 - 3*iqr) | mu_vb > (q3 + 3*iqr)

    cat(sprintf("\nOutlier vs iterations:\n"))
    cat(sprintf("  Outliers:     mean_iter=%.1f, median_iter=%.0f, max_iter=%d, at_max=%d/%d\n",
                mean(iters_vb[is_outlier]), median(iters_vb[is_outlier]), max(iters_vb[is_outlier]),
                sum(iters_vb[is_outlier] >= 200), sum(is_outlier)))
    cat(sprintf("  Non-outliers: mean_iter=%.1f, median_iter=%.0f, max_iter=%d, at_max=%d/%d\n",
                mean(iters_vb[!is_outlier]), median(iters_vb[!is_outlier]), max(iters_vb[!is_outlier]),
                sum(iters_vb[!is_outlier] >= 200), sum(!is_outlier)))
  }
}

# Check grad norms if available
if ("final_grad_norm" %in% names(res)) {
  fgn <- res$final_grad_norm
  valid_g <- !is.na(fgn)
  cat(sprintf("\nFinal grad norm: mean=%.6f, median=%.6f, max=%.6f\n",
              mean(fgn[valid_g]), median(fgn[valid_g]), max(fgn[valid_g])))

  if ("mu_ipw" %in% names(res)) {
    mu <- res$mu_ipw
    valid_both <- !is.na(mu) & !is.na(fgn)
    mu_vb <- mu[valid_both]; fgn_vb <- fgn[valid_both]
    q1 <- quantile(mu_vb, 0.25); q3 <- quantile(mu_vb, 0.75); iqr <- q3 - q1
    is_outlier <- mu_vb < (q1 - 3*iqr) | mu_vb > (q3 + 3*iqr)

    cat(sprintf("\nGrad norm vs outlier:\n"))
    cat(sprintf("  Outliers:     mean=%.6f, median=%.6f, max=%.6f\n",
                mean(fgn_vb[is_outlier]), median(fgn_vb[is_outlier]), max(fgn_vb[is_outlier])))
    cat(sprintf("  Non-outliers: mean=%.6f, median=%.6f, max=%.6f\n",
                mean(fgn_vb[!is_outlier]), median(fgn_vb[!is_outlier]), max(fgn_vb[!is_outlier])))
  }
}

# Check alpha coefficients if available
if ("alpha" %in% names(res)) {
  alpha <- res$alpha
  if (is.matrix(alpha)) {
    cat(sprintf("\nAlpha matrix: %d x %d\n", nrow(alpha), ncol(alpha)))
    for (k in 1:ncol(alpha)) {
      cat(sprintf("  alpha%d: mean=%.4f, sd=%.4f, range=[%.4f, %.4f]\n",
                  k, mean(alpha[,k], na.rm=T), sd(alpha[,k], na.rm=T),
                  min(alpha[,k], na.rm=T), max(alpha[,k], na.rm=T)))
    }
  } else if (is.list(alpha)) {
    cat(sprintf("\nAlpha: list of %d elements\n", length(alpha)))
    alpha_mat <- do.call(rbind, alpha[!sapply(alpha, is.null)])
    cat(sprintf("Alpha matrix: %d x %d\n", nrow(alpha_mat), ncol(alpha_mat)))
    for (k in 1:ncol(alpha_mat)) {
      cat(sprintf("  alpha%d: mean=%.4f, sd=%.4f, range=[%.4f, %.4f]\n",
                  k, mean(alpha_mat[,k], na.rm=T), sd(alpha_mat[,k], na.rm=T),
                  min(alpha_mat[,k], na.rm=T), max(alpha_mat[,k], na.rm=T)))
    }
  }
}

cat("\nDone!\n")
