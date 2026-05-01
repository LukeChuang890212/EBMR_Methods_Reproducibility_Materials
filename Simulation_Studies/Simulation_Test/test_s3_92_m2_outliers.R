## Investigate outliers in setting3, scenario 9-2, model 2, n=500, miss30
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

# Load latest result
sim_file <- "Simulation_Results/EBMR_IPW_setting3-miss30-scenario9-2_2_n500_replicate1000_test59.RDS"
sim_result <- readRDS(sim_file)
mu_true <- get_mu_true("setting3")

cat("=== Investigating outliers: Setting 3, Scenario 9-2, Model 2, n=500, miss30 ===\n\n")
cat(sprintf("True mean: %.4f\n", mu_true))
cat(sprintf("Dimensions: %d rows x %d reps\n\n", nrow(sim_result), ncol(sim_result)))

# Row structure: mu_ipw(1), mu_ipw.true(2), se_ipw(3), se_ipw.true(4), ...
mu_ipw <- sim_result[1, ]
se_ipw <- sim_result[3, ]

# NA check
n_na <- sum(is.na(mu_ipw) | is.na(se_ipw))
cat(sprintf("NAs: %d\n", n_na))

# Remove NAs for analysis
valid <- !is.na(mu_ipw) & !is.na(se_ipw)
mu_v <- mu_ipw[valid]; se_v <- se_ipw[valid]

# Outlier detection (IQR x 3)
Q1 <- quantile(mu_v, 0.25); Q3 <- quantile(mu_v, 0.75)
IQR_val <- Q3 - Q1
is_outlier <- mu_v < (Q1 - 3 * IQR_val) | mu_v > (Q3 + 3 * IQR_val)
cat(sprintf("Outliers (IQR x 3): %d / %d\n", sum(is_outlier), length(mu_v)))

# Distribution of mu
cat(sprintf("\nmu_ipw distribution:\n"))
cat(sprintf("  min=%.4f, Q1=%.4f, med=%.4f, Q3=%.4f, max=%.4f\n",
    min(mu_v), Q1, median(mu_v), Q3, max(mu_v)))
cat(sprintf("  mean=%.4f, sd=%.4f\n", mean(mu_v), sd(mu_v)))
cat(sprintf("  1%%=%.4f, 5%%=%.4f, 95%%=%.4f, 99%%=%.4f\n",
    quantile(mu_v, 0.01), quantile(mu_v, 0.05), quantile(mu_v, 0.95), quantile(mu_v, 0.99)))

# Show outlier values
if (sum(is_outlier) > 0) {
  out_vals <- sort(mu_v[is_outlier])
  cat(sprintf("\nOutlier mu values (%d total):\n", length(out_vals)))
  if (length(out_vals) <= 30) {
    cat("  ", paste(round(out_vals, 3), collapse=", "), "\n")
  } else {
    cat("  First 15: ", paste(round(head(out_vals, 15), 3), collapse=", "), "\n")
    cat("  Last  15: ", paste(round(tail(out_vals, 15), 3), collapse=", "), "\n")
  }
  cat(sprintf("  Lower bound: %.4f, Upper bound: %.4f\n", Q1 - 3*IQR_val, Q3 + 3*IQR_val))
}

# Check if outliers have extreme SE
cat(sprintf("\nSE for outlier reps: mean=%.4f, med=%.4f\n",
    mean(se_v[is_outlier]), median(se_v[is_outlier])))
cat(sprintf("SE for normal reps:  mean=%.4f, med=%.4f\n",
    mean(se_v[!is_outlier]), median(se_v[!is_outlier])))

# Now investigate: what's causing outliers?
# Load the data and refit problematic reps to check alpha, PS, etc.
cat("\n=== Investigating root cause ===\n\n")

ps_spec <- get_ps_spec("9")  # scenario 9-2 uses ps_spec "9"
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
n_val <- 500
data_file <- misspecified_model_all_data_file.list[["setting3"]][["miss30"]]

# Find correct data file for n=500
cat("Data files for setting3 miss30:\n")
for (k in seq_along(data_file)) {
  cat(sprintf("  [%d]: %s\n", k, data_file[[k]]))
}

# Determine which file index has n=500
n_vectors <- n.vector.list$misspecified_model
cat("\nn.vector.list:\n")
print(n_vectors)

# Find outlier rep indices
outlier_idx <- which(is_outlier)
normal_idx <- which(!is_outlier)

# Check a few outlier reps and normal reps
check_reps <- c(head(outlier_idx[order(abs(mu_v[outlier_idx] - median(mu_v)), decreasing=TRUE)], 5),
                head(normal_idx, 3))

cat(sprintf("\nChecking %d reps (%d outlier, %d normal):\n", length(check_reps),
    min(5, sum(is_outlier)), 3))

# Load data
all_data_files <- misspecified_model_all_data_file.list[["setting3"]][["miss30"]]
# Try to find the right file - n=500 data
for (fi in seq_along(all_data_files)) {
  if (file.exists(all_data_files[[fi]])) {
    test_data <- readRDS(all_data_files[[fi]])
    n_per_rep <- nrow(test_data) / 1000
    if (n_per_rep == 500) {
      cat(sprintf("Using data file [%d] with n=%d per rep\n", fi, n_per_rep))
      all_data <- test_data
      break
    }
  }
}

m_idx <- 2  # Model 2: u1_z2 = r ~ y + u1 + z2
single_ps <- list(
  formula.list = list(ps_spec[["formula.list"]][[m_idx]]),
  h_alpha.list = list(ps_spec[["h_alpha.list"]][[m_idx]]),
  inv_link = ps_spec[["inv_link"]], outcome = ps_spec[["outcome"]]
)

for (rep_i in check_reps) {
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
  tag <- if (rep_i %in% outlier_idx) "OUTLIER" else "NORMAL"

  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
    pf <- ebmr$ps_fit.list[[1]]
    ps <- pf$fitted.values
    alpha <- pf$coefficients
    mu <- mean(dat$r * dat$y / ps)

    # Condition number
    g_mat <- (dat$r / ps - 1) * pf$h_x
    W_hat <- tryCatch(solve(crossprod(g_mat) / n_val), error = function(e) diag(ncol(pf$h_x)))
    Gamma <- crossprod(pf$h_x * (dat$r * (1-ps) / ps), pf$design_matrix) / n_val
    H <- crossprod(Gamma, W_hat %*% Gamma)
    eig <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    cond <- max(eig) / max(min(eig), 1e-15)

    cat(sprintf("[%s] Rep %3d: mu=%.4f, alpha=(%s), PS=[%.4f,%.4f], >0.99:%d, <0.01:%d, cond=%.1e\n",
        tag, rep_i, mu,
        paste(round(alpha, 3), collapse=", "),
        min(ps), max(ps), sum(ps > 0.99), sum(ps < 0.01), cond))
  }, error = function(e) {
    cat(sprintf("[%s] Rep %3d: ERROR: %s\n", tag, rep_i, e$message))
  })
}
