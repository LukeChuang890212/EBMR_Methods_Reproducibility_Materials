# AIDS Clinical Trial - Bootstrap Analysis for Arm 1 (ZDV + ddI)
# Parallel computation with progress bar

library(EBMRalgorithmtest)
library(speff2trial)
library(parallel)
library(doSNOW)
library(foreach)

# Load data
data(ACTG175)
dat <- ACTG175

# Prepare data for analysis
dat$y <- dat$cd496
dat$y[dat$r == 0] <- -1

# Subset to Arm 1 only
arm_data <- dat[dat$arms == 1, ]

cat("=============================================================================\n")
cat("Bootstrap Analysis for Arm 1: ZDV + ddI\n")
cat("=============================================================================\n\n")

cat("Sample size:", nrow(arm_data), "\n")
cat("Respondents:", sum(arm_data$r), "(", round(mean(arm_data$r)*100, 1), "%)\n")
cat("Naive estimate (observed mean):", round(mean(arm_data$y[arm_data$r == 1]), 2), "\n\n")

# Define inverse link function
inv_link <- function(eta) 1 / (1 + exp(-eta))

# Define weight matrix function
W <- function(g.matrix) {
  m <- t(g.matrix) %*% g.matrix / nrow(g.matrix)
  tryCatch(solve(m), error = function(e) MASS::ginv(m))
}

# Define h_nu for ensemble estimation
h_nu <- function(data) {
  mat <- cbind(data$cd40, data$cd80, data$cd420, data$cd820, data$age, data$wtkg)
  colnames(mat) <- c("cd40", "cd80", "cd420", "cd820", "age", "wtkg")
  return(mat)
}

# Use same h_alpha for all models
h_alpha1 <- h_alpha2 <- h_alpha3 <- h_nu

# Initial values for alpha
init_vals <- c(0, 0)

# =============================================================================
# Define PS model specifications
# =============================================================================

ps_specs_1 <- list(
  formula.list = list(r ~ o(y) + cd40+cd420+cd820+age),
  h_alpha.list = list(h_alpha1),
  alpha_init.list = list(init_vals),
  inv_link = inv_link
)

ps_specs_2 <- list(
  formula.list = list(r ~ o(y) + cd40+cd820+age+wtkg),
  h_alpha.list = list(h_alpha2),
  alpha_init.list = list(init_vals),
  inv_link = inv_link
)

ps_specs_3 <- list(
  formula.list = list(r ~ o(y) + cd420+cd820+age+wtkg),
  h_alpha.list = list(h_alpha3),
  alpha_init.list = list(init_vals),
  inv_link = inv_link
)

ps_specs_12 <- list(
  formula.list = list(ps_specs_1$formula.list[[1]], ps_specs_2$formula.list[[1]]),
  h_alpha.list = list(h_alpha1, h_alpha2),
  alpha_init.list = list(init_vals, init_vals),
  inv_link = inv_link
)

ps_specs_23 <- list(
  formula.list = list(ps_specs_2$formula.list[[1]], ps_specs_3$formula.list[[1]]),
  h_alpha.list = list(h_alpha2, h_alpha3),
  alpha_init.list = list(init_vals, init_vals),
  inv_link = inv_link
)

ps_specs_123 <- list(
  formula.list = list(ps_specs_1$formula.list[[1]], ps_specs_2$formula.list[[1]], ps_specs_3$formula.list[[1]]),
  h_alpha.list = list(h_alpha1, h_alpha2, h_alpha3),
  alpha_init.list = list(init_vals, init_vals, init_vals),
  inv_link = inv_link
)

ps_specs_list <- list(ps_specs_1, ps_specs_2, ps_specs_3,
                      ps_specs_12, ps_specs_23, ps_specs_123)
combo_names <- c("Model 1 (cd40)", "Model 2 (cd420)", "Model 3 (age)",
                 "Models 1+2", "Models 2+3", "Models 1+2+3")

# =============================================================================
# Point Estimates First
# =============================================================================

cat("=== Computing Point Estimates ===\n\n")

results_list <- list()
for(i in 1:length(ps_specs_list)) {
  cat("Estimating:", combo_names[i], "... ")

  results_list[[i]] <- tryCatch({
    ebmr <- EBMRalgorithmtest::EBMRAlgorithmtest$new(
      y_names = "y",
      ps_specifications = ps_specs_list[[i]],
      data = arm_data,
      W = W
    )

    result <- ebmr$EBMR_IPW(h_nu = h_nu, se.fit = TRUE)
    cat("Success! mu_ipw =", sprintf("%.2f", result$mu_ipw), "\n")

    list(
      mu_ipw = result$mu_ipw,
      se_ipw = result$se_ipw,
      nu.hat = result$nu.hat,
      w.hat = result$w.hat,
      n_models = length(ps_specs_list[[i]]$formula.list)
    )
  }, error = function(e) {
    cat("Error:", conditionMessage(e), "\n")
    list(mu_ipw = NA, se_ipw = NA, nu.hat = NA, w.hat = NA, n_models = NA)
  })
}

# =============================================================================
# Bootstrap Analysis
# =============================================================================

cat("\n")
cat(strrep("=", 80), "\n")
cat("BOOTSTRAP ANALYSIS (B = 1000) - PARALLEL\n")
cat(strrep("=", 80), "\n\n")

set.seed(12345)
B <- 1000
n <- nrow(arm_data)

# Detect number of cores
n_cores <- detectCores() - 1
cat("Using", n_cores, "cores for parallel computation\n\n")

# Function to estimate mu_ipw for a single bootstrap iteration
boot_one_iter <- function(b, arm_data, ps_specs_list, h_nu, W) {
  n <- nrow(arm_data)
  # Bootstrap sample
  boot_idx <- sample(1:n, n, replace = TRUE)
  boot_dat <- arm_data[boot_idx, ]

  # Estimate for each model combination
  estimates <- numeric(6)
  for(i in 1:6) {
    estimates[i] <- tryCatch({
      ebmr <- EBMRalgorithmtest::EBMRAlgorithmtest$new(
        y_names = "y",
        ps_specifications = ps_specs_list[[i]],
        data = boot_dat,
        W = W
      )
      result <- ebmr$EBMR_IPW(h_nu = h_nu, se.fit = FALSE)
      result$mu_ipw
    }, error = function(e) {
      NA
    })
  }
  return(estimates)
}

# Set up cluster with doSNOW for progress bar
cl <- makeCluster(n_cores)
registerDoSNOW(cl)

# Create progress bar
pb <- txtProgressBar(min = 0, max = B, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

cat("Running bootstrap...\n")

# Run parallel bootstrap with progress bar
boot_results <- foreach(b = 1:B, .combine = rbind, .packages = "EBMRalgorithmtest",
                        .options.snow = opts,
                        .export = c("arm_data", "ps_specs_list", "h_nu", "W", "boot_one_iter")) %dopar% {
  boot_one_iter(b, arm_data, ps_specs_list, h_nu, W)
}

close(pb)

# Stop cluster
stopCluster(cl)

# boot_results is already B x 6 matrix from foreach with .combine = rbind
boot_estimates <- boot_results
colnames(boot_estimates) <- combo_names

cat("\nBootstrap complete!\n\n")

# Save bootstrap results
save(boot_estimates, results_list, combo_names, file = "Results/bootstrap_arm1.RData")
cat("Bootstrap results saved to 'Results/bootstrap_arm1.RData'\n\n")

# Calculate bootstrap SE
boot_se <- apply(boot_estimates, 2, sd, na.rm = TRUE)
boot_mean <- apply(boot_estimates, 2, mean, na.rm = TRUE)
boot_valid <- apply(boot_estimates, 2, function(x) sum(!is.na(x)))

# Calculate 95% CI (percentile method)
boot_ci_lower <- apply(boot_estimates, 2, quantile, probs = 0.025, na.rm = TRUE)
boot_ci_upper <- apply(boot_estimates, 2, quantile, probs = 0.975, na.rm = TRUE)

# =============================================================================
# Summary Table
# =============================================================================

cat(strrep("=", 120), "\n")
cat("RESULTS WITH BOOTSTRAP STANDARD ERRORS - ARM 1 (ZDV + ddI)\n")
cat(strrep("=", 120), "\n\n")

cat(sprintf("%-20s %12s %12s %20s %25s %25s %12s\n",
            "Model", "Estimate", "Boot SE", "95% CI", "nu.hat", "w.hat", "Valid Boot"))
cat(strrep("-", 130), "\n")

for(i in 1:6) {
  res <- results_list[[i]]
  if(!is.na(res$mu_ipw)) {
    ci_str <- sprintf("[%.2f, %.2f]", boot_ci_lower[i], boot_ci_upper[i])
    # Format nu.hat and w.hat
    if(length(res$nu.hat) == 1) {
      nu_str <- sprintf("%.4f", res$nu.hat)
      w_str <- sprintf("%.4f", res$w.hat)
    } else {
      nu_str <- paste(sprintf("%.3f", res$nu.hat), collapse = ", ")
      w_str <- paste(sprintf("%.3f", res$w.hat), collapse = ", ")
    }
    cat(sprintf("%-20s %12.2f %12.4f %20s %25s %25s %12d\n",
                combo_names[i], res$mu_ipw, boot_se[i], ci_str, nu_str, w_str, boot_valid[i]))
  } else {
    cat(sprintf("%-20s %12s %12s %20s %25s %25s %12s\n",
                combo_names[i], "ERROR", "N/A", "N/A", "N/A", "N/A", "N/A"))
  }
}

naive_est <- mean(arm_data$y[arm_data$r == 1])
cat(strrep("-", 130), "\n")
cat(sprintf("%-20s %12.2f %12s %20s %25s %25s %12s\n",
            "Naive (observed)", naive_est, "N/A", "N/A", "N/A", "N/A", "N/A"))
cat(strrep("=", 120), "\n")

cat("\nNote: cd496 = CD4 count at 96 weeks (cells/mm^3)\n")
cat("Higher values indicate better immune function.\n")

# =============================================================================
# Save summary to text file
# =============================================================================

sink("Results/bootstrap_arm1_summary.txt")

cat("=============================================================================\n")
cat("Bootstrap Analysis Results - Arm 1: ZDV + ddI\n")
cat("=============================================================================\n\n")

cat("Sample size:", nrow(arm_data), "\n")
cat("Respondents:", sum(arm_data$r), "(", round(mean(arm_data$r)*100, 1), "%)\n")
cat("Naive estimate (observed mean):", round(naive_est, 2), "\n")
cat("Bootstrap replications:", B, "\n")
cat("Cores used:", n_cores, "\n\n")

cat(strrep("=", 120), "\n")
cat("RESULTS WITH BOOTSTRAP STANDARD ERRORS\n")
cat(strrep("=", 120), "\n\n")

cat(sprintf("%-20s %12s %12s %20s %25s %25s %12s\n",
            "Model", "Estimate", "Boot SE", "95% CI", "nu.hat", "w.hat", "Valid Boot"))
cat(strrep("-", 130), "\n")

for(i in 1:6) {
  res <- results_list[[i]]
  if(!is.na(res$mu_ipw)) {
    ci_str <- sprintf("[%.2f, %.2f]", boot_ci_lower[i], boot_ci_upper[i])
    if(length(res$nu.hat) == 1) {
      nu_str <- sprintf("%.4f", res$nu.hat)
      w_str <- sprintf("%.4f", res$w.hat)
    } else {
      nu_str <- paste(sprintf("%.3f", res$nu.hat), collapse = ", ")
      w_str <- paste(sprintf("%.3f", res$w.hat), collapse = ", ")
    }
    cat(sprintf("%-20s %12.2f %12.4f %20s %25s %25s %12d\n",
                combo_names[i], res$mu_ipw, boot_se[i], ci_str, nu_str, w_str, boot_valid[i]))
  } else {
    cat(sprintf("%-20s %12s %12s %20s %25s %25s %12s\n",
                combo_names[i], "ERROR", "N/A", "N/A", "N/A", "N/A", "N/A"))
  }
}

cat(strrep("-", 130), "\n")
cat(sprintf("%-20s %12.2f %12s %20s %25s %25s %12s\n",
            "Naive (observed)", naive_est, "N/A", "N/A", "N/A", "N/A", "N/A"))
cat(strrep("=", 120), "\n")

cat("\nNote: cd496 = CD4 count at 96 weeks (cells/mm^3)\n")
cat("Higher values indicate better immune function.\n")

sink()

cat("\nSummary saved to 'Results/bootstrap_arm1_summary.txt'\n")
