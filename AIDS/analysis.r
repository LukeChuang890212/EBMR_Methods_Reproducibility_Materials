# AIDS Clinical Trial Data Analysis (ACTG175)
# Estimate population mean of cd496 (CD4 count at 96 weeks) using EBMRalgorithmtest
# Analysis by treatment arm

library(EBMRalgorithmtest)
library(speff2trial)

# Load data
data(ACTG175)
dat <- ACTG175

cat("=============================================================================\n")
cat("ACTG175 AIDS Clinical Trial Data Analysis\n")
cat("Goal: Estimate population mean of CD4 count at 96 weeks (cd496)\n")
cat("=============================================================================\n\n")

# Data summary
cat("=== Data Summary ===\n")
cat("Total observations:", nrow(dat), "\n")
cat("Respondents (r=1):", sum(dat$r), "\n")
cat("Non-respondents (r=0):", sum(dat$r == 0), "\n")
cat("Response rate:", round(mean(dat$r) * 100, 2), "%\n\n")

# Treatment arms
cat("=== Treatment Arms ===\n")
cat("arms = 0: zidovudine (ZDV) only\n")
cat("arms = 1: ZDV + didanosine (ddI)\n")
cat("arms = 2: ZDV + zalcitabine (ddC)\n")
cat("arms = 3: didanosine (ddI) only\n\n")

# Summary by arm
cat("=== Summary by Treatment Arm ===\n")
for(arm in 0:3) {
  arm_dat <- dat[dat$arms == arm, ]
  cat(sprintf("Arm %d: n=%d, respondents=%d (%.1f%%), mean cd496 (observed)=%.1f\n",
              arm, nrow(arm_dat), sum(arm_dat$r),
              mean(arm_dat$r) * 100,
              mean(arm_dat$cd496[arm_dat$r == 1], na.rm = TRUE)))
}
cat("\n")

# Prepare data for analysis
# Set y = cd496 for respondents, y = -1 for non-respondents (placeholder)
dat$y <- dat$cd496
dat$y[dat$r == 0] <- -1

# Define inverse link function
inv_link <- function(eta) 1 / (1 + exp(-eta))

# Define weight matrix function
W <- function(g.matrix) {
  m <- t(g.matrix) %*% g.matrix / nrow(g.matrix)
  tryCatch(solve(m), error = function(e) MASS::ginv(m))
}

# Define h_nu for ensemble estimation
# Use baseline CD4 counts and key covariates
h_nu <- function(data) {
  mat <- cbind(data$cd40, data$cd80, data$cd420, data$cd820, data$age, data$wtkg)
  colnames(mat) <- c("cd40", "cd80", "cd420", "cd820", "age", "wtkg")
  return(mat)
}

# Define h_alpha for PS model estimation
h_alpha1 <- function(data) {
  mat <- cbind(data$cd40, data$cd80, data$cd420, data$cd820, data$age, data$wtkg)
  colnames(mat) <- c("cd40", "cd80", "cd420", "cd820", "age", "wtkg")
  return(mat)
}

h_alpha2 <- function(data) {
  mat <- cbind(data$cd40,  data$cd80, data$cd820, data$age, data$wtkg)
  colnames(mat) <- c("cd40", "cd80", "cd820", "age", "wtkg")
  return(mat)
}

h_alpha3 <- function(data) {
  mat <- cbind(data$cd80, data$cd420, data$cd820, data$age, data$wtkg)
  colnames(mat) <- c("cd80", "cd420", "cd820", "age", "wtkg")
  return(mat)
}

h_alpha1 = h_alpha2 = h_alpha3 = h_nu 

# Initial values for alpha: (y_coef, x_coef) - intercept is added automatically
init_vals <- c(0, 0)

# =============================================================================
# Define PS model specifications
# =============================================================================

# Model 1: r ~ o(y) + cd40 (baseline CD4)
ps_specs_1 <- list(
  formula.list = list(r ~ o(y) + cd40+cd420+cd820+age),
  h_alpha.list = list(h_alpha1),
  alpha_init.list = list(init_vals),
  inv_link = inv_link
)

# Model 2: r ~ o(y) + cd420 (CD4 at 20 weeks)
ps_specs_2 <- list(
  formula.list = list(r ~ o(y) + cd40+cd820+age+wtkg),
  h_alpha.list = list(h_alpha2),
  alpha_init.list = list(init_vals),
  inv_link = inv_link
)

# Model 3: r ~ o(y) + age
ps_specs_3 <- list(
  formula.list = list(r ~ o(y) + cd420+cd820+age+wtkg),
  h_alpha.list = list(h_alpha3),
  alpha_init.list = list(init_vals),
  inv_link = inv_link
)

# Combination: Models 1 + 2
ps_specs_12 <- list(
  formula.list = list(
    ps_specs_1$formula.list[[1]],
    ps_specs_2$formula.list[[1]]
  ),
  h_alpha.list = list(h_alpha1, h_alpha2),
  alpha_init.list = list(init_vals, init_vals),
  inv_link = inv_link
)

# Combination: Models 2 + 3
ps_specs_23 <- list(
  formula.list = list(
    ps_specs_2$formula.list[[1]],
    ps_specs_3$formula.list[[1]]
  ),
  h_alpha.list = list(h_alpha2, h_alpha3),
  alpha_init.list = list(init_vals, init_vals),
  inv_link = inv_link
)

# Combination: Models 1 + 2 + 3
ps_specs_123 <- list(
  formula.list = list(
    ps_specs_1$formula.list[[1]],
    ps_specs_2$formula.list[[1]],
    ps_specs_3$formula.list[[1]]
  ),
  h_alpha.list = list(h_alpha1, h_alpha2, h_alpha3),
  alpha_init.list = list(init_vals, init_vals, init_vals),
  inv_link = inv_link
)

ps_specs_list <- list(ps_specs_1, ps_specs_2, ps_specs_3,
                      ps_specs_12, ps_specs_23, ps_specs_123)
combo_names <- c("Model 1 (cd40)", "Model 2 (cd420)", "Model 3 (age)",
                 "Models 1+2", "Models 2+3", "Models 1+2+3")

# =============================================================================
# Function to run analysis for a single arm
# =============================================================================
run_arm_analysis <- function(arm_data, arm_name) {
  cat(strrep("=", 80), "\n")
  cat("Analysis for", arm_name, "\n")
  cat(strrep("=", 80), "\n\n")

  cat("Sample size:", nrow(arm_data), "\n")
  cat("Respondents:", sum(arm_data$r), "(", round(mean(arm_data$r)*100, 1), "%)\n")
  cat("Naive estimate (observed mean):", round(mean(arm_data$y[arm_data$r == 1]), 2), "\n\n")

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

  # Summary table
  cat("\n")
  cat(sprintf("%-20s %12s %12s %15s %25s %25s\n", "Model", "Estimate", "SE", "95% CI", "nu.hat", "w.hat"))
  cat(strrep("-", 120), "\n")

  for(i in 1:length(results_list)) {
    res <- results_list[[i]]
    if(!is.na(res$mu_ipw)) {
      ci_lower <- res$mu_ipw - 1.96 * res$se_ipw
      ci_upper <- res$mu_ipw + 1.96 * res$se_ipw
      ci_str <- sprintf("[%.1f, %.1f]", ci_lower, ci_upper)
      # Format nu.hat and w.hat
      if(length(res$nu.hat) == 1) {
        nu_str <- sprintf("%.4f", res$nu.hat)
        w_str <- sprintf("%.4f", res$w.hat)
      } else {
        nu_str <- paste(sprintf("%.3f", res$nu.hat), collapse = ", ")
        w_str <- paste(sprintf("%.3f", res$w.hat), collapse = ", ")
      }
      cat(sprintf("%-20s %12.2f %12.2f %15s %25s %25s\n",
                  combo_names[i], res$mu_ipw, res$se_ipw, ci_str, nu_str, w_str))
    } else {
      cat(sprintf("%-20s %12s %12s %15s %25s %25s\n", combo_names[i], "ERROR", "N/A", "N/A", "N/A", "N/A"))
    }
  }

  naive_est <- mean(arm_data$y[arm_data$r == 1])
  cat(strrep("-", 60), "\n")
  cat(sprintf("%-20s %12.2f %12s %12s\n", "Naive (observed)", naive_est, "N/A", "N/A"))
  cat("\n")

  return(results_list)
}

# =============================================================================
# Run analysis for each treatment arm
# =============================================================================

results_by_arm <- list()

for(arm in 1) {
  arm_name <- switch(as.character(arm),
                     "0" = "Arm 0: ZDV only",
                     "1" = "Arm 1: ZDV + ddI",
                     "2" = "Arm 2: ZDV + ddC",
                     "3" = "Arm 3: ddI only")

  arm_data <- dat[dat$arms == arm, ]
  results_by_arm[[arm + 1]] <- run_arm_analysis(arm_data, arm_name)
}

# =============================================================================
# Summary comparison across arms
# =============================================================================

cat("\n")
cat(strrep("=", 100), "\n")
cat("SUMMARY: COMPARISON ACROSS TREATMENT ARMS\n")
cat(strrep("=", 100), "\n\n")

# Use the ensemble model (Models 1+2+3) for comparison
cat("Using 'Models 1+2+3' ensemble estimator:\n\n")
cat(sprintf("%-20s %12s %12s %12s %20s\n", "Treatment Arm", "n", "Estimate", "SE", "95% CI"))
cat(strrep("-", 80), "\n")

arm_names <- c("Arm 0: ZDV only", "Arm 1: ZDV + ddI", "Arm 2: ZDV + ddC", "Arm 3: ddI only")
for(arm in 1) {
  arm_data <- dat[dat$arms == arm, ]
  res <- results_by_arm[[arm + 1]][[6]]  # Models 1+2+3
  if(!is.na(res$mu_ipw)) {
    ci_lower <- res$mu_ipw - 1.96 * res$se_ipw
    ci_upper <- res$mu_ipw + 1.96 * res$se_ipw
    ci_str <- sprintf("[%.1f, %.1f]", ci_lower, ci_upper)
    cat(sprintf("%-20s %12d %12.2f %12.2f %20s\n",
                arm_names[arm + 1], nrow(arm_data), res$mu_ipw, res$se_ipw, ci_str))
  } else {
    cat(sprintf("%-20s %12d %12s %12s %20s\n",
                arm_names[arm + 1], nrow(arm_data), "ERROR", "N/A", "N/A"))
  }
}

cat(strrep("-", 80), "\n")

# Naive estimates for comparison
cat("\nNaive estimates (observed mean only):\n")
cat(sprintf("%-20s %12s %12s\n", "Treatment Arm", "n (obs)", "Mean cd496"))
cat(strrep("-", 50), "\n")
for(arm in 0:3) {
  arm_data <- dat[dat$arms == arm, ]
  naive_est <- mean(arm_data$y[arm_data$r == 1])
  cat(sprintf("%-20s %12d %12.2f\n",
              arm_names[arm + 1], sum(arm_data$r), naive_est))
}

cat("\n")
cat("Note: cd496 = CD4 count at 96 weeks (cells/mm^3)\n")
cat("Higher values indicate better immune function.\n")
