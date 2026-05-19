#------------------------------------------------------------------------------#
# Investigate Why MNAR Coefficient on Y is Not Significant
# Focus on the actual variation and identification
#------------------------------------------------------------------------------#

library(EBMRalgorithmFast)
library(tidyverse)

# Source utility functions
source("MHD_functions.R")

#------------------------------------------------------------------------------#
# Read and prepare data
#------------------------------------------------------------------------------#
original_data <- read.csv("data_application.csv")
percent <- original_data$Percentage
n <- 2486
class_n <- round(n * percent / 100)

dat <- gen_data(original_data, class_n, n)

cat("================================================================================\n")
cat("     INVESTIGATION: WHY IS MNAR COEFFICIENT ON Y NOT SIGNIFICANT?              \n")
cat("================================================================================\n\n")

#------------------------------------------------------------------------------#
# Check data structure among OBSERVED cases only
#------------------------------------------------------------------------------#
cat("--------------------------------------------------------------------------------\n")
cat("DATA STRUCTURE AMONG OBSERVED CASES (r = 1)\n")
cat("--------------------------------------------------------------------------------\n\n")

for (health_val in 0:1) {
  subdat <- dat[dat$health == health_val, ]
  obs_dat <- subdat[subdat$r == 1, ]

  cat("Health =", health_val, "\n")
  cat("  Total n:", nrow(subdat), "\n")
  cat("  Observed n:", nrow(obs_dat), "\n")
  cat("  Missing n:", sum(subdat$r == 0), "\n")
  cat("  Missing rate:", round(100 * mean(subdat$r == 0), 1), "%\n\n")

  # Distribution of y among observed
  cat("  Distribution of Y among observed:\n")
  y_tab <- table(obs_dat$teacher_report)
  print(y_tab)
  cat("  Proportions:", round(prop.table(y_tab), 3), "\n\n")

  # Cross-tabulation of y with covariates among observed
  cat("  Cross-tab of Y and Father (among observed):\n")
  print(table(Y = obs_dat$teacher_report, Father = obs_dat$father))
  cat("\n")

  cat("  Cross-tab of Y and Parent_report (among observed):\n")
  print(table(Y = obs_dat$teacher_report, Parent = obs_dat$parent_report))
  cat("\n\n")
}

#------------------------------------------------------------------------------#
# Check the moment conditions more carefully
#------------------------------------------------------------------------------#
cat("--------------------------------------------------------------------------------\n")
cat("MOMENT CONDITIONS ANALYSIS\n")
cat("--------------------------------------------------------------------------------\n\n")

# The moment conditions are: E[(r/pi - 1) * h(x)] = 0
# where pi = inv_link(alpha_0 + alpha_1 * y + alpha_2 * x)
#
# For observed (r=1): contribution is (1/pi - 1) * h(x)
# For missing (r=0): contribution is (0/pi - 1) * h(x) = -h(x)
#
# The identification comes from: among observed cases, does y predict pi?

cat("Key insight: The moment conditions are:\n")
cat("  E[(r/pi - 1) * h(x)] = 0\n\n")
cat("For r=1 (observed): (1/pi - 1) * h(x) = ((1-pi)/pi) * h(x)\n")
cat("For r=0 (missing):  (0/pi - 1) * h(x) = -h(x)\n\n")

cat("The alpha_1 (coef on y) is identified through variation in y AMONG OBSERVED.\n")
cat("But identification also requires that h(x) varies with y.\n\n")

#------------------------------------------------------------------------------#
# Check covariance structure
#------------------------------------------------------------------------------#
cat("--------------------------------------------------------------------------------\n")
cat("COVARIANCE STRUCTURE AMONG OBSERVED\n")
cat("--------------------------------------------------------------------------------\n\n")

for (health_val in 0:1) {
  subdat <- dat[dat$health == health_val, ]
  obs_dat <- subdat[subdat$r == 1, ]

  cat("Health =", health_val, "\n")

  # Correlation between y and h(x) variables among observed
  cat("  Correlation matrix (Y, Father, Parent_report) among observed:\n")
  cor_mat <- cor(obs_dat[, c("teacher_report", "father", "parent_report")])
  print(round(cor_mat, 4))
  cat("\n")

  # Mean of h(x) by y value
  cat("  Mean of covariates by Y value (among observed):\n")
  for (y_val in unique(obs_dat$teacher_report)) {
    y_subset <- obs_dat[obs_dat$teacher_report == y_val, ]
    cat(sprintf("    Y = %d: father = %.3f, parent_report = %.3f, n = %d\n",
                y_val, mean(y_subset$father), mean(y_subset$parent_report), nrow(y_subset)))
  }
  cat("\n")
}

#------------------------------------------------------------------------------#
# Fit MNAR models and examine in detail
#------------------------------------------------------------------------------#
cat("--------------------------------------------------------------------------------\n")
cat("DETAILED MNAR MODEL ANALYSIS\n")
cat("--------------------------------------------------------------------------------\n\n")

full_ps_specifications <- list(
  formula.list = list(
    r ~ o(teacher_report) + father,
    r ~ o(teacher_report) + parent_report,
    r ~ father + parent_report
  ),
  h_x_names.list = list(
    c("father", "parent_report"),
    c("father", "parent_report"),
    c("father", "parent_report")
  ),
  inv_link = function(eta) 1 / (1 + exp(eta))
)

W <- function(g.matrix) {
  solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
}

init.list <- list(
  alpha_init.list0 = list(
    c(-0.3834341, 1.07945126, -0.1318883),
    c(-4.2670736, 5.9582123, -0.5972129),
    c(0.5, -0.2, -0.5)
  ),
  alpha_init.list1 = list(
    c(0.09019492, -1.511548, -0.3946616),
    c(-1.969056, 3.5823611, -1.37232734),
    c(0.5, -0.2, -0.5)
  )
)

model_names <- c("Model 1 (o(y) + father)", "Model 2 (o(y) + parent_report)")

for (health_val in 0:1) {
  subdat <- dat[dat$health == health_val, ]

  cat("================================================================================\n")
  cat("Health =", health_val, "\n")
  cat("================================================================================\n\n")

  ps_specifications <- list(
    formula.list = full_ps_specifications$formula.list[1:2],
    h_x_names.list = full_ps_specifications$h_x_names.list[1:2],
    alpha_init.list = init.list[[health_val + 1]][1:2],
    inv_link = full_ps_specifications$inv_link
  )

  ebmr <- EBMRAlgorithmFast$new("teacher_report", ps_specifications, subdat, W)

  for (j in 1:2) {
    ps_fit <- ebmr$ps_fit.list[[j]]

    cat(model_names[j], "\n")
    cat(paste(rep("-", 60), collapse = ""), "\n")

    alpha <- ps_fit$coefficients
    se <- ps_fit$se

    cat("Coefficients:\n")
    cat(sprintf("  Intercept:     %.4f (SE: %.4f, z: %.2f, p: %.4f)\n",
                alpha[1], se[1], alpha[1]/se[1], 2*pnorm(-abs(alpha[1]/se[1]))))
    cat(sprintf("  o(y):          %.4f (SE: %.4f, z: %.2f, p: %.4f)\n",
                alpha[2], se[2], alpha[2]/se[2], 2*pnorm(-abs(alpha[2]/se[2]))))
    cat(sprintf("  covariate:     %.4f (SE: %.4f, z: %.2f, p: %.4f)\n",
                alpha[3], se[3], alpha[3]/se[3], 2*pnorm(-abs(alpha[3]/se[3]))))

    # Propensity score analysis
    pi_hat <- ps_fit$fitted.values
    cat("\nPropensity score summary:\n")
    cat(sprintf("  Range: [%.4f, %.4f]\n", min(pi_hat), max(pi_hat)))
    cat(sprintf("  Mean:  %.4f\n", mean(pi_hat)))
    cat(sprintf("  SD:    %.4f\n", sd(pi_hat)))

    # PS by response status
    cat("\nPS by response status:\n")
    cat(sprintf("  Mean PS (r=1): %.4f\n", mean(pi_hat[subdat$r == 1])))
    cat(sprintf("  Mean PS (r=0): %.4f\n", mean(pi_hat[subdat$r == 0])))

    # PS by y value (among observed)
    obs_idx <- subdat$r == 1
    cat("\nPS by Y value (among observed):\n")
    for (y_val in sort(unique(subdat$teacher_report[obs_idx]))) {
      y_idx <- obs_idx & (subdat$teacher_report == y_val)
      cat(sprintf("  Y = %d: mean PS = %.4f, n = %d\n",
                  y_val, mean(pi_hat[y_idx]), sum(y_idx)))
    }

    # Check moment conditions at solution
    g_matrix <- ps_fit$gmm_fit$g.matrix
    cat("\nMoment conditions at solution:\n")
    cat(sprintf("  E[g] = (%.6f, %.6f, %.6f)\n",
                mean(g_matrix[,1]), mean(g_matrix[,2]), mean(g_matrix[,3])))

    # Convergence info
    cat("\nConvergence:\n")
    cat(sprintf("  Code: %d\n", ps_fit$gmm_fit$opt$convergence))
    cat(sprintf("  Objective: %.6e\n", ps_fit$gmm_fit$opt$value))

    cat("\n\n")
  }
}

#------------------------------------------------------------------------------#
# Compare with what we would expect if MNAR is true
#------------------------------------------------------------------------------#
cat("================================================================================\n")
cat("EXPECTED MNAR PATTERN\n")
cat("================================================================================\n\n")

cat("If the data is truly MNAR with y influencing missingness:\n")
cat("  - The coefficient on o(y) should be non-zero\n")
cat("  - PS should vary with y among observed cases\n")
cat("  - The direction depends on whether higher y -> more or less likely to respond\n\n")

cat("What we observe:\n")
cat("  - PS does vary with y (see above)\n")
cat("  - But the SE is large, making the coefficient not significant\n\n")

cat("Possible explanations:\n")
cat("  1. Limited variation in y (binary: 0/1)\n")
cat("  2. Weak correlation between y and h(x) instruments\n")
cat("  3. Small effective sample size for identification\n")
cat("  4. Model may be weakly identified\n\n")

#------------------------------------------------------------------------------#
# Check effective sample for identification
#------------------------------------------------------------------------------#
cat("--------------------------------------------------------------------------------\n")
cat("EFFECTIVE SAMPLE FOR IDENTIFICATION\n")
cat("--------------------------------------------------------------------------------\n\n")

cat("The MNAR coefficient alpha_1 is identified through:\n")
cat("  - Variation in (1/pi - 1) * h(x) as a function of y\n")
cat("  - Only observed cases (r=1) contribute to this variation\n\n")

for (health_val in 0:1) {
  subdat <- dat[dat$health == health_val, ]
  obs_dat <- subdat[subdat$r == 1, ]

  cat("Health =", health_val, ":\n")
  cat("  n observed:", nrow(obs_dat), "\n")
  cat("  n with y=0:", sum(obs_dat$teacher_report == 0), "\n")
  cat("  n with y=1:", sum(obs_dat$teacher_report == 1), "\n")

  # Unique covariate patterns by y
  cat("  Unique (father, parent_report) patterns:\n")
  for (y_val in 0:1) {
    y_sub <- obs_dat[obs_dat$teacher_report == y_val, ]
    n_patterns <- nrow(unique(y_sub[, c("father", "parent_report")]))
    cat(sprintf("    Y = %d: %d patterns\n", y_val, n_patterns))
  }
  cat("\n")
}

cat("================================================================================\n")
cat("                         INVESTIGATION COMPLETE                                 \n")
cat("================================================================================\n")
