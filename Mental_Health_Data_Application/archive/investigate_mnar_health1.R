#------------------------------------------------------------------------------#
# Investigate MNAR Model Coefficients for Health = 1 Group
# Why is the coefficient on y not significant despite prior MNAR evidence?
#------------------------------------------------------------------------------#

library(EBMRalgorithmFast)
library(tidyverse)
library(numDeriv)

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
dat$y <- dat$teacher_report

# Focus on health = 1 group
subdat <- dat[dat$health == 1, ]
subdat$y <- subdat$teacher_report

cat("================================================================================\n")
cat("     INVESTIGATION: MNAR MODEL COEFFICIENTS FOR HEALTH = 1 GROUP               \n")
cat("================================================================================\n\n")

cat("Sample size: n =", nrow(subdat), "\n")
cat("Missing rate:", round(100 * (1 - mean(subdat$r)), 1), "%\n")
cat("Number observed (r=1):", sum(subdat$r), "\n")
cat("Number missing (r=0):", sum(1 - subdat$r), "\n\n")

#------------------------------------------------------------------------------#
# Descriptive statistics
#------------------------------------------------------------------------------#
cat("--------------------------------------------------------------------------------\n")
cat("DESCRIPTIVE STATISTICS\n")
cat("--------------------------------------------------------------------------------\n\n")

cat("Outcome (teacher_report) among observed (r=1):\n")
cat("  Mean:", round(mean(subdat$teacher_report[subdat$r == 1]), 4), "\n")
cat("  SD:  ", round(sd(subdat$teacher_report[subdat$r == 1]), 4), "\n")
cat("  Range:", round(range(subdat$teacher_report[subdat$r == 1]), 4), "\n\n")

cat("Covariates:\n")
cat("  father (observed):      ", round(mean(subdat$father[subdat$r == 1]), 4), "\n")
cat("  father (missing):       ", round(mean(subdat$father[subdat$r == 0]), 4), "\n")
cat("  parent_report (observed):", round(mean(subdat$parent_report[subdat$r == 1]), 4), "\n")
cat("  parent_report (missing): ", round(mean(subdat$parent_report[subdat$r == 0]), 4), "\n\n")

#------------------------------------------------------------------------------#
# Check correlation between y and r among observed
#------------------------------------------------------------------------------#
cat("--------------------------------------------------------------------------------\n")
cat("RELATIONSHIP BETWEEN Y AND MISSINGNESS\n")
cat("--------------------------------------------------------------------------------\n\n")

# Cross-tabulation of y and r
cat("Cross-tabulation of teacher_report (y) and response indicator (r):\n")
print(table(y = subdat$teacher_report, r = subdat$r))
cat("\n")

# Row percentages (missingness rate by y value)
tab <- table(y = subdat$teacher_report, r = subdat$r)
miss_rate_by_y <- tab[, "0"] / (tab[, "0"] + tab[, "1"])
cat("Missing rate by y value:\n")
print(round(miss_rate_by_y, 3))
cat("\n")

# Chi-square test
cat("Chi-square test for independence of y and r:\n")
chisq_test <- chisq.test(table(subdat$teacher_report, subdat$r))
print(chisq_test)
cat("\n")

#------------------------------------------------------------------------------#
# Fit standard logistic regression (treating y as observed for all)
# This is just exploratory - not valid for MNAR!
#------------------------------------------------------------------------------#
cat("--------------------------------------------------------------------------------\n")
cat("EXPLORATORY: STANDARD LOGISTIC REGRESSION (using observed y only)\n")
cat("--------------------------------------------------------------------------------\n\n")

# Among observed cases only
obs_dat <- subdat[subdat$r == 1, ]

cat("Fitting logistic regression for r ~ y + father (observed cases only):\n")
glm_fit1 <- glm(r ~ teacher_report + father, data = subdat, family = binomial())
print(summary(glm_fit1))

cat("\nFitting logistic regression for r ~ y + parent_report (observed cases only):\n")
glm_fit2 <- glm(r ~ teacher_report + parent_report, data = subdat, family = binomial())
print(summary(glm_fit2))

#------------------------------------------------------------------------------#
# Fit MNAR models using EBMRalgorithmFast
#------------------------------------------------------------------------------#
cat("--------------------------------------------------------------------------------\n")
cat("MNAR MODEL FITS (Wang-Shao-Kim 2014)\n")
cat("--------------------------------------------------------------------------------\n\n")

full_ps_specifications <- list(
  formula.list = list(
    r ~ o(teacher_report) + father,
    r ~ o(teacher_report) + parent_report,
    r ~ father + parent_report  # MAR model
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

# Try different initial values
init_sets <- list(
  "default" = list(
    c(0, 1, 0),
    c(0, 1, 0),
    c(0.5, -0.2, -0.5)
  ),
  "from_glm" = list(
    c(coef(glm_fit1)[1], coef(glm_fit1)[2], coef(glm_fit1)[3]),
    c(coef(glm_fit2)[1], coef(glm_fit2)[2], coef(glm_fit2)[3]),
    c(0.5, -0.2, -0.5)
  ),
  "positive_y" = list(
    c(0, 2, 0),
    c(0, 2, 0),
    c(0.5, -0.2, -0.5)
  ),
  "negative_y" = list(
    c(0, -2, 0),
    c(0, -2, 0),
    c(0.5, -0.2, -0.5)
  ),
  "strong_positive" = list(
    c(-1, 5, -0.5),
    c(-1, 5, -0.5),
    c(0.5, -0.2, -0.5)
  ),
  "strong_negative" = list(
    c(1, -5, 0.5),
    c(1, -5, 0.5),
    c(0.5, -0.2, -0.5)
  )
)

model_names <- c("Model 1 (o(y) + father)", "Model 2 (o(y) + parent_report)", "Model 3 (MAR)")

cat("Fitting MNAR models with different initial values:\n\n")

results_summary <- data.frame()

for (init_name in names(init_sets)) {
  cat("=== Initial values:", init_name, "===\n")

  ps_specifications <- list(
    formula.list = full_ps_specifications$formula.list[1:3],
    h_x_names.list = full_ps_specifications$h_x_names.list[1:3],
    alpha_init.list = init_sets[[init_name]],
    inv_link = full_ps_specifications$inv_link
  )

  tryCatch({
    ebmr <- EBMRAlgorithmFast$new("teacher_report", ps_specifications, subdat, W)

    for (j in 1:2) {  # Only MNAR models
      ps_fit <- ebmr$ps_fit.list[[j]]

      alpha <- ps_fit$coefficients
      se <- ps_fit$se
      conv <- ps_fit$gmm_fit$opt$convergence
      obj <- ps_fit$gmm_fit$opt$value

      # Compute z-score and p-value for coefficient on y
      z_y <- alpha[2] / se[2]
      p_y <- 2 * pnorm(-abs(z_y))

      cat(sprintf("  %s:\n", model_names[j]))
      cat(sprintf("    alpha = (%.4f, %.4f, %.4f)\n", alpha[1], alpha[2], alpha[3]))
      cat(sprintf("    SE    = (%.4f, %.4f, %.4f)\n", se[1], se[2], se[3]))
      cat(sprintf("    Coef on y: %.4f (SE: %.4f, z: %.2f, p: %.4f)\n",
                  alpha[2], se[2], z_y, p_y))
      cat(sprintf("    Convergence: %d, Obj: %.2e\n", conv, obj))
      cat(sprintf("    PS range: [%.4f, %.4f]\n\n", min(ps_fit$fitted.values), max(ps_fit$fitted.values)))

      results_summary <- rbind(results_summary, data.frame(
        init = init_name,
        model = j,
        alpha_intercept = alpha[1],
        alpha_y = alpha[2],
        alpha_x = alpha[3],
        se_y = se[2],
        z_y = z_y,
        p_y = p_y,
        conv = conv,
        obj = obj,
        ps_min = min(ps_fit$fitted.values),
        ps_max = max(ps_fit$fitted.values)
      ))
    }
  }, error = function(e) {
    cat("  Error:", e$message, "\n\n")
  })

  cat("\n")
}

#------------------------------------------------------------------------------#
# Summary table
#------------------------------------------------------------------------------#
cat("--------------------------------------------------------------------------------\n")
cat("SUMMARY: COEFFICIENT ON Y ACROSS DIFFERENT INITIAL VALUES\n")
cat("--------------------------------------------------------------------------------\n\n")

cat("Model 1 (o(y) + father):\n")
m1_results <- results_summary[results_summary$model == 1, ]
print(m1_results[, c("init", "alpha_y", "se_y", "z_y", "p_y", "conv", "obj")], row.names = FALSE)

cat("\nModel 2 (o(y) + parent_report):\n")
m2_results <- results_summary[results_summary$model == 2, ]
print(m2_results[, c("init", "alpha_y", "se_y", "z_y", "p_y", "conv", "obj")], row.names = FALSE)

#------------------------------------------------------------------------------#
# Check the moment conditions
#------------------------------------------------------------------------------#
cat("\n--------------------------------------------------------------------------------\n")
cat("CHECKING MOMENT CONDITIONS AT ESTIMATED ALPHA\n")
cat("--------------------------------------------------------------------------------\n\n")

# Re-fit with default init
ps_specifications <- list(
  formula.list = full_ps_specifications$formula.list[1:2],
  h_x_names.list = full_ps_specifications$h_x_names.list[1:2],
  alpha_init.list = init_sets[["default"]][1:2],
  inv_link = full_ps_specifications$inv_link
)

ebmr <- EBMRAlgorithmFast$new("teacher_report", ps_specifications, subdat, W)

for (j in 1:2) {
  ps_fit <- ebmr$ps_fit.list[[j]]
  g_matrix <- ps_fit$gmm_fit$g.matrix

  cat(model_names[j], ":\n")
  cat("  Moment conditions (should be ~0):\n")
  cat("    E[g] =", round(colMeans(g_matrix), 6), "\n")
  cat("    Var[g] diag =", round(diag(var(g_matrix)), 4), "\n\n")
}

#------------------------------------------------------------------------------#
# Compare with health = 0 group
#------------------------------------------------------------------------------#
cat("--------------------------------------------------------------------------------\n")
cat("COMPARISON WITH HEALTH = 0 GROUP\n")
cat("--------------------------------------------------------------------------------\n\n")

subdat0 <- dat[dat$health == 0, ]
subdat0$y <- subdat0$teacher_report

cat("Health = 0 group:\n")
cat("  Sample size: n =", nrow(subdat0), "\n")
cat("  Missing rate:", round(100 * (1 - mean(subdat0$r)), 1), "%\n\n")

# Cross-tabulation
cat("Missing rate by y value (health = 0):\n")
tab0 <- table(y = subdat0$teacher_report, r = subdat0$r)
miss_rate_by_y0 <- tab0[, "0"] / (tab0[, "0"] + tab0[, "1"])
print(round(miss_rate_by_y0, 3))

cat("\nMissing rate by y value (health = 1):\n")
print(round(miss_rate_by_y, 3))

cat("\n\nDifference in missing rate pattern suggests:\n")
cat("  Health = 0: Higher y -> ", ifelse(miss_rate_by_y0[2] > miss_rate_by_y0[1], "MORE", "LESS"), " missing\n")
cat("  Health = 1: Higher y -> ", ifelse(miss_rate_by_y[2] > miss_rate_by_y[1], "MORE", "LESS"), " missing\n")

#------------------------------------------------------------------------------#
# Fit MNAR models for health = 0 group for comparison
#------------------------------------------------------------------------------#
cat("\n--------------------------------------------------------------------------------\n")
cat("MNAR MODEL FIT FOR HEALTH = 0 GROUP (COMPARISON)\n")
cat("--------------------------------------------------------------------------------\n\n")

ps_specifications0 <- list(
  formula.list = full_ps_specifications$formula.list[1:2],
  h_x_names.list = full_ps_specifications$h_x_names.list[1:2],
  alpha_init.list = init_sets[["default"]][1:2],
  inv_link = full_ps_specifications$inv_link
)

ebmr0 <- EBMRAlgorithmFast$new("teacher_report", ps_specifications0, subdat0, W)

for (j in 1:2) {
  ps_fit <- ebmr0$ps_fit.list[[j]]

  alpha <- ps_fit$coefficients
  se <- ps_fit$se
  conv <- ps_fit$gmm_fit$opt$convergence

  z_y <- alpha[2] / se[2]
  p_y <- 2 * pnorm(-abs(z_y))

  cat(sprintf("%s:\n", model_names[j]))
  cat(sprintf("  alpha = (%.4f, %.4f, %.4f)\n", alpha[1], alpha[2], alpha[3]))
  cat(sprintf("  SE    = (%.4f, %.4f, %.4f)\n", se[1], se[2], se[3]))
  cat(sprintf("  Coef on y: %.4f (SE: %.4f, z: %.2f, p: %.4f)\n",
              alpha[2], se[2], z_y, p_y))
  cat(sprintf("  Convergence: %d\n\n", conv))
}

cat("================================================================================\n")
cat("                         INVESTIGATION COMPLETE                                 \n")
cat("================================================================================\n")
