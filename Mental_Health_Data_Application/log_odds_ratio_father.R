#------------------------------------------------------------------------------#
# Log Odds Ratio Estimation: teacher_report for father = 1 vs 0
# Using EBMR with 3 propensity score models
#------------------------------------------------------------------------------#

library(EBMRalgorithmFast)
library(tidyverse)
library(numDeriv)
library(Matrix)

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

n <- nrow(dat)
missing_rate <- 1 - mean(dat$r)

cat("================================================================================\n")
cat("     LOG ODDS RATIO ESTIMATION: teacher_report ~ father                         \n")
cat("================================================================================\n\n")

cat("Sample size:", n, "\n")
cat("Missing rate:", round(missing_rate * 100, 1), "%\n\n")

# Summary statistics
cat("Data summary:\n")
cat("  N(father = 0):", sum(dat$father == 0), "\n")
cat("  N(father = 1):", sum(dat$father == 1), "\n")
cat("  N(R = 1):", sum(dat$r == 1), "\n\n")

#------------------------------------------------------------------------------#
# Complete Case Analysis
#------------------------------------------------------------------------------#
dat_cc <- dat[dat$r == 1, ]
n_cc <- nrow(dat_cc)

# Log odds ratio from complete case
tab_cc <- table(dat_cc$father, dat_cc$teacher_report)
cat("Complete case contingency table:\n")
print(tab_cc)
cat("\n")

# Calculate log odds ratio
a <- tab_cc[2, 2]  # father=1, teacher=1
b <- tab_cc[2, 1]  # father=1, teacher=0
c <- tab_cc[1, 2]  # father=0, teacher=1
d <- tab_cc[1, 1]  # father=0, teacher=0

log_or_cc <- log((a * d) / (b * c))
se_log_or_cc <- sqrt(1/a + 1/b + 1/c + 1/d)

cat("Complete Case Log Odds Ratio:\n")
cat("  Log OR:", round(log_or_cc, 4), "\n")
cat("  SE:", round(se_log_or_cc, 4), "\n")
cat("  95% CI: [", round(log_or_cc - 1.96 * se_log_or_cc, 4), ", ",
    round(log_or_cc + 1.96 * se_log_or_cc, 4), "]\n")
cat("  OR:", round(exp(log_or_cc), 4), "\n\n")

#------------------------------------------------------------------------------#
# EBMR Estimation
#------------------------------------------------------------------------------#
cat("================================================================================\n")
cat("                           EBMR ESTIMATION                                      \n")
cat("================================================================================\n\n")

# Optimal weighting matrix
W <- function(g.matrix) {
  solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
}

# First fit MAR model to get initial values
# For father subgroups, use health and parent_report as covariates
mar_fit <- glm(r ~ health + parent_report, data = dat, family = binomial)
mar_coef <- -coef(mar_fit)  # Negate for inv_link = 1/(1+exp(eta))

cat("MAR Model (for initial values): r ~ health + parent_report\n")
cat("  Coefficients (negated):", round(mar_coef, 4), "\n\n")

# PS specifications: 3 models
# pi_1: r ~ o(y) + health
# pi_2: r ~ o(y) + parent_report
# pi_3: r ~ health + parent_report (MAR)
ps_specifications <- list(
  formula.list = list(
    r ~ o(teacher_report) + health,
    r ~ o(teacher_report) + parent_report,
    r ~ health + parent_report
  ),
  h_x_names.list = list(
    c("health", "parent_report"),
    c("health", "parent_report"),
    c("health", "parent_report")
  ),
  alpha_init.list = list(
    c(mar_coef["(Intercept)"], 0, mar_coef["health"]),
    c(mar_coef["(Intercept)"], 0, mar_coef["parent_report"]),
    c(mar_coef["(Intercept)"], mar_coef["health"], mar_coef["parent_report"])
  ),
  inv_link = function(eta) 1 / (1 + exp(eta))
)

# h_nu function for ensemble step
h_nu <- function(data) {
  cbind(
    health = data$health,
    parent_report = data$parent_report,
    hp = data$health * data$parent_report
  )
}

# Model names for output
model_names <- c("r ~ o(y) + health", "r ~ o(y) + parent_report", "r ~ health + parent_report")
var_names_list <- list(
  c("Intercept", "o(y)", "health"),
  c("Intercept", "o(y)", "parent_report"),
  c("Intercept", "health", "parent_report")
)

#------------------------------------------------------------------------------#
# Estimate for father = 1 subgroup
#------------------------------------------------------------------------------#
cat("--------------------------------------------------------------------------------\n")
cat("  Subgroup: father = 1                                                          \n")
cat("--------------------------------------------------------------------------------\n\n")

dat_f1 <- dat[dat$father == 1, ]
n_f1 <- nrow(dat_f1)
cat("  N:", n_f1, "\n")
cat("  Missing rate:", round((1 - mean(dat_f1$r)) * 100, 1), "%\n\n")

ebmr_f1 <- EBMRAlgorithmFast$new("teacher_report", ps_specifications, dat_f1, W)
result_f1 <- ebmr_f1$EBMR_IPW(h_nu = h_nu, true_ps = NULL)

# Model summaries
for (i in 1:3) {
  ps_fit <- ebmr_f1$ps_fit.list[[i]]
  cat(sprintf("  Model %d: %s\n", i, model_names[i]))

  # Check if GMM fit exists (MNAR models) or GLM (MAR model)
  if (!is.null(ps_fit$gmm_fit) && !is.null(ps_fit$gmm_fit$opt)) {
    cat(sprintf("    Conv: %d, Obj: %.2e\n",
                ps_fit$gmm_fit$opt$convergence, ps_fit$gmm_fit$opt$value))
  } else {
    cat("    Fit: GLM (MAR model)\n")
  }

  coef <- ps_fit$coefficients
  se <- ps_fit$se
  z_val <- coef / se
  p_val <- 2 * pnorm(-abs(z_val))

  cat(sprintf("    %15s %10s %10s %10s %10s\n", "Variable", "Estimate", "SE", "z", "p"))
  for (k in seq_along(coef)) {
    stars <- ""
    if (p_val[k] < 0.001) stars <- "***"
    else if (p_val[k] < 0.01) stars <- "**"
    else if (p_val[k] < 0.05) stars <- "*"
    else if (p_val[k] < 0.1) stars <- "."
    cat(sprintf("    %15s %10.4f %10.4f %10.3f %10.4f%s\n",
                var_names_list[[i]][k], coef[k], se[k], z_val[k], p_val[k], stars))
  }
  cat("\n")
}

cat("  nu.hat:", round(result_f1$nu.hat, 4), "\n")
cat("  w.hat:", round(result_f1$w.hat, 4), "\n")
cat("  mu_f1 (P(y=1|father=1)):", round(result_f1$mu_ipw, 4), "\n")
cat("  SE:", round(result_f1$se_ipw, 4), "\n\n")

#------------------------------------------------------------------------------#
# Estimate for father = 0 subgroup
#------------------------------------------------------------------------------#
cat("--------------------------------------------------------------------------------\n")
cat("  Subgroup: father = 0                                                          \n")
cat("--------------------------------------------------------------------------------\n\n")

dat_f0 <- dat[dat$father == 0, ]
n_f0 <- nrow(dat_f0)
cat("  N:", n_f0, "\n")
cat("  Missing rate:", round((1 - mean(dat_f0$r)) * 100, 1), "%\n\n")

ebmr_f0 <- EBMRAlgorithmFast$new("teacher_report", ps_specifications, dat_f0, W)
result_f0 <- ebmr_f0$EBMR_IPW(h_nu = h_nu, true_ps = NULL)

for (i in 1:3) {
  ps_fit <- ebmr_f0$ps_fit.list[[i]]
  cat(sprintf("  Model %d: %s\n", i, model_names[i]))

  if (!is.null(ps_fit$gmm_fit) && !is.null(ps_fit$gmm_fit$opt)) {
    cat(sprintf("    Conv: %d, Obj: %.2e\n",
                ps_fit$gmm_fit$opt$convergence, ps_fit$gmm_fit$opt$value))
  } else {
    cat("    Fit: GLM (MAR model)\n")
  }

  coef <- ps_fit$coefficients
  se <- ps_fit$se
  z_val <- coef / se
  p_val <- 2 * pnorm(-abs(z_val))

  cat(sprintf("    %15s %10s %10s %10s %10s\n", "Variable", "Estimate", "SE", "z", "p"))
  for (k in seq_along(coef)) {
    stars <- ""
    if (p_val[k] < 0.001) stars <- "***"
    else if (p_val[k] < 0.01) stars <- "**"
    else if (p_val[k] < 0.05) stars <- "*"
    else if (p_val[k] < 0.1) stars <- "."
    cat(sprintf("    %15s %10.4f %10.4f %10.3f %10.4f%s\n",
                var_names_list[[i]][k], coef[k], se[k], z_val[k], p_val[k], stars))
  }
  cat("\n")
}

cat("  nu.hat:", round(result_f0$nu.hat, 4), "\n")
cat("  w.hat:", round(result_f0$w.hat, 4), "\n")
cat("  mu_f0 (P(y=1|father=0)):", round(result_f0$mu_ipw, 4), "\n")
cat("  SE:", round(result_f0$se_ipw, 4), "\n\n")

#------------------------------------------------------------------------------#
# Calculate Log Odds Ratio
#------------------------------------------------------------------------------#
cat("================================================================================\n")
cat("                     LOG ODDS RATIO ESTIMATION RESULT                           \n")
cat("================================================================================\n\n")

mu_1 <- result_f1$mu_ipw
mu_0 <- result_f0$mu_ipw
se_1 <- result_f1$se_ipw
se_0 <- result_f0$se_ipw

# Log odds for each group
log_odds_1 <- log(mu_1 / (1 - mu_1))
log_odds_0 <- log(mu_0 / (1 - mu_0))

# Log odds ratio
log_or_ebmr <- log_odds_1 - log_odds_0

# SE via delta method
# Var(log(p/(1-p))) = Var(p) / (p(1-p))^2
var_log_odds_1 <- se_1^2 / (mu_1 * (1 - mu_1))^2
var_log_odds_0 <- se_0^2 / (mu_0 * (1 - mu_0))^2
se_log_or_ebmr <- sqrt(var_log_odds_1 + var_log_odds_0)

cat("  Probability Estimates:\n")
cat(sprintf("    P(y=1 | father=1) = %.4f (SE: %.4f)\n", mu_1, se_1))
cat(sprintf("    P(y=1 | father=0) = %.4f (SE: %.4f)\n", mu_0, se_0))
cat("\n")

cat("  Log Odds:\n")
cat(sprintf("    log(odds | father=1) = %.4f\n", log_odds_1))
cat(sprintf("    log(odds | father=0) = %.4f\n", log_odds_0))
cat("\n")

cat("  Estimator                  Log OR       SE           95% CI              OR\n")
cat("  ", paste(rep("-", 78), collapse = ""), "\n")

ci_cc <- c(log_or_cc - 1.96 * se_log_or_cc, log_or_cc + 1.96 * se_log_or_cc)
cat(sprintf("  Complete Case          %10.4f   %.4f   [%7.4f, %7.4f]   %7.4f\n",
            log_or_cc, se_log_or_cc, ci_cc[1], ci_cc[2], exp(log_or_cc)))

ci_ebmr <- c(log_or_ebmr - 1.96 * se_log_or_ebmr, log_or_ebmr + 1.96 * se_log_or_ebmr)
cat(sprintf("  EBMR (3 models)        %10.4f   %.4f   [%7.4f, %7.4f]   %7.4f\n",
            log_or_ebmr, se_log_or_ebmr, ci_ebmr[1], ci_ebmr[2], exp(log_or_ebmr)))

cat("  ", paste(rep("-", 78), collapse = ""), "\n")

cat("\n  Difference (EBMR - CC):", round(log_or_ebmr - log_or_cc, 4), "\n")

#------------------------------------------------------------------------------#
# Significance Test
#------------------------------------------------------------------------------#
cat("\n================================================================================\n")
cat("                           SIGNIFICANCE TEST                                    \n")
cat("================================================================================\n\n")

z_cc <- log_or_cc / se_log_or_cc
p_cc <- 2 * pnorm(-abs(z_cc))

z_ebmr <- log_or_ebmr / se_log_or_ebmr
p_ebmr <- 2 * pnorm(-abs(z_ebmr))

cat("  H0: Log OR = 0 (no association between father and teacher_report)\n\n")
cat(sprintf("  Complete Case:   z = %.3f, p = %.4f %s\n",
            z_cc, p_cc, ifelse(p_cc < 0.05, "*", "")))
cat(sprintf("  EBMR:            z = %.3f, p = %.4f %s\n",
            z_ebmr, p_ebmr, ifelse(p_ebmr < 0.05, "*", "")))

#------------------------------------------------------------------------------#
# Convergence Check Summary
#------------------------------------------------------------------------------#
cat("\n================================================================================\n")
cat("                         CONVERGENCE CHECK SUMMARY                              \n")
cat("================================================================================\n\n")

cat("  Subgroup       Model   Conv   Objective        Status\n")
cat("  ", paste(rep("-", 60), collapse = ""), "\n")

for (subgroup in c("f1", "f0")) {
  ebmr_obj <- if (subgroup == "f1") ebmr_f1 else ebmr_f0
  subgroup_label <- if (subgroup == "f1") "father=1" else "father=0"

  for (i in 1:3) {
    ps_fit <- ebmr_obj$ps_fit.list[[i]]

    if (!is.null(ps_fit$gmm_fit) && !is.null(ps_fit$gmm_fit$opt)) {
      conv <- ps_fit$gmm_fit$opt$convergence
      obj <- ps_fit$gmm_fit$opt$value

      status <- "OK"
      if (conv != 0) status <- "FAILED"
      else if (any(abs(ps_fit$coefficients) > 10)) status <- "EXTREME"
      else if (any(ps_fit$se > 10)) status <- "UNSTABLE"
      else if (obj > 0.01) status <- "HIGH OBJ"

      cat(sprintf("  %-12s Model %d    %d    %.2e    %s\n",
                  subgroup_label, i, conv, obj, status))
    } else {
      cat(sprintf("  %-12s Model %d    -    GLM          OK\n",
                  subgroup_label, i))
    }
  }
}

cat("  ", paste(rep("-", 60), collapse = ""), "\n")

cat("\n================================================================================\n")
