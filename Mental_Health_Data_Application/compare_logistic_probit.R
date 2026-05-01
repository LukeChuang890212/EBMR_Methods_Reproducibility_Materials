#------------------------------------------------------------------------------#
# Compare Logistic Complement vs Probit Link for Population Mean Estimation
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
mu_cc <- mean(dat$teacher_report[dat$r == 1])
se_cc <- sd(dat$teacher_report[dat$r == 1]) / sqrt(sum(dat$r))

cat("================================================================================\n")
cat("     COMPARISON: LOGISTIC COMPLEMENT vs PROBIT LINK                            \n")
cat("================================================================================\n\n")

cat("Sample size:", n, "\n")
cat("Missing rate:", round(missing_rate * 100, 1), "%\n")
cat("Complete case mean:", round(mu_cc, 4), "\n\n")

#------------------------------------------------------------------------------#
# Common Setup
#------------------------------------------------------------------------------#
W <- function(g.matrix) {
  solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
}

h_nu <- function(data) {
  cbind(
    health = data$health,
    father = data$father,
    hf = data$health * data$father
  )
}

model_names <- c("r ~ o(y) + health", "r ~ o(y) + father", "r ~ health + father")
var_names_list <- list(
  c("Intercept", "o(y)", "health"),
  c("Intercept", "o(y)", "father"),
  c("Intercept", "health", "father")
)

#------------------------------------------------------------------------------#
# LOGISTIC COMPLEMENT LINK: pi = 1/(1+exp(eta))
#------------------------------------------------------------------------------#
cat("================================================================================\n")
cat("                     LOGISTIC COMPLEMENT LINK                                   \n")
cat("                     pi = 1/(1+exp(eta))                                        \n")
cat("================================================================================\n\n")

ps_spec_logistic <- list(
  formula.list = list(
    r ~ o(teacher_report) + health,
    r ~ o(teacher_report) + father,
    r ~ health + father
  ),
  h_x_names.list = list(
    c("health", "father", "parent_report"),
    c("health", "father", "parent_report"),
    c("health", "father", "parent_report")
  ),
  alpha_init.list = list(
    c(-1.0529, 2.2007, -0.3050),
    NULL,
    NULL
  ),
  inv_link = function(eta) 1 / (1 + exp(eta))
)

ebmr_logistic <- EBMRAlgorithmFast$new("teacher_report", ps_spec_logistic, dat, W)
result_logistic <- ebmr_logistic$EBMR_IPW(h_nu = h_nu, true_ps = NULL)

cat("Model Coefficients:\n\n")
for (i in 1:3) {
  ps_fit <- ebmr_logistic$ps_fit.list[[i]]
  cat(sprintf("  Model %d: %s\n", i, model_names[i]))
  cat(sprintf("    Link type: %s\n", ps_fit$link_type))
  cat(sprintf("    Coefficients: %s\n", paste(round(ps_fit$coefficients, 4), collapse = ", ")))
  cat(sprintf("    Conv: %d, Obj: %.2e\n\n",
              ps_fit$gmm_fit$opt$convergence, ps_fit$gmm_fit$opt$value))
}

cat("Ensemble Results:\n")
cat("  nu.hat:", round(result_logistic$nu.hat, 4), "\n")
cat("  w.hat:", round(result_logistic$w.hat, 4), "\n")
cat("  mu_ipw:", round(result_logistic$mu_ipw, 4), "\n")
cat("  se_ipw:", round(result_logistic$se_ipw, 4), "\n\n")

#------------------------------------------------------------------------------#
# PROBIT LINK: pi = pnorm(eta)
#------------------------------------------------------------------------------#
cat("================================================================================\n")
cat("                     PROBIT LINK                                                \n")
cat("                     pi = pnorm(eta)                                            \n")
cat("================================================================================\n\n")

# For probit, we need different initial values since the scale is different
# probit coefficients are typically smaller than logit by factor of ~1.6
# Also use probit_complement: pi = pnorm(-eta) = 1 - pnorm(eta) to match logistic_complement convention
ps_spec_probit <- list(
  formula.list = list(
    r ~ o(teacher_report) + health,
    r ~ o(teacher_report) + father,
    r ~ health + father
  ),
  h_x_names.list = list(
    c("health", "father", "parent_report"),
    c("health", "father", "parent_report"),
    c("health", "father", "parent_report")
  ),
  alpha_init.list = list(
    c(-0.15, 0.1, -0.12),  # Smaller values for probit scale
    NULL,
    NULL
  ),
  inv_link = function(eta) pnorm(-eta)  # probit_complement to match logistic_complement
)

# Check that link is detected correctly
cat("Testing probit_complement link detection:\n")
test_link <- function(eta) pnorm(-eta)
cat("  test_link(1) =", test_link(1), "\n")
cat("  test_link(0) =", test_link(0), "\n")
cat("  Expected: ~0.159 at eta=1, 0.5 at eta=0\n\n")

ebmr_probit <- EBMRAlgorithmFast$new("teacher_report", ps_spec_probit, dat, W)
result_probit <- ebmr_probit$EBMR_IPW(h_nu = h_nu, true_ps = NULL)

cat("Model Coefficients:\n\n")
for (i in 1:3) {
  ps_fit <- ebmr_probit$ps_fit.list[[i]]
  cat(sprintf("  Model %d: %s\n", i, model_names[i]))
  cat(sprintf("    Link type: %s\n", ps_fit$link_type))
  cat(sprintf("    Coefficients: %s\n", paste(round(ps_fit$coefficients, 4), collapse = ", ")))
  cat(sprintf("    Conv: %d, Obj: %.2e\n\n",
              ps_fit$gmm_fit$opt$convergence, ps_fit$gmm_fit$opt$value))
}

cat("Ensemble Results:\n")
cat("  nu.hat:", round(result_probit$nu.hat, 4), "\n")
cat("  w.hat:", round(result_probit$w.hat, 4), "\n")
cat("  mu_ipw:", round(result_probit$mu_ipw, 4), "\n")
cat("  se_ipw:", round(result_probit$se_ipw, 4), "\n\n")

#------------------------------------------------------------------------------#
# COMPARISON SUMMARY
#------------------------------------------------------------------------------#
cat("================================================================================\n")
cat("                     COMPARISON SUMMARY                                         \n")
cat("================================================================================\n\n")

cat(sprintf("  %-25s %12s %12s %12s\n", "Estimator", "Estimate", "SE", "95% CI"))
cat("  ", paste(rep("-", 65), collapse = ""), "\n")

ci_cc <- c(mu_cc - 1.96 * se_cc, mu_cc + 1.96 * se_cc)
cat(sprintf("  %-25s %12.4f %12.4f     [%.4f, %.4f]\n",
            "Complete Case", mu_cc, se_cc, ci_cc[1], ci_cc[2]))

ci_logistic <- c(result_logistic$mu_ipw - 1.96 * result_logistic$se_ipw,
                 result_logistic$mu_ipw + 1.96 * result_logistic$se_ipw)
cat(sprintf("  %-25s %12.4f %12.4f     [%.4f, %.4f]\n",
            "EBMR (Logistic Compl.)", result_logistic$mu_ipw, result_logistic$se_ipw,
            ci_logistic[1], ci_logistic[2]))

ci_probit <- c(result_probit$mu_ipw - 1.96 * result_probit$se_ipw,
               result_probit$mu_ipw + 1.96 * result_probit$se_ipw)
cat(sprintf("  %-25s %12.4f %12.4f     [%.4f, %.4f]\n",
            "EBMR (Probit)", result_probit$mu_ipw, result_probit$se_ipw,
            ci_probit[1], ci_probit[2]))

cat("  ", paste(rep("-", 65), collapse = ""), "\n")

cat("\n  Model Weights Comparison:\n\n")
cat(sprintf("  %-10s %15s %15s\n", "Model", "Logistic Compl.", "Probit"))
cat("  ", paste(rep("-", 45), collapse = ""), "\n")
for (i in 1:3) {
  cat(sprintf("  %-10s %15.4f %15.4f\n",
              paste0("Model ", i), result_logistic$w.hat[i], result_probit$w.hat[i]))
}

cat("\n================================================================================\n")
