#------------------------------------------------------------------------------#
# Compare alpha1 estimates with different h_x_names specifications
#------------------------------------------------------------------------------#

library(EBMRalgorithmFast)
library(tidyverse)
library(Matrix)

source("MHD_functions.R")

# Read and prepare data
original_data <- read.csv("data_application.csv")
percent <- original_data$Percentage
n <- 2486
class_n <- round(n * percent / 100)
dat <- gen_data(original_data, class_n, n)

W <- function(g.matrix) {
  solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
}

cat("================================================================================\n")
cat("  COMPARING ALPHA1 ESTIMATES WITH DIFFERENT h_x_names\n")
cat("================================================================================\n\n")

cat("Model 1: r ~ o(teacher_report) + health\n")
cat("  Parameters: intercept, o(y), health  (3 parameters)\n")
cat("Initial values: c(-1.0529, 2.2007, -0.3050)\n\n")

#------------------------------------------------------------------------------#
# Case 1: h_x_names = c("health", "father") - just-identified (3 moments, 3 params)
#------------------------------------------------------------------------------#
cat("--------------------------------------------------------------------------------\n")
cat("Case 1: h_x_names = c('health', 'father')\n")
cat("--------------------------------------------------------------------------------\n\n")

ps_spec_1 <- list(
  formula.list = list(r ~ o(teacher_report) + health),
  h_x_names.list = list(c("health", "father")),
  alpha_init.list = list(c(-1.0529, 2.2007, -0.3050)),
  inv_link = function(eta) 1 / (1 + exp(eta))
)

h_nu_1 <- function(data) {
  cbind(health = data$health, father = data$father)
}

ebmr_1 <- EBMRAlgorithmFast$new("teacher_report", ps_spec_1, dat, W)
result_1 <- ebmr_1$EBMR_IPW(h_nu = h_nu_1, true_ps = NULL)

ps_fit_1 <- ebmr_1$ps_fit.list[[1]]
g_matrix_1 <- ps_fit_1$gmm_fit$g_matrix
cat("  Moment conditions structure:\n")
cat("    - E[(R - pi) * 1] = 0       (intercept moment)\n")
cat("    - E[(R - pi) * health] = 0  (health moment)\n")
cat("    - E[(R - pi) * father] = 0  (father moment, from h_x_names)\n")
cat("  Number of moment conditions:", ncol(g_matrix_1), "\n")
cat("  Number of parameters:", length(ps_fit_1$coefficients), "\n")
cat("  Identification:", ifelse(ncol(g_matrix_1) == length(ps_fit_1$coefficients), 
                                 "Just-identified", "Over-identified"), "\n\n")

cat("  Coefficients:\n")
cat("    Intercept:", round(ps_fit_1$coefficients[1], 4), "\n")
cat("    o(y):     ", round(ps_fit_1$coefficients[2], 4), "\n")
cat("    health:   ", round(ps_fit_1$coefficients[3], 4), "\n")
cat("  Convergence:", ps_fit_1$gmm_fit$opt$convergence, "\n")
cat("  Objective:  ", format(ps_fit_1$gmm_fit$opt$value, scientific = TRUE), "\n\n")

# Check moment values at solution
cat("  Sample moments at solution (should be ~0):\n")
g_means_1 <- colMeans(g_matrix_1)
cat("    E[(R-pi)*1]:     ", format(g_means_1[1], scientific = TRUE), "\n")
cat("    E[(R-pi)*health]:", format(g_means_1[2], scientific = TRUE), "\n")
cat("    E[(R-pi)*father]:", format(g_means_1[3], scientific = TRUE), "\n\n")

#------------------------------------------------------------------------------#
# Case 2: h_x_names = c("health", "father", "parent_report") - over-identified
#------------------------------------------------------------------------------#
cat("--------------------------------------------------------------------------------\n")
cat("Case 2: h_x_names = c('health', 'father', 'parent_report')\n")
cat("--------------------------------------------------------------------------------\n\n")

ps_spec_2 <- list(
  formula.list = list(r ~ o(teacher_report) + health),
  h_x_names.list = list(c("health", "father", "parent_report")),
  alpha_init.list = list(c(-1.0529, 2.2007, -0.3050)),
  inv_link = function(eta) 1 / (1 + exp(eta))
)

h_nu_2 <- function(data) {
  cbind(health = data$health, father = data$father, parent_report = data$parent_report)
}

ebmr_2 <- EBMRAlgorithmFast$new("teacher_report", ps_spec_2, dat, W)
result_2 <- ebmr_2$EBMR_IPW(h_nu = h_nu_2, true_ps = NULL)

ps_fit_2 <- ebmr_2$ps_fit.list[[1]]
g_matrix_2 <- ps_fit_2$gmm_fit$g_matrix
cat("  Moment conditions structure:\n")
cat("    - E[(R - pi) * 1] = 0             (intercept moment)\n")
cat("    - E[(R - pi) * health] = 0        (health moment)\n")
cat("    - E[(R - pi) * father] = 0        (father moment)\n")
cat("    - E[(R - pi) * parent_report] = 0 (parent_report moment)\n")
cat("  Number of moment conditions:", ncol(g_matrix_2), "\n")
cat("  Number of parameters:", length(ps_fit_2$coefficients), "\n")
cat("  Identification:", ifelse(ncol(g_matrix_2) == length(ps_fit_2$coefficients), 
                                 "Just-identified", "Over-identified"), "\n\n")

cat("  Coefficients:\n")
cat("    Intercept:", round(ps_fit_2$coefficients[1], 4), "\n")
cat("    o(y):     ", round(ps_fit_2$coefficients[2], 4), "\n")
cat("    health:   ", round(ps_fit_2$coefficients[3], 4), "\n")
cat("  Convergence:", ps_fit_2$gmm_fit$opt$convergence, "\n")
cat("  Objective:  ", format(ps_fit_2$gmm_fit$opt$value, scientific = TRUE), "\n\n")

# Check moment values at solution
cat("  Sample moments at solution:\n")
g_means_2 <- colMeans(g_matrix_2)
cat("    E[(R-pi)*1]:            ", format(g_means_2[1], scientific = TRUE), "\n")
cat("    E[(R-pi)*health]:       ", format(g_means_2[2], scientific = TRUE), "\n")
cat("    E[(R-pi)*father]:       ", format(g_means_2[3], scientific = TRUE), "\n")
cat("    E[(R-pi)*parent_report]:", format(g_means_2[4], scientific = TRUE), "\n\n")

#------------------------------------------------------------------------------#
# Summary
#------------------------------------------------------------------------------#
cat("================================================================================\n")
cat("  EXPLANATION\n")
cat("================================================================================\n\n")

cat("The GMM solves: min_alpha || E[g(alpha)] ||_W^2\n\n")

cat("Case 1 (Just-identified):\n")
cat("  - 3 parameters, 3 moment conditions\n")
cat("  - Unique solution exists that makes ALL moments exactly zero\n")
cat("  - Objective value ~0 (1.57e-11)\n")
cat("  - The 'father' moment E[(R-pi)*father]=0 constrains the intercept/o(y) estimates\n\n")

cat("Case 2 (Over-identified):\n")
cat("  - 3 parameters, 4 moment conditions\n")
cat("  - CANNOT make all 4 moments exactly zero simultaneously\n")
cat("  - GMM finds compromise that minimizes weighted sum of squared moments\n")
cat("  - Objective value = 2.7e-03 (not zero!)\n")
cat("  - Different coefficients because it's balancing 4 conditions with 3 params\n\n")

cat("Key insight:\n")
cat("  The additional 'parent_report' moment condition introduces constraints that\n")
cat("  the just-identified solution doesn't satisfy. This forces GMM to find a\n")
cat("  different compromise solution.\n\n")

# Check if Case 1 solution satisfies Case 2 moments
cat("--------------------------------------------------------------------------------\n")
cat("  Does Case 1 solution satisfy the parent_report moment?\n")
cat("--------------------------------------------------------------------------------\n\n")

# Compute propensity scores using Case 1 coefficients
alpha_1 <- ps_fit_1$coefficients
eta_1 <- alpha_1[1] + alpha_1[2] * dat$teacher_report + alpha_1[3] * dat$health
pi_1 <- 1 / (1 + exp(eta_1))

# Check parent_report moment with Case 1 solution
parent_moment_case1 <- mean((dat$r - pi_1) * dat$parent_report)
cat("  E[(R - pi_case1) * parent_report] =", format(parent_moment_case1, scientific = TRUE), "\n")
cat("  This is NOT zero, explaining why Case 2 (which includes this moment)\n")
cat("  must find a different solution.\n\n")

# Check all moments with Case 1 solution
cat("  All moments evaluated at Case 1 solution:\n")
cat("    E[(R-pi)*1]:            ", format(mean(dat$r - pi_1), scientific = TRUE), "\n")
cat("    E[(R-pi)*health]:       ", format(mean((dat$r - pi_1) * dat$health), scientific = TRUE), "\n")
cat("    E[(R-pi)*father]:       ", format(mean((dat$r - pi_1) * dat$father), scientific = TRUE), "\n")
cat("    E[(R-pi)*parent_report]:", format(parent_moment_case1, scientific = TRUE), " <-- NOT ZERO\n")

cat("\n================================================================================\n")
