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

cat("  Coefficients:\n")
cat("    Intercept:", round(ps_fit_1$coefficients[1], 4), "\n")
cat("    o(y):     ", round(ps_fit_1$coefficients[2], 4), "\n")
cat("    health:   ", round(ps_fit_1$coefficients[3], 4), "\n")
cat("  Convergence:", ps_fit_1$gmm_fit$opt$convergence, "\n")
cat("  Objective:  ", format(ps_fit_1$gmm_fit$opt$value, scientific = TRUE), "\n\n")

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

cat("  Coefficients:\n")
cat("    Intercept:", round(ps_fit_2$coefficients[1], 4), "\n")
cat("    o(y):     ", round(ps_fit_2$coefficients[2], 4), "\n")
cat("    health:   ", round(ps_fit_2$coefficients[3], 4), "\n")
cat("  Convergence:", ps_fit_2$gmm_fit$opt$convergence, "\n")
cat("  Objective:  ", format(ps_fit_2$gmm_fit$opt$value, scientific = TRUE), "\n\n")

#------------------------------------------------------------------------------#
# Verify the difference by computing moments manually
#------------------------------------------------------------------------------#
cat("================================================================================\n")
cat("  VERIFYING THE DIFFERENCE\n")
cat("================================================================================\n\n")

# Case 1 solution
alpha_1 <- ps_fit_1$coefficients
eta_1 <- alpha_1[1] + alpha_1[2] * dat$teacher_report + alpha_1[3] * dat$health
pi_1 <- 1 / (1 + exp(eta_1))

cat("Case 1 solution moments (h_x = health, father):\n")
cat("  E[(R-pi)*1]:            ", format(mean(dat$r - pi_1), scientific = TRUE), "\n")
cat("  E[(R-pi)*health]:       ", format(mean((dat$r - pi_1) * dat$health), scientific = TRUE), "\n")
cat("  E[(R-pi)*father]:       ", format(mean((dat$r - pi_1) * dat$father), scientific = TRUE), "\n")
cat("  E[(R-pi)*parent_report]:", format(mean((dat$r - pi_1) * dat$parent_report), scientific = TRUE), " <-- NOT included in Case 1\n\n")

# Case 2 solution
alpha_2 <- ps_fit_2$coefficients
eta_2 <- alpha_2[1] + alpha_2[2] * dat$teacher_report + alpha_2[3] * dat$health
pi_2 <- 1 / (1 + exp(eta_2))

cat("Case 2 solution moments (h_x = health, father, parent_report):\n")
cat("  E[(R-pi)*1]:            ", format(mean(dat$r - pi_2), scientific = TRUE), "\n")
cat("  E[(R-pi)*health]:       ", format(mean((dat$r - pi_2) * dat$health), scientific = TRUE), "\n")
cat("  E[(R-pi)*father]:       ", format(mean((dat$r - pi_2) * dat$father), scientific = TRUE), "\n")
cat("  E[(R-pi)*parent_report]:", format(mean((dat$r - pi_2) * dat$parent_report), scientific = TRUE), "\n\n")

#------------------------------------------------------------------------------#
# Summary comparison table
#------------------------------------------------------------------------------#
cat("================================================================================\n")
cat("  SUMMARY COMPARISON\n")
cat("================================================================================\n\n")

cat(sprintf("  %-35s %12s %12s %12s %12s\n", "h_x_names", "Intercept", "o(y)", "health", "Objective"))
cat("  ", paste(rep("-", 85), collapse = ""), "\n")
cat(sprintf("  %-35s %12.4f %12.4f %12.4f %12.2e\n", 
            "c('health', 'father')", 
            ps_fit_1$coefficients[1], ps_fit_1$coefficients[2], ps_fit_1$coefficients[3],
            ps_fit_1$gmm_fit$opt$value))
cat(sprintf("  %-35s %12.4f %12.4f %12.4f %12.2e\n", 
            "c('health', 'father', 'parent_report')", 
            ps_fit_2$coefficients[1], ps_fit_2$coefficients[2], ps_fit_2$coefficients[3],
            ps_fit_2$gmm_fit$opt$value))
cat("  ", paste(rep("-", 85), collapse = ""), "\n\n")

cat("================================================================================\n")
cat("  EXPLANATION\n")
cat("================================================================================\n\n")

cat("The Wang-Shao-Kim (2014) GMM method uses moment conditions:\n")
cat("  E[h(X) * (R - pi(X,Y; alpha))] = 0\n\n")

cat("For Model 1: r ~ o(y) + health, the PS is:\n")
cat("  pi(X,Y) = 1 / (1 + exp(alpha0 + alpha1*o(y) + alpha2*health))\n\n")

cat("Case 1: h_x_names = c('health', 'father')\n")
cat("  - Uses h(X) = (1, health, father)  -->  3 moment conditions\n")
cat("  - 3 parameters to estimate  -->  JUST-IDENTIFIED\n")
cat("  - Unique solution makes all moments exactly zero\n")
cat("  - Objective = 1.57e-11 (essentially 0)\n\n")

cat("Case 2: h_x_names = c('health', 'father', 'parent_report')\n")
cat("  - Uses h(X) = (1, health, father, parent_report)  -->  4 moment conditions\n")
cat("  - 3 parameters to estimate  -->  OVER-IDENTIFIED\n")
cat("  - Cannot make all 4 moments exactly zero with only 3 parameters\n")
cat("  - GMM finds weighted least squares compromise\n")
cat("  - Objective = 2.71e-03 (not zero - there's remaining moment imbalance)\n\n")

cat("Why the coefficients differ:\n")
cat("  The Case 1 solution satisfies E[(R-pi)*parent_report] = 0.0163, which is NOT zero.\n")
cat("  When we add parent_report as a moment condition (Case 2), GMM must adjust\n")
cat("  the coefficients to reduce this imbalance, but this worsens the fit to the\n")
cat("  other moments. The result is a compromise solution with different coefficients.\n\n")

cat("Key insight:\n")
cat("  Adding more h_x variables provides MORE information for estimation but\n")
cat("  leads to DIFFERENT estimates if the model is not correctly specified.\n")
cat("  If the true PS model is correct, over-identification should give similar\n")
cat("  but more efficient estimates. Large differences suggest model misspecification.\n")

cat("\n================================================================================\n")
