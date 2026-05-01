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
cat("Initial values: c(-1.0529, 2.2007, -0.3050)\n\n")

#------------------------------------------------------------------------------#
# Case 1: h_x_names = c("health", "father")
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
cat("  Number of moment conditions:", ncol(ps_fit_1$gmm_fit$g_matrix), "\n")
cat("  Number of parameters:", length(ps_fit_1$coefficients), "\n")
cat("  Identification:", ifelse(ncol(ps_fit_1$gmm_fit$g_matrix) == length(ps_fit_1$coefficients), 
                                 "Just-identified", "Over-identified"), "\n\n")

cat("  Coefficients:\n")
cat("    Intercept:", round(ps_fit_1$coefficients[1], 4), "\n")
cat("    o(y):     ", round(ps_fit_1$coefficients[2], 4), "\n")
cat("    health:   ", round(ps_fit_1$coefficients[3], 4), "\n")
cat("  Convergence:", ps_fit_1$gmm_fit$opt$convergence, "\n")
cat("  Objective:  ", format(ps_fit_1$gmm_fit$opt$value, scientific = TRUE), "\n\n")

#------------------------------------------------------------------------------#
# Case 2: h_x_names = c("health", "father", "parent_report")
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
cat("  Number of moment conditions:", ncol(ps_fit_2$gmm_fit$g_matrix), "\n")
cat("  Number of parameters:", length(ps_fit_2$coefficients), "\n")
cat("  Identification:", ifelse(ncol(ps_fit_2$gmm_fit$g_matrix) == length(ps_fit_2$coefficients), 
                                 "Just-identified", "Over-identified"), "\n\n")

cat("  Coefficients:\n")
cat("    Intercept:", round(ps_fit_2$coefficients[1], 4), "\n")
cat("    o(y):     ", round(ps_fit_2$coefficients[2], 4), "\n")
cat("    health:   ", round(ps_fit_2$coefficients[3], 4), "\n")
cat("  Convergence:", ps_fit_2$gmm_fit$opt$convergence, "\n")
cat("  Objective:  ", format(ps_fit_2$gmm_fit$opt$value, scientific = TRUE), "\n\n")

#------------------------------------------------------------------------------#
# Case 3: h_x_names = c("health") - just-identified
#------------------------------------------------------------------------------#
cat("--------------------------------------------------------------------------------\n")
cat("Case 3: h_x_names = c('health') [Just-identified]\n")
cat("--------------------------------------------------------------------------------\n\n")

ps_spec_3 <- list(
  formula.list = list(r ~ o(teacher_report) + health),
  h_x_names.list = list(c("health")),
  alpha_init.list = list(c(-1.0529, 2.2007, -0.3050)),
  inv_link = function(eta) 1 / (1 + exp(eta))
)

h_nu_3 <- function(data) {
  cbind(health = data$health)
}

ebmr_3 <- EBMRAlgorithmFast$new("teacher_report", ps_spec_3, dat, W)
result_3 <- ebmr_3$EBMR_IPW(h_nu = h_nu_3, true_ps = NULL)

ps_fit_3 <- ebmr_3$ps_fit.list[[1]]
cat("  Number of moment conditions:", ncol(ps_fit_3$gmm_fit$g_matrix), "\n")
cat("  Number of parameters:", length(ps_fit_3$coefficients), "\n")
cat("  Identification:", ifelse(ncol(ps_fit_3$gmm_fit$g_matrix) == length(ps_fit_3$coefficients), 
                                 "Just-identified", "Over-identified"), "\n\n")

cat("  Coefficients:\n")
cat("    Intercept:", round(ps_fit_3$coefficients[1], 4), "\n")
cat("    o(y):     ", round(ps_fit_3$coefficients[2], 4), "\n")
cat("    health:   ", round(ps_fit_3$coefficients[3], 4), "\n")
cat("  Convergence:", ps_fit_3$gmm_fit$opt$convergence, "\n")
cat("  Objective:  ", format(ps_fit_3$gmm_fit$opt$value, scientific = TRUE), "\n\n")

#------------------------------------------------------------------------------#
# Summary comparison
#------------------------------------------------------------------------------#
cat("================================================================================\n")
cat("  SUMMARY COMPARISON\n")
cat("================================================================================\n\n")

cat(sprintf("  %-30s %12s %12s %12s\n", "h_x_names", "Intercept", "o(y)", "health"))
cat("  ", paste(rep("-", 70), collapse = ""), "\n")
cat(sprintf("  %-30s %12.4f %12.4f %12.4f\n", 
            "c('health')", 
            ps_fit_3$coefficients[1], ps_fit_3$coefficients[2], ps_fit_3$coefficients[3]))
cat(sprintf("  %-30s %12.4f %12.4f %12.4f\n", 
            "c('health', 'father')", 
            ps_fit_1$coefficients[1], ps_fit_1$coefficients[2], ps_fit_1$coefficients[3]))
cat(sprintf("  %-30s %12.4f %12.4f %12.4f\n", 
            "c('health', 'father', 'parent_report')", 
            ps_fit_2$coefficients[1], ps_fit_2$coefficients[2], ps_fit_2$coefficients[3]))
cat("  ", paste(rep("-", 70), collapse = ""), "\n\n")

cat("  The h_x_names parameter determines the auxiliary variables h(X) used in\n")
cat("  the moment conditions: E[h(X) * (R - pi(X,Y))] = 0\n\n")

cat("  More h_x variables means more moment conditions, leading to over-identification.\n")
cat("  With different moment conditions, the GMM minimizes different objective functions,\n")
cat("  potentially leading to different coefficient estimates.\n")

cat("\n================================================================================\n")
