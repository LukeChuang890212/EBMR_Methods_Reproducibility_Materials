library(EBMRalgorithmFast)
library(tidyverse)
library(numDeriv)
library(Matrix)

# Source utility functions
source("MHD_functions.R")

# Read and prepare data
original_data <- read.csv("data_application.csv")
percent <- original_data$Percentage
n <- 2486
class_n <- round(n * percent / 100)

dat <- gen_data(original_data, class_n, n)
dat$y <- dat$teacher_report

# Complete Case
n <- nrow(dat)
missing_rate <- 1 - mean(dat$r)
mu_cc <- mean(dat$teacher_report[dat$r == 1])
se_cc <- sd(dat$teacher_report[dat$r == 1]) / sqrt(sum(dat$r))

cat("================================================================================\n")
cat("       POPULATION MEAN ESTIMATION (Updated h_nu with 4 terms)                   \n")
cat("================================================================================\n\n")

cat("Sample size:", n, "\n")
cat("Missing rate:", round(missing_rate * 100, 1), "%\n")
cat("Complete case mean:", round(mu_cc, 4), "\n\n")

# PS specifications
full_ps_specifications <- list(
  formula.list = list(
    r ~ o(teacher_report) + health + father,
    r ~ o(teacher_report) + health + parent_report,
    r ~ o(teacher_report) + father + parent_report
  ),
  h_x_names.list = list(
    c("health", "father", "parent_report"),
    c("health", "father", "parent_report"),
    c("health", "father", "parent_report")
  ),
  inv_link = function(eta) 1 / (1 + exp(eta))
)

W <- function(g.matrix) {
  solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
}

alpha_init.list <- list(c(0, 1, 0, 0), c(0, 1, 0, 0), c(0, 1, 0, 0))

# UPDATED h_nu function (4 terms instead of 6)
h_nu <- function(data) cbind(
  health = data$health,
  father = data$father,
  parent_report = data$parent_report,
  hf = data$health * data$father
)

# Model combinations
all_model_sets <- list(1, 2, 3, c(1,2), c(1,3), c(2,3), c(1,2,3))
model_set_labels <- c("100", "010", "001", "110", "101", "011", "111")

results_dir <- "MHD_results"
results_list <- list()

cat("=== Running EBMR estimation ===\n\n")

for (j in 1:length(all_model_sets)) {
  model_set <- all_model_sets[[j]]
  label <- model_set_labels[j]
  cat("Processing:", label, "\n")

  ps_specifications <- list(
    formula.list = full_ps_specifications$formula.list[model_set],
    h_x_names.list = full_ps_specifications$h_x_names.list[model_set],
    alpha_init.list = alpha_init.list[model_set],
    inv_link = full_ps_specifications$inv_link
  )

  ebmr <- EBMRAlgorithmFast$new("teacher_report", ps_specifications, dat, W)
  result <- ebmr$EBMR_IPW(h_nu = h_nu, true_ps = NULL)

  saveRDS(result, paste0(results_dir, "/popmean_EBMR_IPW_", label, ".RDS"))
  results_list[[label]] <- result

  if (j == 7) {
    ps_fit.list_full <- ebmr$ps_fit.list
    nu_hat_full <- result$nu.hat
    w_hat_full <- result$w.hat
  }
}

# Results table
cat("\n================================================================================\n")
cat("                         ESTIMATION RESULTS                                     \n")
cat("================================================================================\n\n")

cat(sprintf("%-15s %-20s %10s %10s\n", "Estimator", "Models", "PE", "SE"))
cat(paste(rep("-", 60), collapse = ""), "\n")

for (j in 1:length(all_model_sets)) {
  label <- model_set_labels[j]
  result <- results_list[[label]]
  models <- paste("pi", all_model_sets[[j]], sep = "", collapse = ",")
  cat(sprintf("%-15s %-20s %10.4f %10.4f\n", paste0("mu_", label), models, result$mu_ipw, result$se_ipw))
}
cat(sprintf("%-15s %-20s %10.4f %10.4f\n", "Complete Case", "N/A", mu_cc, se_cc))
cat(paste(rep("-", 60), collapse = ""), "\n")

# Nu and W
cat("\n================================================================================\n")
cat("                     NU AND W FOR EACH MODEL COMBINATION                        \n")
cat("================================================================================\n\n")

cat(sprintf("%-12s %30s %30s\n", "Estimator", "nu.hat", "w.hat"))
cat(paste(rep("-", 75), collapse = ""), "\n")

for (j in 1:length(all_model_sets)) {
  label <- model_set_labels[j]
  result <- results_list[[label]]
  nu_str <- paste(sprintf("%.4f", result$nu.hat), collapse = ", ")
  w_str <- paste(sprintf("%.4f", result$w.hat), collapse = ", ")
  cat(sprintf("%-12s %30s %30s\n", paste0("mu_", label), nu_str, w_str))
}
cat(paste(rep("-", 75), collapse = ""), "\n")

# Key findings
cat("\n================================================================================\n")
cat("                           KEY FINDINGS                                         \n")
cat("================================================================================\n\n")

full_result <- results_list[["111"]]
cat(sprintf("Complete Case:    %.4f (SE: %.4f)\n", mu_cc, se_cc))
cat(sprintf("EBMR mu_111:      %.4f (SE: %.4f)\n", full_result$mu_ipw, full_result$se_ipw))
cat(sprintf("Difference:       %.4f\n", full_result$mu_ipw - mu_cc))

cat("\nModel Weights (full ensemble mu_111):\n")
for (i in 1:3) {
  cat(sprintf("  pi%d: %.4f (%.1f%%)\n", i, w_hat_full[i], w_hat_full[i] * 100))
}
cat("\n")
