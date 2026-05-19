#------------------------------------------------------------------------------#
# Mental Health Data Analysis using EBMRalgorithmFast2
# This script performs the essential analysis from Mental_Health_Data_Application2.Rmd
# using the EBMRalgorithmFast2 package
#------------------------------------------------------------------------------#

library(EBMRalgorithmFast2)
library(tidyverse)
library(numDeriv)
library(Matrix)
library(nnet)
library(foreach)
library(doParallel)
library(doSNOW)

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

#------------------------------------------------------------------------------#
# Primary analysis
#------------------------------------------------------------------------------#
n <- nrow(dat)
missing_rate <- 1 - mean(dat$r)
mu_cc <- mean(dat$teacher_report[dat$r == 1])
se_cc <- sd(dat$teacher_report[dat$r == 1]) / sqrt(sum(dat$r))

cat("Sample size:", n, "\n")
cat("Missing rate:", round(missing_rate * 100, 1), "%\n")
cat("Complete case mean:", round(mu_cc, 4), "\n")
cat("Complete case SE:", round(se_cc, 4), "\n\n")

#------------------------------------------------------------------------------#
# Basic Setup for EBMR analysis
#------------------------------------------------------------------------------#
# EBMRalgorithmFast2 uses explicit formulas instead of o() function
full_ps_specifications <- list(
  formula.list = list(
    r ~ teacher_report + father,
    r ~ teacher_report + parent_report,
    r ~ father + parent_report # MAR model (no teacher_report)
  ),
  h_x_names.list = list(
    c("father", "parent_report"),
    c("father", "parent_report"),
    c("father", "parent_report")
  ),
  outcome = "teacher_report",  # Required for EBMRalgorithmFast2
  inv_link = function(eta) 1 / (1 + exp(eta))
)

# Optimal weighting matrix
W <- function(g.matrix) {
  solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
}

health <- 0:1
subset_names <- c("health0", "health1")

# Effect function (log odds ratio) and delta method are defined in MHD_functions.R

# Initial values for optimization
# Updated for EBMRalgorithmFast2 formula structure
init.list <- list(
  alpha_init.list0 = list(
    c(-0.3834341, 1.07945126, -0.1318883, 0),  # Model 1: intercept, teacher_report, father, teacher_report:father
    c(-4.2670736, 5.9582123, -0.5972129, 0),   # Model 2: intercept, teacher_report, parent_report, teacher_report:parent_report
    c(0.5, -0.2, -0.5)  # MAR model: intercept, father, parent_report
  ),
  alpha_init.list1 = list(
    c(0.09019492, -1.511548, -0.3946616, 0),   # Model 1
    c(-1.969056, 3.5823611, -1.37232734, 0),   # Model 2
    c(0.5, -0.2, -0.5)  # MAR model: intercept, father, parent_report
  )
)

# Helper functions (summary.gmm, rm.extreme) are defined in MHD_functions.R

# Add interaction terms to data for h_nu and h_x functions
dat$fp <- dat$father * dat$parent_report

#------------------------------------------------------------------------------#
# Estimate p_0 and p_1 using all combinations of the candidate models
# Now allowing different model combinations for each health group
#------------------------------------------------------------------------------#

# Generate all possible model sets (7 total: 1, 2, 3, 12, 13, 23, 123)
all_model_sets <- list()
for (model_num in 1:3) {
  model_combinations <- combn(3, model_num)
  for (i in 1:ncol(model_combinations)) {
    all_model_sets <- c(all_model_sets, list(model_combinations[, i]))
  }
}
# all_model_sets: [[1]] = 1, [[2]] = 2, [[3]] = 3, [[4]] = c(1,2), [[5]] = c(1,3), [[6]] = c(2,3), [[7]] = c(1,2,3)

# Create model set labels for naming
model_set_labels <- sapply(all_model_sets, function(x) paste(x, collapse = ""))
# "1", "2", "3", "12", "13", "23", "123"

# Check if all EBMR results already exist (for each health group separately)
all_ebmr_files_exist <- TRUE
for (k in 1:length(health)) {
  for (model_set in all_model_sets) {
    save_file <- paste0("MHD_results/EBMR_IPW_", subset_names[k], "_",
                        paste(model_set, collapse = ""), "_W_optimal_OR_test10.RDS")
    if (!file.exists(save_file)) {
      all_ebmr_files_exist <- FALSE
      break
    }
  }
  if (!all_ebmr_files_exist) break
}

ps_fit.list0 <- list()
ps_fit.list1 <- list()
nu_hat0 <- NULL
nu_hat1 <- NULL
w_hat0 <- NULL
w_hat1 <- NULL

if (all_ebmr_files_exist) {
  cat("=== EBMR results already exist, loading from files ===\n\n")

  # Load the full model results to get nu_hat and w_hat
  result0 <- readRDS("MHD_results/EBMR_IPW_health0_123_W_optimal_OR_test10.RDS")
  result1 <- readRDS("MHD_results/EBMR_IPW_health1_123_W_optimal_OR_test10.RDS")
  nu_hat0 <- result0$nu.hat
  nu_hat1 <- result1$nu.hat
  w_hat0 <- result0$w.hat
  w_hat1 <- result1$w.hat

  # Need to refit to get ps_fit.list for summary (or skip summary)
  for (k in 1:length(health)) {
    subdat <- dat[dat$health == health[k], ]
    subdat$y <- subdat$teacher_report

    ps_specifications <- list(
      formula.list = full_ps_specifications$formula.list[1:3],
      h_x_names.list = full_ps_specifications$h_x_names.list[1:3],
      alpha_init.list = init.list[[k]][1:3],
      outcome = full_ps_specifications$outcome,
      inv_link = full_ps_specifications$inv_link
    )

    ebmr <- EBMRAlgorithmFast2$new("teacher_report", ps_specifications, subdat, W)
    if (k == 1) {
      ps_fit.list0 <- ebmr$ps_fit.list
    } else {
      ps_fit.list1 <- ebmr$ps_fit.list
    }
  }
} else {
  cat("=== Running EBMR estimation for all model combinations ===\n\n")

  for (k in 1:length(health)) {
    subdat <- dat[dat$health == health[k], ]
    subdat$y <- subdat$teacher_report

    for (j in 1:length(all_model_sets)) {
      model_set <- all_model_sets[[j]]
      cat("Processing: health =", health[k], ", model set =", paste(model_set, collapse = ","), "\n")

      ps_specifications <- list(
        formula.list = full_ps_specifications$formula.list[model_set],
        h_x_names.list = full_ps_specifications$h_x_names.list[model_set],
        alpha_init.list = init.list[[k]][model_set],
        outcome = full_ps_specifications$outcome,
        inv_link = full_ps_specifications$inv_link
      )

      ebmr <- EBMRAlgorithmFast2$new("teacher_report", ps_specifications, subdat, W)

      # Define h_nu function for ensemble step
      h_nu <- function(data) cbind(f = data$father, p = data$parent_report, fp = data$fp)

      result <- ebmr$EBMR_IPW(h_nu = h_nu, true_ps = NULL)

      save_file <- paste0("MHD_results/EBMR_IPW_", subset_names[k], "_",
                          paste(model_set, collapse = ""), "_W_optimal_OR_test10.RDS")
      saveRDS(result, save_file)

      # Save full model results for summary
      if (j == 7) {  # Full model (1,2,3)
        if (k == 1) {
          ps_fit.list0 <- ebmr$ps_fit.list
          nu_hat0 <- result$nu.hat
          w_hat0 <- result$w.hat
        } else {
          ps_fit.list1 <- ebmr$ps_fit.list
          nu_hat1 <- result$nu.hat
          w_hat1 <- result$w.hat
        }
      }
    }
  }
}

print_header("PROPENSITY SCORE MODEL SUMMARIES")
print_ps_summary(ps_fit.list0, "Health = 0 (Unexposed)")
print_ps_summary(ps_fit.list1, "Health = 1 (Exposed)")

cat("\nNu estimates (full model 123):\n")
cat("Health = 0:", nu_hat0, "\n")
cat("Health = 1:", nu_hat1, "\n")

cat("\nW estimates (ensemble weights, full model 123):\n")
cat("Health = 0:", round(w_hat0, 4), "\n")
cat("Health = 1:", round(w_hat1, 4), "\n")

# Bootstrap functions (bootstrap, perturbation_bootstrap) are defined in MHD_functions.R

# #------------------------------------------------------------------------------#
# # Run Bootstrap (optional - set boot = TRUE to run)
# #------------------------------------------------------------------------------#
boot <- TRUE
perturb_boot <- TRUE
B <- 1000

# if (boot) {
#   boot_save_file <- "MHD_results/EBMR_IPW_bootstrap_W_optimal_OR_test10.RDS"

#   if (file.exists(boot_save_file)) {
#     cat("\n=== Standard Bootstrap results already exist, loading from file ===\n")
#     boot_result <- readRDS(boot_save_file)
#     cat("Loaded from:", boot_save_file, "\n")
#     cat("Number of bootstrap samples:", ncol(boot_result), "\n")
#     cat("Number of NA results:", sum(is.na(boot_result[1, ])), "\n")
#   } else {
#     cat("\n=== Running Standard Bootstrap (B =", B, ") ===\n")
#     dat$y <- dat$teacher_report

#     set.seed(20200525)
#     nu_init.list <- list(nu_hat0, nu_hat1)
#     boot_result <- bootstrap(full_ps_specifications, W, dat, B, init.list, nu_init.list)

#     saveRDS(boot_result, boot_save_file)

#     cat("\nStandard Bootstrap completed. Results saved to:", boot_save_file, "\n")
#     cat("Number of bootstrap samples:", ncol(boot_result), "\n")
#     cat("Number of NA results:", sum(is.na(boot_result[1, ])), "\n")
#   }
# }

# if (perturb_boot) {
#   perturb_save_file <- "MHD_results/EBMR_IPW_perturbation_bootstrap_W_optimal_OR_test10.RDS"

#   if (file.exists(perturb_save_file)) {
#     cat("\n=== Perturbation Bootstrap results already exist, loading from file ===\n")
#     perturb_result <- readRDS(perturb_save_file)
#     cat("Loaded from:", perturb_save_file, "\n")
#     cat("Number of perturbation samples:", ncol(perturb_result), "\n")
#     cat("Number of NA results:", sum(is.na(perturb_result[1, ])), "\n")
#   } else {
#     cat("\n=== Running Perturbation Bootstrap (B =", B, ") ===\n")
#     dat$y <- dat$teacher_report

#     set.seed(20200526)  # Different seed for perturbation bootstrap
#     nu_init.list <- list(nu_hat0, nu_hat1)
#     perturb_result <- perturbation_bootstrap(full_ps_specifications, W, dat, B, init.list, nu_init.list)

#     saveRDS(perturb_result, perturb_save_file)

#     cat("\nPerturbation Bootstrap completed. Results saved to:", perturb_save_file, "\n")
#     cat("Number of perturbation samples:", ncol(perturb_result), "\n")
#     cat("Number of NA results:", sum(is.na(perturb_result[1, ])), "\n")
#   }
# }

# # Compare standard and perturbation bootstrap results
# if (boot && perturb_boot) {
#   cat("\n=== Bootstrap Comparison ===\n")
#   cat("\nStandard Bootstrap:\n")
#   cat("  Valid samples:", sum(!is.na(boot_result[1, ])), "/", B, "\n")
#   cat("  Mean theta (model 123):", mean(boot_result[7, ], na.rm = TRUE), "\n")
#   cat("  SD theta (model 123):", sd(boot_result[7, ], na.rm = TRUE), "\n")

#   cat("\nPerturbation Bootstrap:\n")
#   cat("  Valid samples:", sum(!is.na(perturb_result[1, ])), "/", B, "\n")
#   cat("  Mean theta (model 123):", mean(perturb_result[7, ], na.rm = TRUE), "\n")
#   cat("  SD theta (model 123):", sd(perturb_result[7, ], na.rm = TRUE), "\n")
# }

# bootstrap_separate and perturbation_bootstrap_separate are defined in MHD_functions.R

# Run standard bootstrap for each health group and model set separately
if (boot) {
  cat("\n=== Running Standard Bootstrap for All Model Combinations ===\n")

  # Check if all files exist
  all_std_boot_files_exist <- TRUE
  for (k in 1:2) {
    for (j in 1:length(all_model_sets)) {
      save_file <- paste0("MHD_results/std_boot_", subset_names[k], "_",
                          model_set_labels[j], "_B", B, ".RDS")
      if (!file.exists(save_file)) {
        all_std_boot_files_exist <- FALSE
        break
      }
    }
    if (!all_std_boot_files_exist) break
  }

  if (all_std_boot_files_exist) {
    cat("All standard bootstrap files already exist, skipping computation.\n")
  } else {
    dat$y <- dat$teacher_report
    nu_init.list <- list(nu_hat0, nu_hat1)

    for (k in 1:2) {
      health_value <- health[k]
      cat("\nHealth =", health_value, ":\n")

      for (j in 1:length(all_model_sets)) {
        save_file <- paste0("MHD_results/std_boot_", subset_names[k], "_",
                            model_set_labels[j], "_B", B, ".RDS")

        if (file.exists(save_file)) {
          cat("  Model set", model_set_labels[j], "- already exists, skipping\n")
        } else {
          cat("  Model set", model_set_labels[j], "- running...\n")
          set.seed(20200528 + k * 10 + j)  # Different seed for each combination

          mu_replicates <- bootstrap_separate(
            full_ps_specifications, W, dat, B, init.list, nu_init.list,
            health_value, j, all_model_sets
          )

          saveRDS(mu_replicates, save_file)
          cat("    Saved:", save_file, "\n")
          cat("    Valid samples:", sum(!is.na(mu_replicates)), "/", B, "\n")
        }
      }
    }
  }
}

# Run perturbation bootstrap for each health group and model set separately
if (perturb_boot) {
  cat("\n=== Running Perturbation Bootstrap for All Model Combinations ===\n")

  # Check if all files exist
  all_perturb_files_exist <- TRUE
  for (k in 1:2) {
    for (j in 1:length(all_model_sets)) {
      save_file <- paste0("MHD_results/perturb_boot_", subset_names[k], "_",
                          model_set_labels[j], "_B", B, ".RDS")
      if (!file.exists(save_file)) {
        all_perturb_files_exist <- FALSE
        break
      }
    }
    if (!all_perturb_files_exist) break
  }

  if (all_perturb_files_exist) {
    cat("All perturbation bootstrap files already exist, skipping computation.\n")
  } else {
    dat$y <- dat$teacher_report
    nu_init.list <- list(nu_hat0, nu_hat1)

    for (k in 1:2) {
      health_value <- health[k]
      cat("\nHealth =", health_value, ":\n")

      for (j in 1:length(all_model_sets)) {
        save_file <- paste0("MHD_results/perturb_boot_", subset_names[k], "_",
                            model_set_labels[j], "_B", B, ".RDS")

        if (file.exists(save_file)) {
          cat("  Model set", model_set_labels[j], "- already exists, skipping\n")
        } else {
          cat("  Model set", model_set_labels[j], "- running...\n")
          set.seed(20200527 + k * 10 + j)  # Different seed for each combination

          mu_replicates <- perturbation_bootstrap_separate(
            full_ps_specifications, W, dat, B, init.list, nu_init.list,
            health_value, j, all_model_sets
          )

          saveRDS(mu_replicates, save_file)
          cat("    Saved:", save_file, "\n")
          cat("    Valid samples:", sum(!is.na(mu_replicates)), "/", B, "\n")
        }
      }
    }
  }
}

#------------------------------------------------------------------------------#
# Summarize results
#------------------------------------------------------------------------------#
# Pretty printing functions (print_header, print_section, format_ci)
# are defined in MHD_functions.R

print_header("EBMR ANALYSIS RESULTS")

# Load all individual health group results
# all_model_sets is already defined: 1, 2, 3, 12, 13, 23, 123
results_health0 <- list()
results_health1 <- list()
for (j in 1:length(all_model_sets)) {
  model_set <- all_model_sets[[j]]
  label <- model_set_labels[j]
  results_health0[[label]] <- readRDS(paste0("MHD_results/EBMR_IPW_health0_", label, "_W_optimal_OR_test10.RDS"))
  results_health1[[label]] <- readRDS(paste0("MHD_results/EBMR_IPW_health1_", label, "_W_optimal_OR_test10.RDS"))
}

# Compute all 49 combinations of theta (7 model sets for health=0 × 7 model sets for health=1)
# Plus complete case = 50 total
n_combinations <- length(all_model_sets)^2
all_results <- matrix(NA, n_combinations + 1, 2)
combination_labels <- character(n_combinations + 1)

idx <- 1
for (j0 in 1:length(all_model_sets)) {
  for (j1 in 1:length(all_model_sets)) {
    label0 <- model_set_labels[j0]
    label1 <- model_set_labels[j1]

    result0 <- results_health0[[label0]]
    result1 <- results_health1[[label1]]

    mu0 <- result0$mu_ipw
    se0 <- result0$se_ipw
    mu1 <- result1$mu_ipw
    se1 <- result1$se_ipw

    effect <- effect.f(mu0, mu1)
    Sigma <- matrix(c(se0^2, 0, 0, se1^2), 2, 2)
    se <- sqrt(delta_method(mu0, mu1) %*% Sigma %*% delta_method(mu0, mu1))

    all_results[idx, ] <- c(effect, se)
    # Create label like "theta_{101}^{110}" for model_set0 = {1,3} and model_set1 = {2,3}
    # Using binary representation: 1=model1, 2=model2, 3=model3
    binary0 <- paste0(as.integer(1 %in% all_model_sets[[j0]]),
                      as.integer(2 %in% all_model_sets[[j0]]),
                      as.integer(3 %in% all_model_sets[[j0]]))
    binary1 <- paste0(as.integer(1 %in% all_model_sets[[j1]]),
                      as.integer(2 %in% all_model_sets[[j1]]),
                      as.integer(3 %in% all_model_sets[[j1]]))
    combination_labels[idx] <- paste0("theta_{", binary0, "}^{", binary1, "}")
    idx <- idx + 1
  }
}

# Complete case results
cc1 <- dat$health == 1 & dat$r == 1
cc0 <- dat$health == 0 & dat$r == 1
mu_cc1 <- mean(dat[cc1, ]$teacher_report)
se_cc1 <- sd(dat[cc1, ]$teacher_report) / sqrt(sum(cc1))
mu_cc0 <- mean(dat[cc0, ]$teacher_report)
se_cc0 <- sd(dat[cc0, ]$teacher_report) / sqrt(sum(cc0))
Sigma <- matrix(c(se_cc0^2, 0, 0, se_cc1^2), 2, 2)
all_results[n_combinations + 1, ] <- c(effect.f(mu_cc0, mu_cc1),
                      sqrt(delta_method(mu_cc0, mu_cc1) %*% Sigma %*% delta_method(mu_cc0, mu_cc1)))
combination_labels[n_combinations + 1] <- "Complete Case"

# Create results data frame
full_results <- as.data.frame(all_results)
colnames(full_results) <- c("PE", "SE")
rownames(full_results) <- combination_labels

# Add 95% CI
full_results$CI_lower <- full_results$PE - qnorm(0.975) * full_results$SE
full_results$CI_upper <- full_results$PE + qnorm(0.975) * full_results$SE

# Identify diagonal combinations (same model set for both health groups)
# These are at indices: 1, 9, 17, 25, 33, 41, 49 (i.e., j0 == j1)
diagonal_indices <- seq(1, n_combinations, by = length(all_model_sets) + 1)
# Actually: 1, 1+8=9 is wrong. Let me recalculate.
# For j0=1,j1=1: idx=1; j0=1,j1=2: idx=2; ... j0=1,j1=7: idx=7
# For j0=2,j1=1: idx=8; j0=2,j1=2: idx=9; ...
# So diagonal is at: 1, 9, 17, 25, 33, 41, 49
diagonal_indices <- c(1, 9, 17, 25, 33, 41, 49)

# Load bootstrap SEs for diagonal combinations (same model set for both health groups)
# Bootstrap only computes diagonal cases (7 model combinations + CC)
boot_effect_se <- rep(NA, 8)
boot_valid_n <- NA
if (file.exists("MHD_results/EBMR_IPW_bootstrap_W_optimal_OR_test10.RDS")) {
  boot_result <- readRDS("MHD_results/EBMR_IPW_bootstrap_W_optimal_OR_test10.RDS")
  boot_valid_n <- sum(!is.na(boot_result[1, ]))
  for (k in 1:7) {
    theta.hat <- paste0("theta.hat", k)
    boot_effect_se[k] <- sd(boot_result[theta.hat, ], na.rm = TRUE)
  }
  boot_effect_se[8] <- sd(rm.extreme(boot_result["theta_cc", !is.na(boot_result["theta_cc", ])]))
}

# Load standard bootstrap SEs for all 49 combinations (if separate files exist)
boot_all_se <- rep(NA, n_combinations + 1)  # 49 combinations + complete case
boot_all_valid_n <- NA

# Check for separate standard bootstrap files (14 files: 7 for health=0, 7 for health=1)
all_separate_boot_exist <- TRUE
for (k in 1:2) {
  for (j in 1:length(all_model_sets)) {
    check_file <- paste0("MHD_results/std_boot_", subset_names[k], "_",
                         model_set_labels[j], "_B", B, ".RDS")
    if (!file.exists(check_file)) {
      all_separate_boot_exist <- FALSE
      break
    }
  }
  if (!all_separate_boot_exist) break
}

if (all_separate_boot_exist) {
  cat("Loading separate standard bootstrap files for all 49 combinations...\n")

  # Load all mu replicates for each health group and model set
  boot_mu_replicates_health0 <- list()
  boot_mu_replicates_health1 <- list()

  for (j in 1:length(all_model_sets)) {
    file0 <- paste0("MHD_results/std_boot_health0_", model_set_labels[j], "_B", B, ".RDS")
    file1 <- paste0("MHD_results/std_boot_health1_", model_set_labels[j], "_B", B, ".RDS")
    boot_mu_replicates_health0[[j]] <- readRDS(file0)
    boot_mu_replicates_health1[[j]] <- readRDS(file1)
  }

  # Compute theta SE for all 49 combinations
  idx <- 1
  for (j0 in 1:length(all_model_sets)) {
    for (j1 in 1:length(all_model_sets)) {
      mu0_reps <- boot_mu_replicates_health0[[j0]]
      mu1_reps <- boot_mu_replicates_health1[[j1]]

      # Compute theta replicate-by-replicate
      # Handle NAs: only use replicates where both are valid
      valid_idx <- !is.na(mu0_reps) & !is.na(mu1_reps)
      if (sum(valid_idx) > 10) {
        theta_reps <- effect.f(mu0_reps[valid_idx], mu1_reps[valid_idx])
        boot_all_se[idx] <- sd(theta_reps, na.rm = TRUE)
      }
      idx <- idx + 1
    }
  }

  # Complete case SE (use diagonal bootstrap if available)
  if (file.exists("MHD_results/EBMR_IPW_bootstrap_W_optimal_OR_test10.RDS")) {
    boot_result <- readRDS("MHD_results/EBMR_IPW_bootstrap_W_optimal_OR_test10.RDS")
    boot_all_se[n_combinations + 1] <- sd(rm.extreme(boot_result["theta_cc", !is.na(boot_result["theta_cc", ])]))
  }

  boot_all_valid_n <- sum(!is.na(boot_mu_replicates_health0[[1]]))
  cat("  Standard bootstrap SE computed for all 49 combinations\n")
  cat("  Valid samples:", boot_all_valid_n, "/", B, "\n")

  # Update boot_valid_n if not already set
  if (is.na(boot_valid_n)) {
    boot_valid_n <- boot_all_valid_n
  }
}

# Load perturbation bootstrap SEs for all 49 combinations
# First check if separate bootstrap files exist for all model combinations
perturb_all_se <- rep(NA, n_combinations + 1)  # 49 combinations + complete case
perturb_valid_n <- NA

# Check for separate perturbation bootstrap files (14 files: 7 for health=0, 7 for health=1)
all_separate_perturb_exist <- TRUE
for (k in 1:2) {
  for (j in 1:length(all_model_sets)) {
    check_file <- paste0("MHD_results/perturb_boot_", subset_names[k], "_",
                         model_set_labels[j], "_B", B, ".RDS")
    if (!file.exists(check_file)) {
      all_separate_perturb_exist <- FALSE
      break
    }
  }
  if (!all_separate_perturb_exist) break
}

if (all_separate_perturb_exist) {
  cat("Loading separate perturbation bootstrap files for all 49 combinations...\n")

  # Load all mu replicates for each health group and model set
  mu_replicates_health0 <- list()
  mu_replicates_health1 <- list()

  for (j in 1:length(all_model_sets)) {
    file0 <- paste0("MHD_results/perturb_boot_health0_", model_set_labels[j], "_B", B, ".RDS")
    file1 <- paste0("MHD_results/perturb_boot_health1_", model_set_labels[j], "_B", B, ".RDS")
    mu_replicates_health0[[j]] <- readRDS(file0)
    mu_replicates_health1[[j]] <- readRDS(file1)
  }

  # Compute theta SE for all 49 combinations
  idx <- 1
  for (j0 in 1:length(all_model_sets)) {
    for (j1 in 1:length(all_model_sets)) {
      mu0_reps <- mu_replicates_health0[[j0]]
      mu1_reps <- mu_replicates_health1[[j1]]

      # Compute theta replicate-by-replicate
      # Handle NAs: only use replicates where both are valid
      valid_idx <- !is.na(mu0_reps) & !is.na(mu1_reps)
      if (sum(valid_idx) > 10) {
        theta_reps <- effect.f(mu0_reps[valid_idx], mu1_reps[valid_idx])
        perturb_all_se[idx] <- sd(theta_reps, na.rm = TRUE)
      }
      idx <- idx + 1
    }
  }

  # Complete case SE (use diagonal bootstrap if available)
  if (file.exists("MHD_results/EBMR_IPW_perturbation_bootstrap_W_optimal_OR_test10.RDS")) {
    perturb_result <- readRDS("MHD_results/EBMR_IPW_perturbation_bootstrap_W_optimal_OR_test10.RDS")
    perturb_all_se[n_combinations + 1] <- sd(rm.extreme(perturb_result["theta_cc", !is.na(perturb_result["theta_cc", ])]))
  }

  perturb_valid_n <- sum(!is.na(mu_replicates_health0[[1]]))
  cat("  Perturbation bootstrap SE computed for all 49 combinations\n")
  cat("  Valid samples:", perturb_valid_n, "/", B, "\n")

} else if (file.exists("MHD_results/EBMR_IPW_perturbation_bootstrap_W_optimal_OR_test10.RDS")) {
  # Fall back to diagonal-only bootstrap
  cat("Loading diagonal-only perturbation bootstrap...\n")
  perturb_result <- readRDS("MHD_results/EBMR_IPW_perturbation_bootstrap_W_optimal_OR_test10.RDS")
  perturb_valid_n <- sum(!is.na(perturb_result[1, ]))

  # Only diagonal indices have SE
  for (k in 1:7) {
    theta.hat <- paste0("theta.hat", k)
    perturb_all_se[diagonal_indices[k]] <- sd(perturb_result[theta.hat, ], na.rm = TRUE)
  }
  perturb_all_se[n_combinations + 1] <- sd(rm.extreme(perturb_result["theta_cc", !is.na(perturb_result["theta_cc", ])]))
}

# For backward compatibility, also create perturb_effect_se for diagonal cases
perturb_effect_se <- rep(NA, 8)
for (k in 1:7) {
  perturb_effect_se[k] <- perturb_all_se[diagonal_indices[k]]
}
perturb_effect_se[8] <- perturb_all_se[n_combinations + 1]

print_section("Log Odds Ratio Estimates (theta = log OR)")

cat("\n  Notation: theta_{abc}^{xyz} where abc = models for health=0, xyz = models for health=1\n")
cat("  Binary encoding: position 1=Model1(MNAR:father), 2=Model2(MNAR:parent), 3=Model3(MAR)\n\n")

# Show bootstrap sample sizes if available
if (!is.na(boot_valid_n) || !is.na(perturb_valid_n)) {
  cat("  Bootstrap sample sizes:\n")
  if (!is.na(boot_valid_n)) {
    cat(sprintf("    Standard Bootstrap: %d / %d valid samples\n", boot_valid_n, B))
  }
  if (!is.na(perturb_valid_n)) {
    cat(sprintf("    Perturbation Bootstrap: %d / %d valid samples\n", perturb_valid_n, B))
  }
  cat("\n")
}

# Show all 49 model combinations
cat("  === All Model Combinations (49 total) ===\n\n")

# Determine header based on available bootstrap data
has_boot <- !is.na(boot_valid_n)
has_perturb <- !is.na(perturb_valid_n)
has_boot_all <- all_separate_boot_exist  # Check if we have SE for all 49

# Model set descriptions for annotations
model_set_desc <- c(
  "{1}", "{2}", "{3}", "{1,2}", "{1,3}", "{2,3}", "{1,2,3}"
)

# Header with Boot SE and/or Perturb SE columns if available
if (has_boot_all && has_perturb) {
  cat(sprintf("  %-25s %10s %10s %12s %12s   %-22s\n", "Estimator", "Estimate", "SE", "Boot SE", "Perturb SE", "95% CI"))
  cat("  ", paste(rep("-", 105), collapse = ""), "\n")
} else if (has_boot_all) {
  cat(sprintf("  %-25s %10s %10s %12s   %-22s\n", "Estimator", "Estimate", "SE", "Boot SE", "95% CI"))
  cat("  ", paste(rep("-", 90), collapse = ""), "\n")
} else if (has_perturb) {
  cat(sprintf("  %-25s %10s %10s %12s   %-22s\n", "Estimator", "Estimate", "SE", "Perturb SE", "95% CI"))
  cat("  ", paste(rep("-", 90), collapse = ""), "\n")
} else {
  cat(sprintf("  %-25s %10s %10s   %-22s\n", "Estimator", "Estimate", "SE", "95% CI"))
  cat("  ", paste(rep("-", 75), collapse = ""), "\n")
}

# Loop through all 49 combinations
for (idx in 1:n_combinations) {
  ci <- format_ci(full_results$CI_lower[idx], full_results$CI_upper[idx])

  # Check if this is a diagonal case (same model set for both health groups)
  is_diagonal <- idx %in% diagonal_indices
  diagonal_marker <- ifelse(is_diagonal, "*", " ")

  # Get bootstrap SEs (now available for all 49 combinations if separate bootstrap was run)
  if (has_boot_all && has_perturb) {
    boot_se_str <- ifelse(!is.na(boot_all_se[idx]), sprintf("%12.4f", boot_all_se[idx]), "          --")
    perturb_se_str <- ifelse(!is.na(perturb_all_se[idx]), sprintf("%12.4f", perturb_all_se[idx]), "          --")
    cat(sprintf(" %s%-24s %10.4f %10.4f %s %s   %-22s\n",
                diagonal_marker,
                combination_labels[idx],
                full_results$PE[idx],
                full_results$SE[idx],
                boot_se_str,
                perturb_se_str,
                ci))
  } else if (has_boot_all) {
    boot_se_str <- ifelse(!is.na(boot_all_se[idx]), sprintf("%12.4f", boot_all_se[idx]), "          --")
    cat(sprintf(" %s%-24s %10.4f %10.4f %s   %-22s\n",
                diagonal_marker,
                combination_labels[idx],
                full_results$PE[idx],
                full_results$SE[idx],
                boot_se_str,
                ci))
  } else if (has_perturb) {
    perturb_se_str <- ifelse(!is.na(perturb_all_se[idx]), sprintf("%12.4f", perturb_all_se[idx]), "          --")
    cat(sprintf(" %s%-24s %10.4f %10.4f %s   %-22s\n",
                diagonal_marker,
                combination_labels[idx],
                full_results$PE[idx],
                full_results$SE[idx],
                perturb_se_str,
                ci))
  } else {
    cat(sprintf(" %s%-24s %10.4f %10.4f   %-22s\n",
                diagonal_marker,
                combination_labels[idx],
                full_results$PE[idx],
                full_results$SE[idx],
                ci))
  }
}

if (has_boot_all && has_perturb) {
  cat("  ", paste(rep("-", 105), collapse = ""), "\n")
} else if (has_boot_all || has_perturb) {
  cat("  ", paste(rep("-", 90), collapse = ""), "\n")
} else {
  cat("  ", paste(rep("-", 75), collapse = ""), "\n")
}
cat("  * = Diagonal case (same model set for both health groups)\n")
if (all_separate_boot_exist) {
  cat("  Boot SE computed from separate standard bootstrap for each health group and model set.\n")
}
if (all_separate_perturb_exist) {
  cat("  Perturb SE computed from separate perturbation bootstrap for each health group and model set.\n\n")
} else if (!all_separate_boot_exist && !has_perturb) {
  cat("  Note: Bootstrap SE only available for diagonal cases (run with boot=TRUE and std_boot=TRUE to get all 49).\n\n")
} else {
  cat("\n")
}

# Show diagonal cases with bootstrap SE if available
if (has_boot || has_perturb) {
  cat("  === Diagonal Cases with Bootstrap SE ===\n\n")

  if (has_boot && has_perturb) {
    cat(sprintf("  %-25s %10s %10s %10s %10s   %-22s\n", "Estimator", "Estimate", "SE", "Boot SE", "Perturb SE", "95% CI"))
    cat("  ", paste(rep("-", 105), collapse = ""), "\n")
  } else if (has_perturb) {
    cat(sprintf("  %-25s %10s %10s %10s   %-22s\n", "Estimator", "Estimate", "SE", "Perturb SE", "95% CI"))
    cat("  ", paste(rep("-", 90), collapse = ""), "\n")
  } else {
    cat(sprintf("  %-25s %10s %10s %10s   %-22s\n", "Estimator", "Estimate", "SE", "Boot SE", "95% CI"))
    cat("  ", paste(rep("-", 90), collapse = ""), "\n")
  }

  for (i in 1:7) {
    row_idx <- diagonal_indices[i]
    ci <- format_ci(full_results$CI_lower[row_idx], full_results$CI_upper[row_idx])

    if (has_boot && has_perturb) {
      boot_se_str <- ifelse(is.na(boot_effect_se[i]), "        --", sprintf("%10.4f", boot_effect_se[i]))
      perturb_se_str <- ifelse(is.na(perturb_effect_se[i]), "        --", sprintf("%10.4f", perturb_effect_se[i]))
      cat(sprintf("  %-25s %10.4f %10.4f %s %s   %-22s\n",
                  combination_labels[row_idx], full_results$PE[row_idx], full_results$SE[row_idx],
                  boot_se_str, perturb_se_str, ci))
    } else if (has_perturb) {
      perturb_se_str <- ifelse(is.na(perturb_effect_se[i]), "        --", sprintf("%10.4f", perturb_effect_se[i]))
      cat(sprintf("  %-25s %10.4f %10.4f %s   %-22s\n",
                  combination_labels[row_idx], full_results$PE[row_idx], full_results$SE[row_idx],
                  perturb_se_str, ci))
    } else {
      boot_se_str <- ifelse(is.na(boot_effect_se[i]), "        --", sprintf("%10.4f", boot_effect_se[i]))
      cat(sprintf("  %-25s %10.4f %10.4f %s   %-22s\n",
                  combination_labels[row_idx], full_results$PE[row_idx], full_results$SE[row_idx],
                  boot_se_str, ci))
    }
  }
}

# Show highlighted cases: Model 3 for health=0 AND Model 1 for health=1
# Model 3 included (health=0): rows 3, 5, 6, 7 (sets: "3", "13", "23", "123")
# Model 1 included (health=1): cols 1, 4, 5, 7 (sets: "1", "12", "13", "123")
model3_rows <- c(3, 5, 6, 7)
model1_cols <- c(1, 4, 5, 7)

cat("\n  === Highlighted Cases: Model 3 (health=0) + Model 1 (health=1) ===\n")
cat("  (MAR model for unexposed group, MNAR:father for exposed group)\n\n")

if (has_boot_all && has_perturb) {
  cat(sprintf("  %-25s %10s %10s %10s %10s   %-22s\n", "Estimator", "Estimate", "SE", "Boot SE", "Perturb SE", "95% CI"))
  cat("  ", paste(rep("-", 105), collapse = ""), "\n")
} else if (has_boot_all) {
  cat(sprintf("  %-25s %10s %10s %10s   %-22s\n", "Estimator", "Estimate", "SE", "Boot SE", "95% CI"))
  cat("  ", paste(rep("-", 90), collapse = ""), "\n")
} else if (has_perturb) {
  cat(sprintf("  %-25s %10s %10s %10s   %-22s\n", "Estimator", "Estimate", "SE", "Perturb SE", "95% CI"))
  cat("  ", paste(rep("-", 90), collapse = ""), "\n")
} else if (has_boot) {
  cat(sprintf("  %-25s %10s %10s %10s   %-22s\n", "Estimator", "Estimate", "SE", "Boot SE", "95% CI"))
  cat("  ", paste(rep("-", 90), collapse = ""), "\n")
} else {
  cat(sprintf("  %-25s %10s %10s   %-22s\n", "Estimator", "Estimate", "SE", "95% CI"))
  cat("  ", paste(rep("-", 75), collapse = ""), "\n")
}

# Loop through highlighted combinations (16 total: 4 rows x 4 cols)
for (row in model3_rows) {
  for (col in model1_cols) {
    # Calculate linear index: (row-1)*7 + col
    idx <- (row - 1) * 7 + col
    ci <- format_ci(full_results$CI_lower[idx], full_results$CI_upper[idx])

    # Mark diagonal cases
    is_diagonal <- row == col
    diagonal_marker <- ifelse(is_diagonal, "*", " ")

    if (has_boot_all && has_perturb) {
      # Get boot SE from boot_all_se (available for all 49 combinations)
      boot_se_str <- ifelse(!is.na(boot_all_se[idx]), sprintf("%10.4f", boot_all_se[idx]), "        --")
      perturb_se_str <- ifelse(!is.na(perturb_all_se[idx]), sprintf("%10.4f", perturb_all_se[idx]), "        --")
      cat(sprintf(" %s%-24s %10.4f %10.4f %s %s   %-22s\n",
                  diagonal_marker,
                  combination_labels[idx], full_results$PE[idx], full_results$SE[idx],
                  boot_se_str, perturb_se_str, ci))
    } else if (has_boot_all) {
      boot_se_str <- ifelse(!is.na(boot_all_se[idx]), sprintf("%10.4f", boot_all_se[idx]), "        --")
      cat(sprintf(" %s%-24s %10.4f %10.4f %s   %-22s\n",
                  diagonal_marker,
                  combination_labels[idx], full_results$PE[idx], full_results$SE[idx],
                  boot_se_str, ci))
    } else if (has_perturb) {
      perturb_se_str <- ifelse(!is.na(perturb_all_se[idx]), sprintf("%10.4f", perturb_all_se[idx]), "        --")
      cat(sprintf(" %s%-24s %10.4f %10.4f %s   %-22s\n",
                  diagonal_marker,
                  combination_labels[idx], full_results$PE[idx], full_results$SE[idx],
                  perturb_se_str, ci))
    } else if (has_boot) {
      # Fall back to diagonal-only boot SE
      boot_se_str <- "        --"
      if (is_diagonal) {
        diag_k <- which(diagonal_indices == idx)
        if (length(diag_k) > 0 && !is.na(boot_effect_se[diag_k])) {
          boot_se_str <- sprintf("%10.4f", boot_effect_se[diag_k])
        }
      }
      cat(sprintf(" %s%-24s %10.4f %10.4f %s   %-22s\n",
                  diagonal_marker,
                  combination_labels[idx], full_results$PE[idx], full_results$SE[idx],
                  boot_se_str, ci))
    } else {
      cat(sprintf(" %s%-24s %10.4f %10.4f   %-22s\n",
                  diagonal_marker,
                  combination_labels[idx], full_results$PE[idx], full_results$SE[idx], ci))
    }
  }
}

if (has_boot_all || has_boot || has_perturb) {
  if ((has_boot_all && has_perturb) || (has_boot && has_perturb)) {
    cat("  ", paste(rep("-", 105), collapse = ""), "\n")
  } else {
    cat("  ", paste(rep("-", 90), collapse = ""), "\n")
  }
  cat("  * = Diagonal case (same model set for both health groups)\n")
}

# Complete case - show in main table
cat("\n  === Complete Case Analysis ===\n\n")
cc_ci <- format_ci(full_results$CI_lower[n_combinations + 1], full_results$CI_upper[n_combinations + 1])

# Use boot_all_se if available, otherwise fall back to boot_effect_se
cc_boot_se <- ifelse(has_boot_all && !is.na(boot_all_se[n_combinations + 1]),
                     boot_all_se[n_combinations + 1],
                     boot_effect_se[8])
cc_perturb_se <- perturb_all_se[n_combinations + 1]

if ((has_boot_all || has_boot) && has_perturb) {
  boot_se_str <- ifelse(is.na(cc_boot_se), "        --", sprintf("%10.4f", cc_boot_se))
  perturb_se_str <- ifelse(is.na(cc_perturb_se), "        --", sprintf("%10.4f", cc_perturb_se))
  cat(sprintf("  %-25s %10s %10s %10s %10s   %-22s\n", "Estimator", "Estimate", "SE", "Boot SE", "Perturb SE", "95% CI"))
  cat("  ", paste(rep("-", 105), collapse = ""), "\n")
  cat(sprintf("  %-25s %10.4f %10.4f %s %s   %-22s\n",
              "Complete Case",
              full_results$PE[n_combinations + 1],
              full_results$SE[n_combinations + 1],
              boot_se_str, perturb_se_str, cc_ci))
} else if (has_perturb) {
  perturb_se_str <- ifelse(is.na(cc_perturb_se), "        --", sprintf("%10.4f", cc_perturb_se))
  cat(sprintf("  %-25s %10s %10s %10s   %-22s\n", "Estimator", "Estimate", "SE", "Perturb SE", "95% CI"))
  cat("  ", paste(rep("-", 90), collapse = ""), "\n")
  cat(sprintf("  %-25s %10.4f %10.4f %s   %-22s\n",
              "Complete Case",
              full_results$PE[n_combinations + 1],
              full_results$SE[n_combinations + 1],
              perturb_se_str, cc_ci))
} else if (has_boot_all || has_boot) {
  boot_se_str <- ifelse(is.na(cc_boot_se), "        --", sprintf("%10.4f", cc_boot_se))
  cat(sprintf("  %-25s %10s %10s %10s   %-22s\n", "Estimator", "Estimate", "SE", "Boot SE", "95% CI"))
  cat("  ", paste(rep("-", 90), collapse = ""), "\n")
  cat(sprintf("  %-25s %10.4f %10.4f %s   %-22s\n",
              "Complete Case",
              full_results$PE[n_combinations + 1],
              full_results$SE[n_combinations + 1],
              boot_se_str, cc_ci))
} else {
  cat(sprintf("  %-25s %10s %10s   %-22s\n", "Estimator", "Estimate", "SE", "95% CI"))
  cat("  ", paste(rep("-", 75), collapse = ""), "\n")
  cat(sprintf("  %-25s %10.4f %10.4f   %-22s\n",
              "Complete Case",
              full_results$PE[n_combinations + 1],
              full_results$SE[n_combinations + 1], cc_ci))
}

# Show full matrix of all 49 combinations
cat("\n  === Full Results Matrix (all 49 combinations) ===\n")
cat("  Rows: model set for health=0, Columns: model set for health=1\n\n")

# Create matrix view
result_matrix <- matrix(full_results$PE[1:n_combinations], nrow = 7, ncol = 7, byrow = TRUE)
rownames(result_matrix) <- model_set_labels
colnames(result_matrix) <- model_set_labels

cat("  Point Estimates:\n")
cat("  ", sprintf("%8s", ""), " ")
for (j in 1:7) cat(sprintf("%8s", model_set_labels[j]))
cat("\n")
for (i in 1:7) {
  cat("  ", sprintf("%8s", model_set_labels[i]), " ")
  for (j in 1:7) {
    cat(sprintf("%8.4f", result_matrix[i, j]))
  }
  cat("\n")
}

#------------------------------------------------------------------------------#
# Generate heatmap visualization
#------------------------------------------------------------------------------#
cat("\n  Generating heatmap visualization...\n")

# Identify which model sets include Model 3 (for health=0) and Model 1 (for health=1)
# Model sets: 1="1", 2="2", 3="3", 4="12", 5="13", 6="23", 7="123"
# Model 3 included: indices 3, 5, 6, 7 (sets: "3", "13", "23", "123")
# Model 1 included: indices 1, 4, 5, 7 (sets: "1", "12", "13", "123")
model3_included <- c(3, 5, 6, 7)  # Row indices (health=0)
model1_included <- c(1, 4, 5, 7)  # Column indices (health=1)

# Create highlight matrix: 1 = both conditions met, 0.5 = one condition, 0 = neither
highlight_matrix <- matrix(0, 7, 7)
for (i in 1:7) {
  for (j in 1:7) {
    has_model3 <- i %in% model3_included
    has_model1 <- j %in% model1_included
    if (has_model3 && has_model1) {
      highlight_matrix[i, j] <- 2  # Both conditions
    } else if (has_model3 || has_model1) {
      highlight_matrix[i, j] <- 1  # One condition
    }
  }
}

# Save heatmap as PNG
png("MHD_results/theta_heatmap.png", width = 10, height = 8, units = "in", res = 300)

# Set up layout for heatmap with legend
layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1))

# Color palette for non-highlighted cells: muted gray tones
# Highlighted cells (both conditions): vibrant green-to-purple scale
n_colors <- 100
pe_range <- range(result_matrix)
max_abs <- max(abs(pe_range))
breaks <- seq(-max_abs, max_abs, length.out = n_colors + 1)

# Muted color palette for non-highlighted cells (desaturated)
gray_blue <- colorRampPalette(c("#8C96A0", "#A8B0B8", "#C4CCD4", "#E0E4E8", "#F0F0F0"))
gray_red <- colorRampPalette(c("#F0F0F0", "#E8E0DC", "#D4C4BC", "#C0A8A0", "#A08C88"))
muted_colors <- c(gray_blue(n_colors/2), gray_red(n_colors/2))

# Vibrant color palette for highlighted cells (Model3 in health=0 AND Model1 in health=1)
# Use a striking gold/orange color scheme
highlight_blue <- colorRampPalette(c("#1A5276", "#2874A6", "#3498DB", "#85C1E9", "#FFFFFF"))
highlight_orange <- colorRampPalette(c("#FFFFFF", "#FAD7A0", "#F5B041", "#E67E22", "#D35400"))
highlight_colors <- c(highlight_blue(n_colors/2), highlight_orange(n_colors/2))

# Map values to colors
get_color <- function(val, highlighted) {
  idx <- findInterval(val, breaks, all.inside = TRUE)
  if (highlighted) {
    return(highlight_colors[idx])
  } else {
    return(muted_colors[idx])
  }
}

# Labels for axes
row_labels <- c("{1}", "{2}", "{3}", "{1,2}", "{1,3}", "{2,3}", "{1,2,3}")
col_labels <- row_labels

# Main heatmap
par(mar = c(6, 7, 5, 1))
plot(NULL, xlim = c(0.5, 7.5), ylim = c(0.5, 7.5), xlab = "", ylab = "",
     xaxt = "n", yaxt = "n", bty = "n", asp = 1)

# First pass: Draw non-highlighted cells (faded background)
for (i in 1:7) {
  for (j in 1:7) {
    y_pos <- 8 - i  # Flip y-axis
    is_highlighted <- highlight_matrix[i, j] == 2

    if (!is_highlighted) {
      cell_col <- get_color(result_matrix[i, j], FALSE)
      rect(j - 0.5, y_pos - 0.5, j + 0.5, y_pos + 0.5,
           col = cell_col, border = "#CCCCCC", lwd = 0.5)
      # Muted text for non-highlighted
      text(j, y_pos, sprintf("%.3f", result_matrix[i, j]), cex = 0.75, col = "#666666", font = 1)
    }
  }
}

# Second pass: Draw highlighted cells on top with emphasis
for (i in 1:7) {
  for (j in 1:7) {
    y_pos <- 8 - i  # Flip y-axis
    is_highlighted <- highlight_matrix[i, j] == 2

    if (is_highlighted) {
      cell_col <- get_color(result_matrix[i, j], TRUE)

      # Draw outer glow effect (golden border)
      rect(j - 0.52, y_pos - 0.52, j + 0.52, y_pos + 0.52,
           col = NA, border = "#FFD700", lwd = 6)

      # Draw cell with thick black border
      rect(j - 0.5, y_pos - 0.5, j + 0.5, y_pos + 0.5,
           col = cell_col, border = "#000000", lwd = 3)

      # Bold, larger text for highlighted cells
      text_col <- ifelse(abs(result_matrix[i, j]) > max_abs * 0.5, "white", "black")
      text(j, y_pos, sprintf("%.3f", result_matrix[i, j]), cex = 1.0, col = text_col, font = 2)
    }
  }
}

# Add axes labels - highlight rows with Model 3 and columns with Model 1
axis_col_default <- "gray40"
axis_col_highlight <- "#D35400"  # Orange for highlighted

# X-axis (columns - Model 1 highlighted)
for (j in 1:7) {
  col <- ifelse(j %in% model1_included, axis_col_highlight, axis_col_default)
  font_style <- ifelse(j %in% model1_included, 2, 1)
  mtext(col_labels[j], side = 1, at = j, line = 0.5, cex = 1.0, col = col, font = font_style)
}

# Y-axis (rows - Model 3 highlighted)
for (i in 1:7) {
  y_pos <- 8 - i
  col <- ifelse(i %in% model3_included, axis_col_highlight, axis_col_default)
  font_style <- ifelse(i %in% model3_included, 2, 1)
  mtext(row_labels[i], side = 2, at = y_pos, line = 0.5, cex = 1.0, col = col, font = font_style, las = 1)
}

# Add axis titles
mtext("Candidate set for Health = 1", side = 1, line = 4.5, cex = 1.2, font = 2)
mtext("Candidate set for Health = 0", side = 2, line = 5, cex = 1.2, font = 2)

# Title
mtext(expression(paste("Log Odds Ratio Estimates (", hat(theta), ") by Model Combination")),
      side = 3, line = 2.5, cex = 1.4, font = 2)
# mtext("Highlighted: Model 3 (MAR) for Health=0, Model 1 (MNAR:father) for Health=1",
#       side = 3, line = 0.8, cex = 0.95, font = 3)

# Color legend - show highlighted color scale
par(mar = c(6, 1, 5, 3))
legend_y <- seq(0, 1, length.out = n_colors)
image(1, legend_y, matrix(legend_y, nrow = 1), col = highlight_colors,
      xaxt = "n", yaxt = "n", xlab = "", ylab = "")

# Add border to legend
box(lwd = 2)

# Legend axis
legend_ticks <- pretty(c(-max_abs, max_abs), n = 5)
legend_ticks <- legend_ticks[legend_ticks >= -max_abs & legend_ticks <= max_abs]
legend_pos <- (legend_ticks + max_abs) / (2 * max_abs)
axis(4, at = legend_pos, labels = sprintf("%.2f", legend_ticks), las = 1, cex.axis = 0.9)
mtext(expression(hat(theta)), side = 4, line = 2, cex = 1.1)
# mtext("(Highlighted cells)", side = 1, line = 1, cex = 0.8, font = 3)

dev.off()

cat("  Heatmap saved to: MHD_results/theta_heatmap.png\n")

# Summary statistics
print_section("Summary Statistics")
cat(sprintf("\n  %-40s %10.4f\n", "Mean PE (diagonal, same models):", mean(full_results$PE[diagonal_indices])))
cat(sprintf("  %-40s %10.4f\n", "SD of PE (diagonal, same models):", sd(full_results$PE[diagonal_indices])))
cat(sprintf("  %-40s %10.4f\n", "Mean PE (all 49 combinations):", mean(full_results$PE[1:n_combinations])))
cat(sprintf("  %-40s %10.4f\n", "SD of PE (all 49 combinations):", sd(full_results$PE[1:n_combinations])))
cat(sprintf("  %-40s %10.4f\n", "Complete Case PE:", full_results$PE[n_combinations + 1]))
cat(sprintf("  %-40s %10.4f\n", "Range of PE (all combinations):",
            max(full_results$PE[1:n_combinations]) - min(full_results$PE[1:n_combinations])))

#------------------------------------------------------------------------------#
# Final comparison table
#------------------------------------------------------------------------------#
print_section("Comparison: MNAR (EBMR) vs Complete Case")

# Index for full model (1,2,3) for both health groups = theta_{111}^{111}
full_model_idx <- 49  # Last diagonal element

cat(sprintf("\n  %-30s %12s %12s %24s\n", "Method", "Estimate", "SE", "95% CI"))
cat("  ", paste(rep("-", 82), collapse = ""), "\n")
cat(sprintf("  %-30s %12.4f %12.4f %24s\n",
            "EBMR (All models, both)",
            full_results$PE[full_model_idx],
            full_results$SE[full_model_idx],
            format_ci(full_results$CI_lower[full_model_idx], full_results$CI_upper[full_model_idx])))
cat(sprintf("  %-30s %12.4f %12.4f %24s\n",
            "Complete Case",
            full_results$PE[n_combinations + 1],
            full_results$SE[n_combinations + 1],
            format_ci(full_results$CI_lower[n_combinations + 1], full_results$CI_upper[n_combinations + 1])))
cat("  ", paste(rep("-", 82), collapse = ""), "\n")

# Interpretation
print_section("Interpretation")
cat("\n  Note: theta represents the log odds ratio comparing health=1 vs health=0\n")
cat("        for the teacher report outcome.\n\n")
cat("  - EBMR accounts for potential MNAR mechanism\n")
cat("  - Complete Case ignores missing data mechanism\n")

#------------------------------------------------------------------------------#
# Sensitivity Analysis
#------------------------------------------------------------------------------#
print_header("SENSITIVITY ANALYSIS")

cat("Performing sensitivity analysis with locally misspecified models...\n\n")

# Configuration for sensitivity analysis
possibly_true_ps.vector <- c(3, 1)  # Which PS model to perturb for each health group
exp_tilt_x_names <- c("father", "parent_report")

# Define xi range for each health group (based on n^{-1/2})
xi.vector.list <- list(
  xi.vector0 = seq(0, nrow(dat[dat$health == 0, ])^(-1/2), length.out = 30),
  xi.vector1 = seq(0, nrow(dat[dat$health == 1, ])^(-1/2), length.out = 30)
)

cat("Xi range for health=0: [", round(min(xi.vector.list[[1]]), 4), ",",
    round(max(xi.vector.list[[1]]), 4), "]\n")
cat("Xi range for health=1: [", round(min(xi.vector.list[[2]]), 4), ",",
    round(max(xi.vector.list[[2]]), 4), "]\n\n")

# Check if all sensitivity analysis results already exist
all_sa_files_exist <- TRUE
for (k in 1:2) {
  xi.vector <- xi.vector.list[[k]]
  for (xi in xi.vector) {
    sa_file <- paste0("MHD_results/EBMR_IPW_", subset_names[k], "_123_mild",
                      round(xi, 3), "_OR_test10.RDS")
    if (!file.exists(sa_file)) {
      all_sa_files_exist <- FALSE
      break
    }
  }
  if (!all_sa_files_exist) break
}

if (all_sa_files_exist) {
  cat("=== Sensitivity analysis results already exist, skipping computation ===\n")
} else {
  cat("Running sensitivity analysis for 30 xi values...\n")

  for (k in 1:2) {
    subdat <- dat[dat$health == k - 1, ]
    xi.vector <- xi.vector.list[[k]]

    for (xi in xi.vector) {
      ps_specifications <- list(
        formula.list = full_ps_specifications$formula.list[1:3],
        h_x_names.list = full_ps_specifications$h_x_names.list[1:3],
        alpha_init.list = init.list[[k]][1:3],
        outcome = full_ps_specifications$outcome,
        inv_link = full_ps_specifications$inv_link
      )

      ebmr <- EBMRAlgorithmFast2$new("teacher_report", ps_specifications, subdat, W)

      h_nu <- function(data) cbind(f = data$father, p = data$parent_report, fp = data$fp)
      result <- ebmr$EBMR_IPW(h_nu = h_nu, true_ps = NULL)

      # Exponential tilt function: exp(xi * (y + father + parent_report))
      exp_tilt <- function(y, x) exp(xi * as.matrix(cbind(y, x)) %*% c(1, 1, 1))

      sa_result <- ebmr$EBMR_IPW_with_locally_misspecified_model(
        ps.matrix = result$ps.matrix,
        perturb_ps = possibly_true_ps.vector[k],
        exp_tilt = exp_tilt,
        exp_tilt_x_names = exp_tilt_x_names,
        h_nu = h_nu
      )

      save_file <- paste0("MHD_results/EBMR_IPW_", subset_names[k], "_123_mild",
                          round(xi, 3), "_OR_test10.RDS")
      saveRDS(sa_result, save_file)
    }
    cat("  Completed health =", k - 1, "\n")
  }
}

# Read and summarize sensitivity analysis results
print_section("Sensitivity Analysis Results")

sa_results <- matrix(NA, 30, 2)
for (i in 1:30) {
  all_estimates <- matrix(NA, length(health), 2)
  for (k in 1:length(health)) {
    xi <- xi.vector.list[[k]][i]
    read_file <- paste0("MHD_results/EBMR_IPW_", subset_names[k], "_123_mild",
                        round(xi, 3), "_OR_test10.RDS")
    sa_result <- readRDS(read_file)
    all_estimates[k, ] <- unlist(sa_result[1:2])
  }
  effect <- effect.f(all_estimates[1, 1], all_estimates[2, 1])
  Sigma <- matrix(c(all_estimates[1, 2]^2, 0, 0, all_estimates[2, 2]^2), 2, 2)
  se <- sqrt(delta_method(all_estimates[1, 1], all_estimates[2, 1]) %*% Sigma %*%
             delta_method(all_estimates[1, 1], all_estimates[2, 1]))
  sa_results[i, ] <- c(effect, se)
}

cat("\n  Summary of theta estimates across xi values:\n")
cat(sprintf("  %-25s %10.4f\n", "Minimum theta:", min(sa_results[, 1])))
cat(sprintf("  %-25s %10.4f\n", "Maximum theta:", max(sa_results[, 1])))
cat(sprintf("  %-25s %10.4f\n", "Range:", max(sa_results[, 1]) - min(sa_results[, 1])))
cat(sprintf("  %-25s %10.4f\n", "Theta at xi=0:", sa_results[1, 1]))
cat(sprintf("  %-25s %10.4f\n", "Theta at max xi:", sa_results[30, 1]))

# Display results table
print_section("Sensitivity Analysis: Theta by Xi Index")
cat(sprintf("\n  %5s  %12s  %12s  %24s\n", "Index", "Theta", "SE", "95% CI"))
cat("  ", paste(rep("-", 60), collapse = ""), "\n")

# Show first 5, middle, and last 5 values
show_indices <- c(1:5, 15, 26:30)
for (i in show_indices) {
  ci <- format_ci(sa_results[i, 1] - 1.96 * sa_results[i, 2],
                  sa_results[i, 1] + 1.96 * sa_results[i, 2])
  cat(sprintf("  %5d  %12.4f  %12.4f  %24s\n",
              i, sa_results[i, 1], sa_results[i, 2], ci))
  if (i == 5) cat("  %5s  %12s  %12s  %24s\n", "...", "...", "...", "...")
  if (i == 15) cat("  %5s  %12s  %12s  %24s\n", "...", "...", "...", "...")
}

cat("  ", paste(rep("-", 60), collapse = ""), "\n")

#------------------------------------------------------------------------------#
# Sensitivity Analysis Visualization
#------------------------------------------------------------------------------#
cat("\n  Generating sensitivity analysis visualization...\n")

# Create xi values for x-axis (use health=0 as reference, both are similar)
xi_values <- xi.vector.list[[1]]

# Compute 95% CI from baseline theta_{111}^{111} (the full model)
# Load the full model results to get the correct SE
full_res0 <- readRDS("MHD_results/EBMR_IPW_health0_123_W_optimal_OR_test10.RDS")
full_res1 <- readRDS("MHD_results/EBMR_IPW_health1_123_W_optimal_OR_test10.RDS")

# Extract probabilities and SEs
p0_full <- full_res0$mu_ipw
p1_full <- full_res1$mu_ipw
se_p0_full <- full_res0$se_ipw
se_p1_full <- full_res1$se_ipw

# Compute baseline theta
baseline_theta <- sa_results[1, 1]

# Compute SE using delta method
# d(theta)/d(p0) = -1/(p0*(1-p0)), d(theta)/d(p1) = 1/(p1*(1-p1))
grad_p0 <- -1/(p0_full * (1 - p0_full))
grad_p1 <- 1/(p1_full * (1 - p1_full))
baseline_se <- sqrt(grad_p0^2 * se_p0_full^2 + grad_p1^2 * se_p1_full^2)

# 95% CI from baseline estimate (constant band across all xi)
baseline_ci_lower <- baseline_theta - 1.96 * baseline_se
baseline_ci_upper <- baseline_theta + 1.96 * baseline_se

cat(sprintf("  Baseline theta: %.4f, SE: %.4f\n", baseline_theta, baseline_se))
cat(sprintf("  95%% CI: [%.4f, %.4f]\n", baseline_ci_lower, baseline_ci_upper))

# Check which perturbed estimates lie within the 95% CI
in_ci <- sa_results[, 1] >= baseline_ci_lower & sa_results[, 1] <= baseline_ci_upper
n_in_ci <- sum(in_ci)
cat(sprintf("  Perturbed estimates within 95%% CI: %d / %d (%.1f%%)\n",
            n_in_ci, length(in_ci), 100 * n_in_ci / length(in_ci)))

# Save sensitivity analysis plot
png("MHD_results/sensitivity_analysis.png", width = 10, height = 7, units = "in", res = 300)

# Set up margins
par(mar = c(5, 5, 4, 2))

# Compute y-axis limits focused on the actual theta range (zoomed in)
# Add small padding based on the data range
theta_range <- max(sa_results[, 1]) - min(sa_results[, 1])
y_padding <- max(theta_range * 0.15, 0.005)  # At least 0.005 padding
y_min <- min(sa_results[, 1]) - y_padding
y_max <- max(sa_results[, 1]) + y_padding

# Create the plot with zoomed y-axis
n_xi <- length(xi_values)
plot(xi_values, sa_results[, 1], type = "n",
     xlim = c(0, max(xi_values) * 1.02),
     ylim = c(y_min, y_max),
     xlab = expression(paste("Sensitivity Parameter (", xi, ")")),
     ylab = expression(paste("Log Odds Ratio (", hat(theta), ")")),
     main = "",
     las = 1,
     cex.lab = 1.3,
     cex.axis = 1.1)

# Add grid
grid(col = "gray85", lty = 1, lwd = 0.5)

# Add title
mtext("Sensitivity Analysis: Effect of Local Model Misspecification",
      side = 3, line = 2.5, cex = 1.5, font = 2)
mtext("EBMR estimator under exponential tilt perturbation",
      side = 3, line = 1, cex = 1.0, font = 3)

# Add main estimate line
lines(xi_values, sa_results[, 1], col = "#2874A6", lwd = 3)
points(xi_values, sa_results[, 1], pch = 19, col = "#2874A6", cex = 1.0)

# Add reference lines
# Baseline (xi = 0)
abline(h = sa_results[1, 1], lty = 3, col = "#666666", lwd = 1.5)

# Mark key points
# Starting point (xi = 0)
points(xi_values[1], sa_results[1, 1], pch = 21, bg = "#F39C12", col = "black", cex = 2.0, lwd = 2)

# End point (max xi)
points(xi_values[n_xi], sa_results[n_xi, 1], pch = 21, bg = "#E74C3C", col = "black", cex = 2.0, lwd = 2)

# Add annotations for key values
text(xi_values[1] + max(xi_values) * 0.08, sa_results[1, 1],
     sprintf("%.4f", sa_results[1, 1]),
     cex = 0.9, col = "#F39C12", font = 2)

text(xi_values[n_xi] - max(xi_values) * 0.08, sa_results[n_xi, 1],
     sprintf("%.4f", sa_results[n_xi, 1]),
     cex = 0.9, col = "#E74C3C", font = 2)

# Legend (simplified)
legend("topright",
       legend = c(
         expression(paste("EBMR ", hat(theta), " (", xi, ")")),
         expression(paste("Baseline (", xi, " = 0)"))
       ),
       lty = c(1, 3),
       lwd = c(3, 1.5),
       pch = c(19, NA),
       col = c("#2874A6", "#666666"),
       pt.cex = c(1.0, NA),
       bg = "white",
       box.lwd = 1,
       cex = 0.95)

# Add interpretation text
par(xpd = TRUE)
text_x <- max(xi_values) * 0.5
text_y_base <- y_min - y_padding * 0.3
text(text_x, text_y_base,
     sprintf("Range of theta: [%.4f, %.4f], Change: %.4f",
             min(sa_results[, 1]), max(sa_results[, 1]), theta_range),
     cex = 0.85, font = 1)
par(xpd = FALSE)

dev.off()

cat("  Sensitivity analysis plot saved to: MHD_results/sensitivity_analysis.png\n")

#------------------------------------------------------------------------------#
# Additional: Sensitivity Analysis with separate p0 and p1 trajectories
#------------------------------------------------------------------------------#
cat("  Generating detailed sensitivity plot with separate health groups...\n")

# Collect p0 and p1 estimates at each xi
p0_estimates <- numeric(30)
p1_estimates <- numeric(30)
p0_se <- numeric(30)
p1_se <- numeric(30)

for (i in 1:30) {
  for (k in 1:length(health)) {
    xi <- xi.vector.list[[k]][i]
    read_file <- paste0("MHD_results/EBMR_IPW_", subset_names[k], "_123_mild",
                        round(xi, 3), "_OR_test10.RDS")
    sa_result <- readRDS(read_file)
    if (k == 1) {
      p0_estimates[i] <- sa_result$mu_ipw
      p0_se[i] <- sa_result$se_ipw
    } else {
      p1_estimates[i] <- sa_result$mu_ipw
      p1_se[i] <- sa_result$se_ipw
    }
  }
}

# Save detailed sensitivity plot
png("MHD_results/sensitivity_analysis_detailed.png", width = 12, height = 5, units = "in", res = 300)

par(mfrow = c(1, 3), mar = c(5, 5, 4, 2))

# Panel 1: p0 (health = 0) estimates
p0_has_se <- !any(is.na(p0_se))
plot(xi_values, p0_estimates, type = "n",
     xlab = expression(xi),
     ylab = expression(paste(hat(p)[0], " (Health = 0)")),
     main = "Unexposed Group (Health = 0)",
     las = 1, cex.lab = 1.2)
if (p0_has_se) {
  polygon(c(xi_values, rev(xi_values)),
          c(p0_estimates - 1.96 * p0_se, rev(p0_estimates + 1.96 * p0_se)),
          col = rgb(0.2, 0.6, 0.9, 0.2), border = NA)
}
lines(xi_values, p0_estimates, col = "#3498DB", lwd = 2)
points(xi_values, p0_estimates, pch = 19, col = "#3498DB", cex = 0.8)
abline(h = p0_estimates[1], lty = 2, col = "gray50")
grid(col = "gray90")

# Panel 2: p1 (health = 1) estimates
p1_has_se <- !any(is.na(p1_se))
plot(xi_values, p1_estimates, type = "n",
     xlab = expression(xi),
     ylab = expression(paste(hat(p)[1], " (Health = 1)")),
     main = "Exposed Group (Health = 1)",
     las = 1, cex.lab = 1.2)
if (p1_has_se) {
  polygon(c(xi_values, rev(xi_values)),
          c(p1_estimates - 1.96 * p1_se, rev(p1_estimates + 1.96 * p1_se)),
          col = rgb(0.9, 0.3, 0.2, 0.2), border = NA)
}
lines(xi_values, p1_estimates, col = "#E74C3C", lwd = 2)
points(xi_values, p1_estimates, pch = 19, col = "#E74C3C", cex = 0.8)
abline(h = p1_estimates[1], lty = 2, col = "gray50")
grid(col = "gray90")

# Panel 3: Log Odds Ratio (theta)
plot(xi_values, sa_results[, 1], type = "n",
     xlab = expression(xi),
     ylab = expression(hat(theta)),
     main = expression(paste("Log Odds Ratio (", hat(theta), ")")),
     las = 1, cex.lab = 1.2)
if (has_se) {
  polygon(c(xi_values, rev(xi_values)),
          c(sa_ci_lower, rev(sa_ci_upper)),
          col = rgb(0.2, 0.2, 0.3, 0.2), border = NA)
}
lines(xi_values, sa_results[, 1], col = "#2C3E50", lwd = 2)
points(xi_values, sa_results[, 1], pch = 19, col = "#2C3E50", cex = 0.8)
abline(h = 0, lty = 1, col = "#E74C3C", lwd = 1.5)
abline(h = sa_results[1, 1], lty = 2, col = "gray50")
grid(col = "gray90")

# Add overall title
mtext("Sensitivity Analysis: Effect of Exponential Tilt Perturbation",
      side = 3, outer = TRUE, line = -1.5, cex = 1.3, font = 2)

dev.off()

cat("  Detailed sensitivity plot saved to: MHD_results/sensitivity_analysis_detailed.png\n")

print_header("ANALYSIS COMPLETE")
