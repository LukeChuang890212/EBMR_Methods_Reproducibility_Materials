#------------------------------------------------------------------------------#
# Population Mean Analysis using EBMRalgorithmFast
# This script estimates the overall population mean of teacher_report (y)
# using EBMR with multiple propensity score models
#------------------------------------------------------------------------------#

library(EBMRalgorithmFast)
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
dat$y <- dat$teacher_report

#------------------------------------------------------------------------------#
# Primary analysis - Complete Case
#------------------------------------------------------------------------------#
n <- nrow(dat)
missing_rate <- 1 - mean(dat$r)
mu_cc <- mean(dat$teacher_report[dat$r == 1])
se_cc <- sd(dat$teacher_report[dat$r == 1]) / sqrt(sum(dat$r))

cat("================================================================================\n")
cat("                    POPULATION MEAN ESTIMATION ANALYSIS                         \n")
cat("================================================================================\n\n")

cat("Sample size:", n, "\n")
cat("Missing rate:", round(missing_rate * 100, 1), "%\n")
cat("Complete case mean:", round(mu_cc, 4), "\n")
cat("Complete case SE:", round(se_cc, 4), "\n\n")

#------------------------------------------------------------------------------#
# Basic Setup for EBMR analysis
# Models:
#   pi1: r ~ o(y) + health + father
#   pi2: r ~ o(y) + health + parent_report
#   pi3: r ~ o(y) + father + parent_report
#------------------------------------------------------------------------------#
full_ps_specifications <- list(
  formula.list = list(
    r ~ o(teacher_report) + health + father,          # Model 1: MNAR with health, father
    r ~ o(teacher_report) + health + parent_report,   # Model 2: MNAR with health, parent_report
    r ~ o(teacher_report) + father + parent_report    # Model 3: MNAR with father, parent_report
  ),
  h_x_names.list = list(
    c("health", "father", "parent_report"),
    c("health", "father", "parent_report"),
    c("health", "father", "parent_report")
  ),
  inv_link = function(eta) 1 / (1 + exp(eta))
)

# Optimal weighting matrix
W <- function(g.matrix) {
  solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
}

# Initial values for optimization (will be refined by glm)
alpha_init.list <- list(
  c(0, 1, 0, 0),      # Model 1: intercept, o(y), health, father
  c(0, 1, 0, 0),      # Model 2: intercept, o(y), health, parent_report
  c(0, 1, 0, 0)       # Model 3: intercept, o(y), father, parent_report
)

# h_nu function for ensemble step
h_nu <- function(data) cbind(
  health = data$health,
  father = data$father,
  parent_report = data$parent_report,
  hf = data$health * data$father
)

#------------------------------------------------------------------------------#
# Generate all model combinations
#------------------------------------------------------------------------------#
all_model_sets <- list()
for (model_num in 1:3) {
  model_combinations <- combn(3, model_num)
  for (i in 1:ncol(model_combinations)) {
    all_model_sets <- c(all_model_sets, list(model_combinations[, i]))
  }
}
# all_model_sets: [[1]] = 1, [[2]] = 2, [[3]] = 3, [[4]] = c(1,2), [[5]] = c(1,3), [[6]] = c(2,3), [[7]] = c(1,2,3)

# Create binary labels for model sets
model_set_to_binary <- function(model_set) {
  paste0(as.integer(1 %in% model_set),
         as.integer(2 %in% model_set),
         as.integer(3 %in% model_set))
}

model_set_labels <- sapply(all_model_sets, model_set_to_binary)
# "100", "010", "001", "110", "101", "011", "111"

#------------------------------------------------------------------------------#
# Run EBMR estimation for all model combinations
#------------------------------------------------------------------------------#
results_dir <- "MHD_results"
if (!dir.exists(results_dir)) dir.create(results_dir)

# Check if all EBMR results already exist
all_ebmr_files_exist <- TRUE
for (j in 1:length(all_model_sets)) {
  save_file <- paste0(results_dir, "/popmean_EBMR_IPW_", model_set_labels[j], ".RDS")
  if (!file.exists(save_file)) {
    all_ebmr_files_exist <- FALSE
    break
  }
}

# Store results
results_list <- list()
ps_fit.list_full <- NULL
nu_hat_full <- NULL
w_hat_full <- NULL

if (all_ebmr_files_exist) {
  cat("=== EBMR results already exist, loading from files ===\n\n")

  for (j in 1:length(all_model_sets)) {
    save_file <- paste0(results_dir, "/popmean_EBMR_IPW_", model_set_labels[j], ".RDS")
    results_list[[model_set_labels[j]]] <- readRDS(save_file)
  }

  # Load full model to get ps_fit.list for summary
  full_result <- results_list[["111"]]
  nu_hat_full <- full_result$nu.hat
  w_hat_full <- full_result$w.hat

  # Re-fit to get ps_fit.list for summary
  ps_specifications <- list(
    formula.list = full_ps_specifications$formula.list[1:3],
    h_x_names.list = full_ps_specifications$h_x_names.list[1:3],
    alpha_init.list = alpha_init.list[1:3],
    inv_link = full_ps_specifications$inv_link
  )
  ebmr <- EBMRAlgorithmFast$new("teacher_report", ps_specifications, dat, W)
  ps_fit.list_full <- ebmr$ps_fit.list

} else {
  cat("=== Running EBMR estimation for all model combinations ===\n\n")

  for (j in 1:length(all_model_sets)) {
    model_set <- all_model_sets[[j]]
    label <- model_set_labels[j]
    cat("Processing model set:", label, "(models:", paste(model_set, collapse = ","), ")\n")

    ps_specifications <- list(
      formula.list = full_ps_specifications$formula.list[model_set],
      h_x_names.list = full_ps_specifications$h_x_names.list[model_set],
      alpha_init.list = alpha_init.list[model_set],
      inv_link = full_ps_specifications$inv_link
    )

    ebmr <- EBMRAlgorithmFast$new("teacher_report", ps_specifications, dat, W)
    result <- ebmr$EBMR_IPW(h_nu = h_nu, true_ps = NULL)

    save_file <- paste0(results_dir, "/popmean_EBMR_IPW_", label, ".RDS")
    saveRDS(result, save_file)
    results_list[[label]] <- result

    if (j == 7) {  # Full model (1,2,3)
      ps_fit.list_full <- ebmr$ps_fit.list
      nu_hat_full <- result$nu.hat
      w_hat_full <- result$w.hat
    }
  }
}

# Update alpha_init.list with fitted values for bootstrap
for (j in 1:length(all_model_sets)) {
  model_set <- all_model_sets[[j]]
  if (j == 7) {
    for (m in 1:3) {
      alpha_init.list[[m]] <- ps_fit.list_full[[m]]$coefficients
    }
  }
}

#------------------------------------------------------------------------------#
# Propensity Score Model Summaries
#------------------------------------------------------------------------------#
cat("\n")
cat("================================================================================\n")
cat("                       PROPENSITY SCORE MODEL SUMMARIES                         \n")
cat("================================================================================\n")

# Model descriptions
model_names <- c(
  "Model 1 (pi1): r ~ o(y) + health + father",
  "Model 2 (pi2): r ~ o(y) + health + parent_report",
  "Model 3 (pi3): r ~ o(y) + father + parent_report"
)

for (i in 1:3) {
  ps_fit <- ps_fit.list_full[[i]]
  cat("\n--------------------------------------------------------------------------------\n")
  cat(" ", model_names[i], "\n")
  cat("--------------------------------------------------------------------------------\n\n")

  # Get coefficients and SEs
  coef <- ps_fit$coefficients
  se <- ps_fit$se
  z_val <- coef / se
  p_val <- (1 - pnorm(abs(z_val))) * 2

  # Determine variable names based on model
  if (i == 1) {
    var_names <- c("Intercept", "o(teacher_report)", "health", "father")
  } else if (i == 2) {
    var_names <- c("Intercept", "o(teacher_report)", "health", "parent_report")
  } else {
    var_names <- c("Intercept", "o(teacher_report)", "father", "parent_report")
  }

  # Print table header
  cat(sprintf("  %-20s %12s %12s %12s %12s\n", "Variable", "Estimate", "Std.Error", "z value", "Pr(>|z|)"))
  cat("  ", paste(rep("-", 72), collapse = ""), "\n")

  # Print each coefficient
  for (k in 1:length(coef)) {
    stars <- ""
    if (p_val[k] < 0.001) stars <- "***"
    else if (p_val[k] < 0.01) stars <- "**"
    else if (p_val[k] < 0.05) stars <- "*"
    else if (p_val[k] < 0.1) stars <- "."

    cat(sprintf("  %-20s %12.4f %12.4f %12.3f %10.4f%s\n",
                var_names[k], coef[k], se[k], z_val[k], p_val[k], stars))
  }
  cat("  ", paste(rep("-", 72), collapse = ""), "\n")
}

cat("\n  Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")

#------------------------------------------------------------------------------#
# Nu and W estimates
#------------------------------------------------------------------------------#
cat("\n")
cat("================================================================================\n")
cat("                           NU AND W ESTIMATES                                   \n")
cat("================================================================================\n\n")

cat("Full model (all 3 PS models combined):\n\n")
cat("  nu.hat (ensemble parameters):\n")
for (i in 1:length(nu_hat_full)) {
  cat(sprintf("    Model %d: %10.4f\n", i, nu_hat_full[i]))
}

cat("\n  w.hat (model weights, sum to 1):\n")
for (i in 1:length(w_hat_full)) {
  cat(sprintf("    Model %d: %10.4f (%.1f%%)\n", i, w_hat_full[i], w_hat_full[i] * 100))
}

#------------------------------------------------------------------------------#
# Summary table: Point Estimates and SEs
#------------------------------------------------------------------------------#
cat("\n")
cat("================================================================================\n")
cat("                         ESTIMATION RESULTS                                     \n")
cat("================================================================================\n\n")

cat("  Population Mean of Teacher Report (mu)\n")
cat("  Notation: mu_xxx where each position indicates model inclusion (1=yes, 0=no)\n")
cat("           Position 1: Model 1 (pi1), Position 2: Model 2 (pi2), Position 3: Model 3 (pi3)\n\n")

# Create results table
results_table <- data.frame(
  Estimator = character(),
  Models = character(),
  PE = numeric(),
  SE = numeric(),
  stringsAsFactors = FALSE
)

for (j in 1:length(all_model_sets)) {
  label <- model_set_labels[j]
  result <- results_list[[label]]
  models_used <- paste("pi", all_model_sets[[j]], sep = "", collapse = ", ")

  results_table <- rbind(results_table, data.frame(
    Estimator = paste0("mu_", label),
    Models = models_used,
    PE = result$mu_ipw,
    SE = result$se_ipw,
    stringsAsFactors = FALSE
  ))
}

# Add complete case
results_table <- rbind(results_table, data.frame(
  Estimator = "Complete Case",
  Models = "N/A",
  PE = mu_cc,
  SE = se_cc,
  stringsAsFactors = FALSE
))

#------------------------------------------------------------------------------#
# Bootstrap for SE estimation
#------------------------------------------------------------------------------#
B <- 1000
boot <- TRUE

if (boot) {
  cat("\n=== Standard Bootstrap (B =", B, ") ===\n\n")

  # Check if all bootstrap files exist
  all_boot_files_exist <- TRUE
  for (j in 1:length(all_model_sets)) {
    boot_file <- paste0(results_dir, "/popmean_boot_", model_set_labels[j], "_B", B, ".RDS")
    if (!file.exists(boot_file)) {
      all_boot_files_exist <- FALSE
      break
    }
  }

  boot_se <- rep(NA, length(all_model_sets) + 1)  # +1 for complete case

  if (all_boot_files_exist) {
    cat("All bootstrap files already exist, loading results...\n")

    for (j in 1:length(all_model_sets)) {
      boot_file <- paste0(results_dir, "/popmean_boot_", model_set_labels[j], "_B", B, ".RDS")
      boot_result <- readRDS(boot_file)
      boot_se[j] <- sd(boot_result, na.rm = TRUE)
      cat("  Model set", model_set_labels[j], ": Boot SE =", round(boot_se[j], 4),
          ", Valid:", sum(!is.na(boot_result)), "/", B, "\n")
    }

    # Bootstrap for complete case
    boot_cc_file <- paste0(results_dir, "/popmean_boot_CC_B", B, ".RDS")
    if (file.exists(boot_cc_file)) {
      boot_cc_result <- readRDS(boot_cc_file)
      boot_se[8] <- sd(boot_cc_result, na.rm = TRUE)
    }

  } else {
    cat("Running bootstrap for each model combination...\n\n")

    # Setup parallel cluster
    cores <- parallel::detectCores()
    n_cores <- max(1, cores - 2)

    for (j in 1:length(all_model_sets)) {
      model_set <- all_model_sets[[j]]
      label <- model_set_labels[j]

      boot_file <- paste0(results_dir, "/popmean_boot_", label, "_B", B, ".RDS")

      if (file.exists(boot_file)) {
        cat("  Model set", label, "- already exists, skipping\n")
        boot_result <- readRDS(boot_file)
        boot_se[j] <- sd(boot_result, na.rm = TRUE)
      } else {
        cat("  Model set", label, "- running bootstrap...\n")

        # Get nu_init for this model set
        if (length(model_set) == 1) {
          nu_init <- 1
        } else if (length(model_set) == 2) {
          nu_init <- nu_hat_full[model_set]
        } else {
          nu_init <- nu_hat_full
        }

        cl <- parallel::makeCluster(n_cores)
        doSNOW::registerDoSNOW(cl)

        pb <- txtProgressBar(max = B, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
        parallel_packages <- c("EBMRalgorithmFast", "stringr", "Matrix", "numDeriv")

        ps_spec_boot <- list(
          formula.list = full_ps_specifications$formula.list[model_set],
          h_x_names.list = full_ps_specifications$h_x_names.list[model_set],
          alpha_init.list = alpha_init.list[model_set],
          inv_link = full_ps_specifications$inv_link
        )

        boot_result <- foreach::foreach(
          i = 1:B,
          .combine = 'c',
          .options.snow = opts,
          .packages = parallel_packages
        ) %dopar% {

          tryCatch({
            # Resample with replacement
            boot_indices <- sample.int(n, replace = TRUE)
            boot_dat <- dat[boot_indices, ]

            ebmr_boot <- EBMRAlgorithmFast$new("teacher_report", ps_spec_boot, boot_dat, W)
            result_boot <- ebmr_boot$EBMR_IPW(h_nu = h_nu, nu_init = nu_init,
                                               true_ps = NULL, se.fit = FALSE)
            result_boot$mu_ipw

          }, error = function(e) {
            NA
          })
        }

        close(pb)
        parallel::stopCluster(cl)
        gc()

        saveRDS(boot_result, boot_file)
        boot_se[j] <- sd(boot_result, na.rm = TRUE)
        cat("    Valid samples:", sum(!is.na(boot_result)), "/", B, "\n")
      }
    }

    # Bootstrap for complete case
    boot_cc_file <- paste0(results_dir, "/popmean_boot_CC_B", B, ".RDS")
    if (!file.exists(boot_cc_file)) {
      cat("  Complete Case - running bootstrap...\n")
      set.seed(20241115)
      boot_cc_result <- replicate(B, {
        boot_indices <- sample.int(n, replace = TRUE)
        boot_dat <- dat[boot_indices, ]
        mean(boot_dat$teacher_report[boot_dat$r == 1])
      })
      saveRDS(boot_cc_result, boot_cc_file)
      boot_se[8] <- sd(boot_cc_result, na.rm = TRUE)
    } else {
      boot_cc_result <- readRDS(boot_cc_file)
      boot_se[8] <- sd(boot_cc_result, na.rm = TRUE)
    }
  }

  # Add bootstrap SE to results table
  results_table$Boot_SE <- boot_se
}

#------------------------------------------------------------------------------#
# Print final results table
#------------------------------------------------------------------------------#
cat("\n")
cat("================================================================================\n")
cat("                           FINAL RESULTS TABLE                                  \n")
cat("================================================================================\n\n")

cat(sprintf("  %-15s %-25s %10s %10s %10s   %-22s\n",
            "Estimator", "Models Used", "PE", "SE", "Boot SE", "95% CI"))
cat("  ", paste(rep("-", 100), collapse = ""), "\n")

for (i in 1:nrow(results_table)) {
  ci_lower <- results_table$PE[i] - qnorm(0.975) * results_table$SE[i]
  ci_upper <- results_table$PE[i] + qnorm(0.975) * results_table$SE[i]
  ci_str <- sprintf("[%.4f, %.4f]", ci_lower, ci_upper)

  boot_se_str <- ifelse(is.na(results_table$Boot_SE[i]), "--",
                         sprintf("%.4f", results_table$Boot_SE[i]))

  cat(sprintf("  %-15s %-25s %10.4f %10.4f %10s   %-22s\n",
              results_table$Estimator[i],
              results_table$Models[i],
              results_table$PE[i],
              results_table$SE[i],
              boot_se_str,
              ci_str))
}
cat("  ", paste(rep("-", 100), collapse = ""), "\n")

# Add nu.hat and w.hat for each combination
cat("\n")
cat("================================================================================\n")
cat("                     NU AND W FOR EACH MODEL COMBINATION                        \n")
cat("================================================================================\n\n")

cat(sprintf("  %-15s %35s %35s\n", "Estimator", "nu.hat", "w.hat"))
cat("  ", paste(rep("-", 90), collapse = ""), "\n")

for (j in 1:length(all_model_sets)) {
  label <- model_set_labels[j]
  result <- results_list[[label]]

  nu_str <- paste(sprintf("%.4f", result$nu.hat), collapse = ", ")
  w_str <- paste(sprintf("%.4f", result$w.hat), collapse = ", ")

  cat(sprintf("  %-15s %35s %35s\n",
              paste0("mu_", label), nu_str, w_str))
}
cat("  ", paste(rep("-", 90), collapse = ""), "\n")

#------------------------------------------------------------------------------#
# Sensitivity Analysis
#------------------------------------------------------------------------------#
cat("\n")
cat("================================================================================\n")
cat("                          SENSITIVITY ANALYSIS                                  \n")
cat("================================================================================\n\n")

cat("Perturbation: Model 1 (pi1) with exponential tilt function\n")
cat("  exp(xi * (y + health + father + parent_report))\n")
cat("  for xi in 30 equally spaced values in [0, n^(-1/2)]\n\n")

xi_max <- n^(-1/2)
xi_values <- seq(0, xi_max, length.out = 30)

cat("Xi range: [0, ", round(xi_max, 6), "]\n\n")

# Check if all sensitivity analysis results exist
all_sa_files_exist <- TRUE
for (i in 1:30) {
  sa_file <- paste0(results_dir, "/popmean_sensitivity_xi", i, ".RDS")
  if (!file.exists(sa_file)) {
    all_sa_files_exist <- FALSE
    break
  }
}

sa_results <- matrix(NA, 30, 2)
colnames(sa_results) <- c("mu_ipw", "se_ipw")

if (all_sa_files_exist) {
  cat("=== Sensitivity analysis results already exist, loading from files ===\n")

  for (i in 1:30) {
    sa_file <- paste0(results_dir, "/popmean_sensitivity_xi", i, ".RDS")
    sa_result <- readRDS(sa_file)
    sa_results[i, ] <- c(sa_result$mu_ipw, sa_result$se_ipw)
  }

} else {
  cat("Running sensitivity analysis for 30 xi values...\n")

  # First fit the full model to get baseline
  ps_specifications <- list(
    formula.list = full_ps_specifications$formula.list[1:3],
    h_x_names.list = full_ps_specifications$h_x_names.list[1:3],
    alpha_init.list = alpha_init.list[1:3],
    inv_link = full_ps_specifications$inv_link
  )

  ebmr <- EBMRAlgorithmFast$new("teacher_report", ps_specifications, dat, W)
  baseline_result <- ebmr$EBMR_IPW(h_nu = h_nu, true_ps = NULL)

  # Exponential tilt variables
  exp_tilt_x_names <- c("health", "father", "parent_report")

  for (i in 1:30) {
    xi <- xi_values[i]
    cat("  xi[", i, "] =", round(xi, 6), "\n")

    # Exponential tilt function: exp(xi * (y + health + father + parent_report))
    exp_tilt <- function(y, x) {
      exp(xi * as.matrix(cbind(y, x)) %*% c(1, 1, 1, 1))
    }

    # Perturb Model 1
    sa_result <- ebmr$EBMR_IPW_with_locally_misspecified_model(
      ps.matrix = baseline_result$ps.matrix,
      perturb_ps = 1,  # Perturb Model 1
      exp_tilt = exp_tilt,
      exp_tilt_x_names = exp_tilt_x_names,
      h_nu = h_nu
    )

    sa_results[i, ] <- c(sa_result$mu_ipw, sa_result$se_ipw)

    sa_file <- paste0(results_dir, "/popmean_sensitivity_xi", i, ".RDS")
    saveRDS(sa_result, sa_file)
  }
}

# Display sensitivity analysis results
cat("\n")
cat("  Sensitivity Analysis Results:\n\n")
cat(sprintf("  %5s  %12s  %12s  %24s\n", "Index", "xi", "mu_hat", "95% CI"))
cat("  ", paste(rep("-", 60), collapse = ""), "\n")

# Show first 5, middle, and last 5 values
show_indices <- c(1:5, 15, 26:30)
for (i in show_indices) {
  ci_lower <- sa_results[i, 1] - 1.96 * sa_results[i, 2]
  ci_upper <- sa_results[i, 1] + 1.96 * sa_results[i, 2]
  ci_str <- sprintf("[%.4f, %.4f]", ci_lower, ci_upper)

  cat(sprintf("  %5d  %12.6f  %12.4f  %24s\n", i, xi_values[i], sa_results[i, 1], ci_str))
  if (i == 5 || i == 15) cat("  %5s  %12s  %12s  %24s\n", "...", "...", "...", "...")
}
cat("  ", paste(rep("-", 60), collapse = ""), "\n")

# Summary statistics
cat("\n  Summary:\n")
cat(sprintf("    Baseline mu (xi=0):  %.4f\n", sa_results[1, 1]))
cat(sprintf("    Maximum mu:          %.4f (at xi = %.6f)\n", max(sa_results[, 1]), xi_values[which.max(sa_results[, 1])]))
cat(sprintf("    Minimum mu:          %.4f (at xi = %.6f)\n", min(sa_results[, 1]), xi_values[which.min(sa_results[, 1])]))
cat(sprintf("    Range:               %.4f\n", max(sa_results[, 1]) - min(sa_results[, 1])))

#------------------------------------------------------------------------------#
# Sensitivity Analysis Visualization
#------------------------------------------------------------------------------#
cat("\n  Generating sensitivity analysis plot...\n")

png(paste0(results_dir, "/popmean_sensitivity_analysis.png"), width = 10, height = 7, units = "in", res = 300)

par(mar = c(5, 5, 4, 2))

# Compute y-axis limits
mu_range <- max(sa_results[, 1]) - min(sa_results[, 1])
y_padding <- max(mu_range * 0.15, 0.01)
y_min <- min(sa_results[, 1]) - y_padding
y_max <- max(sa_results[, 1]) + y_padding

# Create plot
plot(xi_values, sa_results[, 1], type = "n",
     xlim = c(0, max(xi_values) * 1.02),
     ylim = c(y_min, y_max),
     xlab = expression(paste("Sensitivity Parameter (", xi, ")")),
     ylab = expression(paste("Population Mean (", hat(mu), ")")),
     main = "",
     las = 1,
     cex.lab = 1.3,
     cex.axis = 1.1)

# Add grid
grid(col = "gray85", lty = 1, lwd = 0.5)

# Title
mtext("Sensitivity Analysis: Population Mean Estimation",
      side = 3, line = 2.5, cex = 1.5, font = 2)
mtext("EBMR estimator under exponential tilt perturbation of Model 1",
      side = 3, line = 1, cex = 1.0, font = 3)

# Add 95% CI band
ci_lower <- sa_results[, 1] - 1.96 * sa_results[, 2]
ci_upper <- sa_results[, 1] + 1.96 * sa_results[, 2]
polygon(c(xi_values, rev(xi_values)),
        c(ci_lower, rev(ci_upper)),
        col = rgb(0.16, 0.45, 0.65, 0.2), border = NA)

# Main estimate line
lines(xi_values, sa_results[, 1], col = "#2874A6", lwd = 3)
points(xi_values, sa_results[, 1], pch = 19, col = "#2874A6", cex = 1.0)

# Baseline reference
abline(h = sa_results[1, 1], lty = 3, col = "#666666", lwd = 1.5)

# Mark key points
points(xi_values[1], sa_results[1, 1], pch = 21, bg = "#F39C12", col = "black", cex = 2.0, lwd = 2)
points(xi_values[30], sa_results[30, 1], pch = 21, bg = "#E74C3C", col = "black", cex = 2.0, lwd = 2)

# Add complete case reference
abline(h = mu_cc, lty = 2, col = "#27AE60", lwd = 2)

# Annotations
text(xi_values[1] + max(xi_values) * 0.08, sa_results[1, 1],
     sprintf("%.4f", sa_results[1, 1]), cex = 0.9, col = "#F39C12", font = 2)
text(xi_values[30] - max(xi_values) * 0.08, sa_results[30, 1],
     sprintf("%.4f", sa_results[30, 1]), cex = 0.9, col = "#E74C3C", font = 2)

# Legend
legend("topright",
       legend = c(
         expression(paste("EBMR ", hat(mu), " (", xi, ")")),
         "95% CI",
         expression(paste("Baseline (", xi, " = 0)")),
         "Complete Case"
       ),
       lty = c(1, NA, 3, 2),
       lwd = c(3, NA, 1.5, 2),
       pch = c(19, 15, NA, NA),
       col = c("#2874A6", rgb(0.16, 0.45, 0.65, 0.5), "#666666", "#27AE60"),
       pt.cex = c(1.0, 2.5, NA, NA),
       bg = "white",
       box.lwd = 1,
       cex = 0.95)

dev.off()

cat("  Sensitivity analysis plot saved to:", paste0(results_dir, "/popmean_sensitivity_analysis.png"), "\n")

#------------------------------------------------------------------------------#
# Final Summary
#------------------------------------------------------------------------------#
cat("\n")
cat("================================================================================\n")
cat("                              ANALYSIS COMPLETE                                 \n")
cat("================================================================================\n\n")

cat("Key Findings:\n\n")

# Find the full model result
full_model_result <- results_list[["111"]]
cat(sprintf("  1. Complete Case Estimate:     %.4f (SE: %.4f)\n", mu_cc, se_cc))
cat(sprintf("  2. EBMR (all 3 models):        %.4f (SE: %.4f)\n",
            full_model_result$mu_ipw, full_model_result$se_ipw))
cat(sprintf("  3. Difference (EBMR - CC):     %.4f\n", full_model_result$mu_ipw - mu_cc))

cat("\n  Model Weights (full model):\n")
for (i in 1:3) {
  cat(sprintf("    pi%d: %.4f (%.1f%%)\n", i, w_hat_full[i], w_hat_full[i] * 100))
}

cat("\n  Sensitivity Analysis Range: [", round(min(sa_results[, 1]), 4), ", ",
    round(max(sa_results[, 1]), 4), "]\n")

cat("\n  Output files saved to:", results_dir, "\n")
cat("    - popmean_EBMR_IPW_*.RDS (estimation results)\n")
cat("    - popmean_boot_*.RDS (bootstrap results)\n")
cat("    - popmean_sensitivity_*.RDS (sensitivity analysis)\n")
cat("    - popmean_sensitivity_analysis.png (visualization)\n")
