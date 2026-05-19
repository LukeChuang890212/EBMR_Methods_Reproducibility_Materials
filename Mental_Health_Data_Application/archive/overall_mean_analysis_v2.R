#------------------------------------------------------------------------------#
# Overall Population Mean Estimation: E[teacher_report]
# Using EBMR with 3 propensity score models
#------------------------------------------------------------------------------#

library(EBMRalgorithmFast3)
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
cat("     OVERALL POPULATION MEAN ESTIMATION: E[teacher_report]                     \n")
cat("================================================================================\n\n")

cat("Sample size:", n, "\n")
cat("Missing rate:", round(missing_rate * 100, 1), "%\n")
cat("Complete case mean:", round(mu_cc, 4), "\n\n")

summary(glm(health~r, data = dat))
summary(glm(father~r, data = dat))
summary(glm(parent_report~r, data = dat))

#------------------------------------------------------------------------------#
# MAR Analysis: IPW with saturated propensity score model
#------------------------------------------------------------------------------#
cat("\n================================================================================\n")
cat("                     MAR ANALYSIS: IPW ESTIMATOR                               \n")
cat("================================================================================\n\n")

# Fit saturated logistic regression for propensity score
# R ~ intercept + father + health + parent_report + all 2-way + 3-way interactions
ps_mar_fit <- glm(r ~ father * health * parent_report, data = dat, family = binomial)

cat("Propensity Score Model (MAR assumption):\n")
cat("  Formula: R ~ father * health * parent_report (saturated model)\n\n")
print(summary(ps_mar_fit))

# Get fitted propensity scores
ps_mar <- fitted(ps_mar_fit)
cat("\nFitted propensity score range: [", round(min(ps_mar), 4), ", ",
    round(max(ps_mar), 4), "]\n", sep = "")

# Compute IPW estimator under MAR
# mu_MAR = (1/n) * sum(r * y / ps)
mu_mar <- sum(dat$r * dat$teacher_report / ps_mar) / n
cat("IPW estimate under MAR:", round(mu_mar, 4), "\n\n")

# Bootstrap for MAR IPW SE
cat("Running bootstrap for MAR IPW standard error...\n")
B_mar <- 1000
boot_mar_file <- "MHD_results/popmean_boot_MAR_B1000.RDS"

if (file.exists(boot_mar_file)) {
  cat("  Loading existing bootstrap results from file...\n")
  boot_mar <- readRDS(boot_mar_file)
  boot_se_mar <- sd(boot_mar, na.rm = TRUE)
} else {
  set.seed(12345)
  boot_mar <- numeric(B_mar)

  pb_mar <- txtProgressBar(min = 0, max = B_mar, style = 3)
  for (b in 1:B_mar) {
    idx <- sample(1:n, n, replace = TRUE)
    dat_b <- dat[idx, ]

    # Fit propensity score model on bootstrap sample
    ps_fit_b <- tryCatch({
      glm(r ~ father * health * parent_report, data = dat_b, family = binomial)
    }, error = function(e) NULL, warning = function(w) NULL)

    if (!is.null(ps_fit_b)) {
      ps_b <- fitted(ps_fit_b)
      # Avoid division by very small propensity scores
      ps_b <- pmax(ps_b, 0.01)
      boot_mar[b] <- sum(dat_b$r * dat_b$teacher_report / ps_b) / nrow(dat_b)
    } else {
      boot_mar[b] <- NA
    }

    setTxtProgressBar(pb_mar, b)
  }
  close(pb_mar)

  # Save bootstrap results
  saveRDS(boot_mar, boot_mar_file)
  cat("\n  Saved:", boot_mar_file, "\n")

  boot_se_mar <- sd(boot_mar, na.rm = TRUE)
}

# MAR results summary
ci_mar <- c(mu_mar - 1.96 * boot_se_mar, mu_mar + 1.96 * boot_se_mar)
cat("\nMAR IPW Results:\n")
cat(sprintf("  Estimate:     %.4f\n", mu_mar))
cat(sprintf("  Bootstrap SE: %.4f\n", boot_se_mar))
cat(sprintf("  95%% CI:       [%.4f, %.4f]\n", ci_mar[1], ci_mar[2]))
cat(sprintf("  Difference from CC: %.4f\n\n", mu_mar - mu_cc))

#------------------------------------------------------------------------------#
# EBMR Setup
#------------------------------------------------------------------------------#
W <- function(g.matrix) {
  solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
}

# 3 models:
# pi_1: r ~ y + health + father + y:health + y:father
# pi_2: r ~ y + health + parent_report + y:health + y:parent_report
# pi_3: r ~ y + parent_report + father + y:parent_report + y:father
ps_specifications <- list(
  formula.list = list(
    r ~ teacher_report + health + father + teacher_report:health + teacher_report:father,
    r ~ teacher_report + health + parent_report + teacher_report:health + teacher_report:parent_report,
    r ~ teacher_report + parent_report + father + teacher_report:parent_report+teacher_report:father
  ),
  h_x_names.list = list(
    c("health", "father", "parent_report", "fp", "fh", "hp"),
    c("health", "father", "parent_report", "fp", "fh", "hp"),
    c("health", "father", "parent_report", "fp", "fh", "hp")
  ),
  outcome = "teacher_report",  # New in EBMRalgorithmFast3: specify outcome variable
  inv_link = function(eta) 1 / (1 + exp(eta))
)

# h_nu function for ensemble step
h_nu <- function(data) {
  cbind(
    health = data$health,
    father = data$father,
    parent_report = data$parent_report,
    fp = data$fp,
    fh = data$fh,
    hp = data$hp,
    fhp = data$fhp
  )
}

#------------------------------------------------------------------------------#
# Fit EBMR
#------------------------------------------------------------------------------#
ebmr <- EBMRAlgorithmFast3$new("teacher_report", ps_specifications, dat, W, method = "L-BFGS-B-prog")
result <- ebmr$EBMR_IPW(h_nu = h_nu, true_ps = NULL, method = "L-BFGS-B-prog")

#------------------------------------------------------------------------------#
# Model Summaries
#------------------------------------------------------------------------------#
cat("================================================================================\n")
cat("                       PROPENSITY SCORE MODEL SUMMARIES                        \n")
cat("================================================================================\n\n")

for (i in seq_along(ebmr$ps_fit.list)) {
  ps_fit <- ebmr$ps_fit.list[[i]]

  # Get formula as string for display
  formula_str <- deparse(ps_specifications$formula.list[[i]])
  cat(sprintf("Model %d: %s\n", i, formula_str))
  cat(sprintf("  Conv: %d, Obj: %.2e\n",
              ps_fit$gmm_fit$opt$convergence, ps_fit$gmm_fit$opt$value))

  coef <- ps_fit$coefficients
  se <- ps_fit$se
  z_val <- coef / se
  p_val <- 2 * pnorm(-abs(z_val))

  # Get variable names from design matrix column names
  var_names <- colnames(ps_fit$design_matrix)

  # Determine max width for variable names
  max_width <- max(15, max(nchar(var_names)) + 2)

  cat(sprintf("  %*s %10s %10s %10s %10s\n",
              max_width, "Variable", "Estimate", "SE", "z", "p"))
  for (k in seq_along(coef)) {
    stars <- ""
    if (p_val[k] < 0.001) stars <- "***"
    else if (p_val[k] < 0.01) stars <- "**"
    else if (p_val[k] < 0.05) stars <- "*"
    else if (p_val[k] < 0.1) stars <- "."
    cat(sprintf("  %*s %10.4f %10.4f %10.3f %10.4f%s\n",
                max_width, var_names[k], coef[k], se[k], z_val[k], p_val[k], stars))
  }

  # Report range of fitted propensity scores
  fitted_ps <- ps_fit$fitted.values
  cat(sprintf("  Fitted PS range: [%.4f, %.4f]\n", min(fitted_ps), max(fitted_ps)))
  cat("\n")
}

#------------------------------------------------------------------------------#
# Nu and W Estimates
#------------------------------------------------------------------------------#
cat("================================================================================\n")
cat("                           NU AND W ESTIMATES                                  \n")
cat("================================================================================\n\n")

cat("nu.hat:", round(result$nu.hat, 4), "\n")
cat("w.hat: ", round(result$w.hat, 4), "\n\n")

for (i in 1:3) {
  cat(sprintf("  Model %d weight: %.1f%%\n", i, result$w.hat[i] * 100))
}

#------------------------------------------------------------------------------#
# Population Mean Estimation Result
#------------------------------------------------------------------------------#
cat("\n================================================================================\n")
cat("                     POPULATION MEAN ESTIMATION RESULT                         \n")
cat("================================================================================\n\n")

cat(sprintf("  %-20s %12s %12s %20s\n", "Estimator", "Estimate", "SE", "95% CI"))
cat("  ", paste(rep("-", 70), collapse = ""), "\n")

ci_cc <- c(mu_cc - 1.96 * se_cc, mu_cc + 1.96 * se_cc)
cat(sprintf("  %-20s %12.4f %12.4f     [%.4f, %.4f]\n",
            "Complete Case", mu_cc, se_cc, ci_cc[1], ci_cc[2]))

ci_ebmr <- c(result$mu_ipw - 1.96 * result$se_ipw,
             result$mu_ipw + 1.96 * result$se_ipw)
cat(sprintf("  %-20s %12.4f %12.4f     [%.4f, %.4f]\n",
            "EBMR (3 models)", result$mu_ipw, result$se_ipw,
            ci_ebmr[1], ci_ebmr[2]))

cat("  ", paste(rep("-", 70), collapse = ""), "\n")
cat("\n  Difference (EBMR - CC):", round(result$mu_ipw - mu_cc, 4), "\n")

#------------------------------------------------------------------------------#
# All Model Combinations
#------------------------------------------------------------------------------#
cat("\n================================================================================\n")
cat("                     ALL MODEL COMBINATIONS                                    \n")
cat("================================================================================\n\n")

# Define all combinations
all_model_sets <- list(
  c(1),       # Model 1 only
  c(2),       # Model 2 only
  c(3),       # Model 3 only
  c(1, 2),    # Models 1 and 2
  c(1, 3),    # Models 1 and 3
  c(2, 3),    # Models 2 and 3
  c(1, 2, 3)  # All 3 models
)
model_set_labels <- c("100", "010", "001", "110", "101", "011", "111")
model_set_names <- c(
  "pi1 only",
  "pi2 only",
  "pi3 only",
  "pi1 + pi2",
  "pi1 + pi3",
  "pi2 + pi3",
  "pi1 + pi2 + pi3"
)

results_list <- list()

for (j in seq_along(all_model_sets)) {
  model_set <- all_model_sets[[j]]
  label <- model_set_labels[j]

  ps_spec_subset <- list(
    formula.list = ps_specifications$formula.list[model_set],
    h_x_names.list = ps_specifications$h_x_names.list[model_set],
    outcome = ps_specifications$outcome,  # Include outcome for EBMRalgorithmFast3
    inv_link = ps_specifications$inv_link
  )

  ebmr_sub <- EBMRAlgorithmFast3$new("teacher_report", ps_spec_subset, dat, W, method = "L-BFGS-B-prog")
  result_sub <- ebmr_sub$EBMR_IPW(h_nu = h_nu, true_ps = NULL, method = "L-BFGS-B-prog")

  results_list[[label]] <- list(
    mu = result_sub$mu_ipw,
    se = result_sub$se_ipw,
    nu = result_sub$nu.hat,
    w = result_sub$w.hat
  )
}

# Print results table
cat(sprintf("  %-20s %-18s %10s %10s %20s\n",
            "Label", "Models", "Estimate", "SE", "95% CI"))
cat("  ", paste(rep("-", 82), collapse = ""), "\n")

# Complete case first
cat(sprintf("  %-20s %-18s %10.4f %10.4f     [%.4f, %.4f]\n",
            "CC", "Complete Case", mu_cc, se_cc,
            mu_cc - 1.96 * se_cc, mu_cc + 1.96 * se_cc))

for (j in seq_along(all_model_sets)) {
  label <- model_set_labels[j]
  name <- model_set_names[j]
  res <- results_list[[label]]
  ci <- c(res$mu - 1.96 * res$se, res$mu + 1.96 * res$se)
  cat(sprintf("  %-20s %-18s %10.4f %10.4f     [%.4f, %.4f]\n",
              label, name, res$mu, res$se, ci[1], ci[2]))
}

cat("  ", paste(rep("-", 82), collapse = ""), "\n")

# Nu and W for each combination
cat("\n  Model Weights (w.hat) for each combination:\n\n")
cat(sprintf("  %-10s %-18s %30s\n", "Label", "Models", "w.hat"))
cat("  ", paste(rep("-", 60), collapse = ""), "\n")

for (j in seq_along(all_model_sets)) {
  label <- model_set_labels[j]
  name <- model_set_names[j]
  res <- results_list[[label]]
  w_str <- paste(sprintf("%.4f", res$w), collapse = ", ")
  cat(sprintf("  %-10s %-18s %30s\n", label, name, w_str))
}

#------------------------------------------------------------------------------#
# Bootstrap Standard Errors (Parallel with Progress Bar)
#------------------------------------------------------------------------------#
library(parallel)
library(foreach)
library(doSNOW)

B <- 1000
n_cores <- detectCores() - 1  # Leave one core free

# Check if standard bootstrap results already exist
boot_file_exists <- file.exists("MHD_results/popmean_boot_CC_B1000.RDS")

if (boot_file_exists) {
  cat("\n================================================================================\n")
  cat("                     BOOTSTRAP STANDARD ERRORS (FROM FILE)                     \n")
  cat("================================================================================\n\n")
  cat("  Bootstrap results already exist. Loading from files...\n\n")

  # Load existing results
  boot_cc <- readRDS("MHD_results/popmean_boot_CC_B1000.RDS")
  boot_se_cc <- sd(boot_cc, na.rm = TRUE)
  boot_se <- numeric(length(model_set_labels))
  names(boot_se) <- model_set_labels
  boot_results <- list()
  for (label in model_set_labels) {
    boot_file <- sprintf("MHD_results/popmean_boot_%s_B1000.RDS", label)
    if (file.exists(boot_file)) {
      boot_results[[label]] <- readRDS(boot_file)
      boot_se[label] <- sd(boot_results[[label]], na.rm = TRUE)
    } else {
      boot_se[label] <- NA
    }
  }
} else {
  cat("\n================================================================================\n")
  cat("                     BOOTSTRAP STANDARD ERRORS (PARALLEL)                      \n")
  cat("================================================================================\n\n")

  cat("Number of bootstrap samples:", B, "\n")
  cat("Number of cores:", n_cores, "\n\n")

  # Set up parallel backend with doSNOW for progress bar
  cl <- makeCluster(n_cores)
  registerDoSNOW(cl)

  # Progress bar
  pb <- txtProgressBar(max = B, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  cat("Running parallel bootstrap:\n")

  # Run bootstrap with foreach and progress bar
  boot_results_raw <- foreach(
    b = 1:B,
    .options.snow = opts,
    .packages = c("EBMRalgorithmFast3", "Matrix"),
    .export = c("dat", "n", "ps_specifications", "all_model_sets",
                "model_set_labels", "W", "h_nu")
  ) %dopar% {
    # Bootstrap sample
    set.seed(12345 + b)
    idx <- sample(1:n, n, replace = TRUE)
    dat_b <- dat[idx, ]

    # Complete case
    mu_cc_b <- mean(dat_b$teacher_report[dat_b$r == 1])

    # Each model combination
    mu_list <- list()
    for (j in seq_along(all_model_sets)) {
      model_set <- all_model_sets[[j]]
      label <- model_set_labels[j]

      ps_spec_b <- list(
        formula.list = ps_specifications$formula.list[model_set],
        h_x_names.list = ps_specifications$h_x_names.list[model_set],
        outcome = ps_specifications$outcome,
        inv_link = ps_specifications$inv_link
      )

      mu_list[[label]] <- tryCatch({
        ebmr_b <- EBMRAlgorithmFast3$new("teacher_report", ps_spec_b, dat_b, W, method = "L-BFGS-B-prog")
        result_b <- ebmr_b$EBMR_IPW(h_nu = h_nu, true_ps = NULL, se.fit = FALSE, method = "L-BFGS-B-prog")
        result_b$mu_ipw
      }, error = function(e) {
        NA
      })
    }

    c(CC = mu_cc_b, unlist(mu_list))
  }

  close(pb)
  stopCluster(cl)
  cat("\n")

  # Convert to matrix and extract results
  boot_matrix <- do.call(rbind, boot_results_raw)
  boot_cc <- boot_matrix[, "CC"]
  boot_results <- list()
  for (label in model_set_labels) {
    boot_results[[label]] <- boot_matrix[, label]
  }

  # Calculate bootstrap SEs
  boot_se <- sapply(boot_results, function(x) sd(x, na.rm = TRUE))
  boot_se_cc <- sd(boot_cc, na.rm = TRUE)

  # Save bootstrap results immediately
  cat("\n  Saving bootstrap results to MHD_results/...\n")
  saveRDS(boot_cc, "MHD_results/popmean_boot_CC_B1000.RDS")
  cat("    Saved: MHD_results/popmean_boot_CC_B1000.RDS\n")
  for (label in model_set_labels) {
    filename <- sprintf("MHD_results/popmean_boot_%s_B1000.RDS", label)
    saveRDS(boot_results[[label]], filename)
    cat(sprintf("    Saved: %s\n", filename))
  }
}

# Print bootstrap results table
cat(sprintf("\n  %-20s %-18s %10s %12s %12s\n",
            "Label", "Models", "Estimate", "Analytical", "Bootstrap"))
cat("  ", paste(rep("-", 76), collapse = ""), "\n")

cat(sprintf("  %-20s %-18s %10.4f %12.4f %12.4f\n",
            "CC", "Complete Case", mu_cc, se_cc, boot_se_cc))

for (j in seq_along(all_model_sets)) {
  label <- model_set_labels[j]
  name <- model_set_names[j]
  res <- results_list[[label]]
  cat(sprintf("  %-20s %-18s %10.4f %12.4f %12.4f\n",
              label, name, res$mu, res$se, boot_se[label]))
}

cat("  ", paste(rep("-", 76), collapse = ""), "\n")

#------------------------------------------------------------------------------#
# Perturbation Bootstrap Standard Errors (Parallel with Progress Bar)
#------------------------------------------------------------------------------#

# Check if perturbation bootstrap results already exist
perturb_file_exists <- file.exists("MHD_results/popmean_perturb_CC_B1000.RDS")

if (perturb_file_exists) {
  cat("\n================================================================================\n")
  cat("             PERTURBATION BOOTSTRAP STANDARD ERRORS (FROM FILE)                \n")
  cat("================================================================================\n\n")
  cat("  Perturbation bootstrap results already exist. Loading from files...\n\n")

  # Load existing results
  perturb_cc <- readRDS("MHD_results/popmean_perturb_CC_B1000.RDS")
  perturb_se_cc <- sd(perturb_cc, na.rm = TRUE)
  perturb_se <- numeric(length(model_set_labels))
  names(perturb_se) <- model_set_labels
  perturb_results <- list()
  for (label in model_set_labels) {
    perturb_file <- sprintf("MHD_results/popmean_perturb_%s_B1000.RDS", label)
    if (file.exists(perturb_file)) {
      perturb_results[[label]] <- readRDS(perturb_file)
      perturb_se[label] <- sd(perturb_results[[label]], na.rm = TRUE)
    } else {
      perturb_se[label] <- NA
    }
  }
} else {
  cat("\n================================================================================\n")
  cat("                PERTURBATION BOOTSTRAP STANDARD ERRORS (PARALLEL)              \n")
  cat("================================================================================\n\n")

  B_perturb <- 1000
  cat("Number of perturbation bootstrap samples:", B_perturb, "\n")
  cat("Number of cores:", n_cores, "\n\n")

  # Set up parallel backend with doSNOW for progress bar
  cl_perturb <- makeCluster(n_cores)
  registerDoSNOW(cl_perturb)

  # Progress bar
  pb_perturb <- txtProgressBar(max = B_perturb, style = 3)
  progress_perturb <- function(n) setTxtProgressBar(pb_perturb, n)
  opts_perturb <- list(progress = progress_perturb)

  cat("Running parallel perturbation bootstrap:\n")

  # Run perturbation bootstrap with foreach and progress bar
  perturb_results_raw <- foreach(
    b = 1:B_perturb,
    .options.snow = opts_perturb,
    .packages = c("EBMRalgorithmFast3", "Matrix"),
    .export = c("dat", "n", "ps_specifications", "all_model_sets",
                "model_set_labels", "W", "h_nu")
  ) %dopar% {
    # Generate exponential weights
    set.seed(12345 + b)
    wt <- rexp(n, rate = 1)

    # Complete case (weighted)
    r_idx <- dat$r == 1
    mu_cc_b <- sum(wt[r_idx] * dat$teacher_report[r_idx]) / sum(wt[r_idx])

    # Each model combination
    mu_list <- list()
    for (j in seq_along(all_model_sets)) {
      model_set <- all_model_sets[[j]]
      label <- model_set_labels[j]

      ps_spec_b <- list(
        formula.list = ps_specifications$formula.list[model_set],
        h_x_names.list = ps_specifications$h_x_names.list[model_set],
        outcome = ps_specifications$outcome,
        inv_link = ps_specifications$inv_link
      )

      mu_list[[label]] <- tryCatch({
        ebmr_b <- EBMRAlgorithmFast3$new("teacher_report", ps_spec_b, dat, W, wt = wt, method = "L-BFGS-B-prog")
        result_b <- ebmr_b$EBMR_IPW(h_nu = h_nu, true_ps = NULL, se.fit = FALSE, wt = wt, method = "L-BFGS-B-prog")
        result_b$mu_ipw
      }, error = function(e) {
        NA
      })
    }

    c(CC = mu_cc_b, unlist(mu_list))
  }

  close(pb_perturb)
  stopCluster(cl_perturb)
  cat("\n")

  # Convert to matrix and extract results
  perturb_matrix <- do.call(rbind, perturb_results_raw)
  perturb_cc <- perturb_matrix[, "CC"]
  perturb_results <- list()
  for (label in model_set_labels) {
    perturb_results[[label]] <- perturb_matrix[, label]
  }

  # Calculate perturbation bootstrap SEs
  perturb_se <- sapply(perturb_results, function(x) sd(x, na.rm = TRUE))
  perturb_se_cc <- sd(perturb_cc, na.rm = TRUE)

  # Save perturbation bootstrap results immediately
  cat("\n  Saving perturbation bootstrap results to MHD_results/...\n")
  saveRDS(perturb_cc, "MHD_results/popmean_perturb_CC_B1000.RDS")
  cat("    Saved: MHD_results/popmean_perturb_CC_B1000.RDS\n")
  for (label in model_set_labels) {
    filename <- sprintf("MHD_results/popmean_perturb_%s_B1000.RDS", label)
    saveRDS(perturb_results[[label]], filename)
    cat(sprintf("    Saved: %s\n", filename))
  }
}

# Print perturbation bootstrap results table
cat(sprintf("\n  %-20s %-18s %10s %12s %12s\n",
            "Label", "Models", "Estimate", "Analytical", "Perturb"))
cat("  ", paste(rep("-", 76), collapse = ""), "\n")

cat(sprintf("  %-20s %-18s %10.4f %12.4f %12.4f\n",
            "CC", "Complete Case", mu_cc, se_cc, perturb_se_cc))

for (j in seq_along(all_model_sets)) {
  label <- model_set_labels[j]
  name <- model_set_names[j]
  res <- results_list[[label]]
  cat(sprintf("  %-20s %-18s %10.4f %12.4f %12.4f\n",
              label, name, res$mu, res$se, perturb_se[label]))
}

cat("  ", paste(rep("-", 76), collapse = ""), "\n")

#------------------------------------------------------------------------------#
# Summary: Create only when both bootstraps are complete
#------------------------------------------------------------------------------#
cat("\n================================================================================\n")
cat("                              FINAL SUMMARY                                     \n")
cat("================================================================================\n\n")

# Print combined comparison table
cat(sprintf("  %-20s %-18s %10s %12s %12s %12s\n",
            "Label", "Models", "Estimate", "Analytical", "Bootstrap", "Perturb"))
cat("  ", paste(rep("-", 88), collapse = ""), "\n")

cat(sprintf("  %-20s %-18s %10.4f %12.4f %12.4f %12.4f\n",
            "CC", "Complete Case", mu_cc, se_cc, boot_se_cc, perturb_se_cc))

cat(sprintf("  %-20s %-18s %10.4f %12s %12.4f %12s\n",
            "MAR", "IPW (MAR)", mu_mar, "-", boot_se_mar, "-"))

for (j in seq_along(all_model_sets)) {
  label <- model_set_labels[j]
  name <- model_set_names[j]
  res <- results_list[[label]]
  cat(sprintf("  %-20s %-18s %10.4f %12.4f %12.4f %12.4f\n",
              label, name, res$mu, res$se, boot_se[label], perturb_se[label]))
}

cat("  ", paste(rep("-", 88), collapse = ""), "\n")

# Create and save summary with all SE estimates
popmean_summary <- data.frame(
  label = c("CC", "MAR", model_set_labels),
  name = c("Complete Case", "IPW (MAR)", model_set_names),
  estimate = c(mu_cc, mu_mar, sapply(results_list, function(x) x$mu)),
  analytical_se = c(se_cc, NA, sapply(results_list, function(x) x$se)),
  bootstrap_se = c(boot_se_cc, boot_se_mar, boot_se[model_set_labels]),
  perturbation_se = c(perturb_se_cc, NA, perturb_se[model_set_labels])
)
saveRDS(popmean_summary, "MHD_results/popmean_summary.RDS")
cat("\n  Saved: MHD_results/popmean_summary.RDS\n")

#------------------------------------------------------------------------------#
# Sensitivity Analysis: Perturb Model 1 with Exponential Tilt
#------------------------------------------------------------------------------#
cat("\n================================================================================\n")
cat("                         SENSITIVITY ANALYSIS                                  \n")
cat("================================================================================\n\n")

# Sensitivity analysis perturbs model 1 using exponential tilt function:
# exp(xi * (teacher_report + health + father + parent_report))
# where xi ranges from 0 to n^(-1/2)

n_xi <- 30
xi_max <- 10*n^(-1/2)
xi_values <- seq(0, xi_max, length.out = n_xi)

cat("Sensitivity analysis parameters:\n")
cat("  xi range: [0, ", round(xi_max, 6), "]\n", sep = "")
cat("  Number of xi values:", n_xi, "\n")
cat("  Perturbing model: 1\n")
cat("  Exponential tilt variables: teacher_report, health, father, parent_report\n\n")

# Use full 3-model result as baseline
baseline_ps_matrix <- result$ps.matrix

# exp_tilt_x_names for the exponential tilt function
exp_tilt_x_names <- c("health", "father")

# Store sensitivity analysis results
sensitivity_results <- data.frame(
  xi = xi_values,
  mu_ipw = NA,
  se_ipw = NA
)

cat("Running sensitivity analysis:\n")
pb_sa <- txtProgressBar(min = 0, max = n_xi, style = 3)

for (i in 1:n_xi) {
  xi <- xi_values[i]

  # Exponential tilt function: exp(xi * (y + health + father))
  exp_tilt <- function(y, x) {
    exp(xi * as.matrix(cbind(y, x)) %*% c(1, 1, 1))
  }

  # Run sensitivity analysis
  sa_result <- tryCatch({
    ebmr$EBMR_IPW_with_locally_misspecified_model(
      ps.matrix = baseline_ps_matrix,
      perturb_ps = 1,
      exp_tilt = exp_tilt,
      exp_tilt_x_names = exp_tilt_x_names,
      h_nu = h_nu,
      nu_init = result$nu.hat,
      se.fit = FALSE,
      method = "L-BFGS-B-prog"
    )
  }, error = function(e) {
    NULL
  })

  if (!is.null(sa_result)) {
    sensitivity_results$mu_ipw[i] <- sa_result$mu_ipw
  }

  setTxtProgressBar(pb_sa, i)
}
close(pb_sa)

# # Save sensitivity analysis results
# for (i in 1:n_xi) {
#   saveRDS(sensitivity_results[i, ], sprintf("MHD_results/popmean_sensitivity_xi%d.RDS", i))
# }

cat("\n\nSensitivity Analysis Results:\n")
cat(sprintf("  %-10s %12s\n", "xi", "mu_ipw"))
cat("  ", paste(rep("-", 24), collapse = ""), "\n")
for (i in seq(1, n_xi, by = 5)) {
  cat(sprintf("  %-10.6f %12.4f\n", sensitivity_results$xi[i], sensitivity_results$mu_ipw[i]))
}
cat("  ... (showing every 5th value)\n\n")

# Report mu_ipw range
cat("Summary:\n")
cat("  mu_ipw range: [", round(min(sensitivity_results$mu_ipw, na.rm = TRUE), 4),
    ", ", round(max(sensitivity_results$mu_ipw, na.rm = TRUE), 4), "]\n", sep = "")
cat("  Baseline (xi=0):", round(sensitivity_results$mu_ipw[1], 4), "\n")
cat("  At xi_max:", round(sensitivity_results$mu_ipw[n_xi], 4), "\n")
cat("  Change from baseline:", round(sensitivity_results$mu_ipw[n_xi] - sensitivity_results$mu_ipw[1], 4), "\n\n")

# Create sensitivity analysis plot
cat("Creating sensitivity analysis plot...\n")
png("MHD_results/popmean_sensitivity_analysis.png", width = 800, height = 600)
plot(sensitivity_results$xi, sensitivity_results$mu_ipw,
     type = "b", pch = 19, col = "blue",
     xlab = expression(xi),
     ylab = expression(hat(mu)[IPW]),
     main = "Sensitivity Analysis: Population Mean Estimation\nPerturbing Model 1 with Exponential Tilt",
     ylim = range(c(sensitivity_results$mu_ipw, mu_cc), na.rm = TRUE))
abline(h = mu_cc, col = "red", lty = 2, lwd = 2)
abline(h = result$mu_ipw, col = "darkgreen", lty = 3, lwd = 2)
legend("topright",
       legend = c("Sensitivity Analysis", "Complete Case", "EBMR (xi=0)"),
       col = c("blue", "red", "darkgreen"),
       lty = c(1, 2, 3), pch = c(19, NA, NA), lwd = 2)
dev.off()
cat("  Saved: MHD_results/popmean_sensitivity_analysis.png\n")

#------------------------------------------------------------------------------#
# Generate LaTeX Table (with outlier-trimmed Bootstrap SE)
#------------------------------------------------------------------------------#
cat("\n================================================================================\n")
cat("                           LATEX TABLE OUTPUT                                  \n")
cat("================================================================================\n\n")

# Function to format CI string
format_ci <- function(est, se) {
  lower <- est - 1.96 * se
  upper <- est + 1.96 * se
  sprintf("{}[%.3f, %.3f]", lower, upper)
}

# Function to compute trimmed SD (IQR x 2 method) and report outliers
compute_trimmed_se <- function(x, label) {
  x <- x[!is.na(x)]
  n_total <- length(x)
  if (n_total < 10) return(list(se = NA, n_outlier = NA, pct_outlier = NA))

  Q1 <- quantile(x, 0.25)
  Q3 <- quantile(x, 0.75)
  IQR_val <- Q3 - Q1
  lower <- Q1 - 2 * IQR_val
  upper <- Q3 + 2 * IQR_val

  outlier_mask <- x < lower | x > upper
  n_outlier <- sum(outlier_mask)
  pct_outlier <- round(100 * n_outlier / n_total, 1)

  x_trim <- x[!outlier_mask]
  se_trim <- sd(x_trim)

  list(se = se_trim, n_outlier = n_outlier, pct_outlier = pct_outlier)
}

# Compute trimmed bootstrap SE for each estimator
cat("Bootstrap SE with outlier removal (IQR x 2 method):\n")
cat(sprintf("  %-10s %12s %12s %15s\n", "Label", "Boot SE", "Outliers", "Pct"))
cat("  ", paste(rep("-", 52), collapse = ""), "\n")

# Complete case
cc_trim <- compute_trimmed_se(boot_cc, "CC")
boot_se_trim_cc <- cc_trim$se
cat(sprintf("  %-10s %12.4f %12d %14.1f%%\n",
            "CC", boot_se_trim_cc, cc_trim$n_outlier, cc_trim$pct_outlier))

# MAR IPW
mar_trim <- compute_trimmed_se(boot_mar, "MAR")
boot_se_trim_mar <- mar_trim$se
cat(sprintf("  %-10s %12.4f %12d %14.1f%%\n",
            "MAR", boot_se_trim_mar, mar_trim$n_outlier, mar_trim$pct_outlier))

# Each model combination
boot_se_trim <- numeric(length(model_set_labels))
names(boot_se_trim) <- model_set_labels
for (label in model_set_labels) {
  trim_result <- compute_trimmed_se(boot_results[[label]], label)
  boot_se_trim[label] <- trim_result$se
  cat(sprintf("  %-10s %12.4f %12d %14.1f%%\n",
              label, trim_result$se, trim_result$n_outlier, trim_result$pct_outlier))
}
cat("  ", paste(rep("-", 52), collapse = ""), "\n\n")

# Build LaTeX table
latex_lines <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Point estimates (PE), estimated standard errors (SE), bootstrap standard errors (Bootstrap SE) based on 1,000 bootstrap replicates, and 95\\% confidence intervals (CI) for the proposed ensemble estimators of $\\mu$.}",
  "\\centering",
  "\\begin{threeparttable}",
  "\\begin{tabularx}{\\textwidth}{l *{4}{>{\\centering\\arraybackslash}X}} ",
  "\\toprule",
  " & PE & SE & Bootstrap SE & 95\\% CI \\\\",
  "\\midrule"
)

# Add rows for each model combination (100, 010, 001, 110, 101, 011, 111)
for (j in seq_along(model_set_labels)) {
  label <- model_set_labels[j]
  res <- results_list[[label]]
  ci_str <- format_ci(res$mu, res$se)
  latex_lines <- c(latex_lines, sprintf(
    "$\\hat{\\mu}_{%s}$ & %.3f & %.3f & %.3f & %s\\\\",
    label, res$mu, res$se, boot_se_trim[label], ci_str
  ))
}

# Add MAR IPW row
ci_mar_str <- format_ci(mu_mar, boot_se_trim_mar)
latex_lines <- c(latex_lines, sprintf(
  "$\\hat{\\mu}_\\text{MAR}$ & %.3f & - & %.3f & %s\\\\",
  mu_mar, boot_se_trim_mar, ci_mar_str
))

# Add complete case row
ci_cc_str <- format_ci(mu_cc, se_cc)
latex_lines <- c(latex_lines, sprintf(
  "$\\hat{\\mu}_\\text{CC}$ & ~%.3f & %.3f & %.3f & %s\\\\",
  mu_cc, se_cc, boot_se_trim_cc, ci_cc_str
))

# Close table
latex_lines <- c(latex_lines,
  "\\bottomrule",
  "\\end{tabularx}",
  "\\end{threeparttable}",
  "\\label{tab:popmean}",
  "\\end{table}"
)

# Print LaTeX table
cat("LaTeX Table:\n\n")
cat(paste(latex_lines, collapse = "\n"))
cat("\n\n")

# Save LaTeX table to file
latex_file <- "MHD_results/popmean_latex_table.tex"
writeLines(latex_lines, latex_file)
cat("  Saved:", latex_file, "\n")

cat("\n================================================================================")
