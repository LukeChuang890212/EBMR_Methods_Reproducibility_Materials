# Logistic Regression Analysis: teacher_report ~ health
# Using EBMR-IPW with All 7 Model Combinations

# Load the EBMRalgorithmFast2 package
library(EBMRalgorithmFast2)

# Source the MHD_functions for gen_data
source("MHD_functions.R")

# Load and generate the mental health data
original_data <- read.csv("data_application.csv")
percent <- original_data$Percentage
n <- 2486
class_n <- round(n * percent / 100)

dat <- gen_data(original_data, class_n, n)
dat$y <- dat$teacher_report

cat("================================================================================\n")
cat("     LOGISTIC REGRESSION ANALYSIS: teacher_report ~ health\n")
cat("     Using EBMR-IPW with All 7 Model Combinations\n")
cat("================================================================================\n\n")

cat("Data Summary:\n")
cat("  Sample size:", nrow(dat), "\n")
cat("  Missing rate:", round(100 * mean(dat$r == 0), 1), "%\n")
cat("  Complete cases:", sum(dat$r == 1), "\n\n")

cat("teacher_report distribution (complete cases):\n")
tr_tab <- table(dat$teacher_report[dat$r == 1])
cat("  Normal (0):", tr_tab[1], sprintf("(%.1f%%)\n", 100*tr_tab[1]/sum(tr_tab)))
cat("  Abnormal (1):", tr_tab[2], sprintf("(%.1f%%)\n", 100*tr_tab[2]/sum(tr_tab)))

cat("\nhealth distribution:\n")
h_tab <- table(dat$health)
cat("  No problems (0):", h_tab[1], sprintf("(%.1f%%)\n", 100*h_tab[1]/sum(h_tab)))
cat("  Health problems (1):", h_tab[2], sprintf("(%.1f%%)\n", 100*h_tab[2]/sum(h_tab)))

# Define full ps_specifications
full_ps_specifications <- list(
  formula.list = list(
    r ~ teacher_report + health + father + teacher_report:health + teacher_report:father,
    r ~ teacher_report + health + parent_report + teacher_report:health + teacher_report:parent_report,
    r ~ teacher_report + parent_report + father + teacher_report:parent_report + teacher_report:father
  ),
  h_x_names.list = list(
    c("health", "father", "parent_report", "fp", "fh", "hp"),
    c("health", "father", "parent_report", "fp", "fh", "hp"),
    c("health", "father", "parent_report", "fp", "fh", "hp")
  ),
  outcome = "teacher_report",
  inv_link = function(eta) 1 / (1 + exp(eta))
)

# Define h_nu function
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

# Compute optimal W matrix
W <- function(g.matrix) {
  solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
}

# Define all 7 model combinations
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

#------------------------------------------------------------------------------#
# Complete Case Analysis (baseline)
#------------------------------------------------------------------------------#
cc_fit <- glm(teacher_report ~ health, data = dat[dat$r == 1, ], family = binomial)
cc_coef <- coef(cc_fit)
cc_se <- summary(cc_fit)$coefficients[, "Std. Error"]
cc_or <- exp(cc_coef[2])

cat("\n================================================================================\n")
cat("                 COMPLETE CASE ANALYSIS (Baseline)\n")
cat("================================================================================\n\n")

cat("Complete Case Logistic Regression:\n")
cat(sprintf("  %-15s %12s %12s %12s %12s\n", "Variable", "Estimate", "Std.Error", "z value", "Pr(>|z|)"))
cat(sprintf("  %s\n", paste(rep("-", 67), collapse="")))
var_names <- c("(Intercept)", "health")
for (i in 1:2) {
  z_val <- cc_coef[i] / cc_se[i]
  p_val <- 2 * pnorm(-abs(z_val))
  stars <- ""
  if (p_val < 0.001) stars <- "***"
  else if (p_val < 0.01) stars <- "**"
  else if (p_val < 0.05) stars <- "*"
  else if (p_val < 0.1) stars <- "."
  cat(sprintf("  %-15s %12.4f %12.4f %12.4f %12.4f %s\n",
              var_names[i], cc_coef[i], cc_se[i], z_val, p_val, stars))
}
cat(sprintf("\nComplete Case OR = %.4f (95%% CI: %.4f - %.4f)\n",
            cc_or, exp(cc_coef[2] - 1.96 * cc_se[2]), exp(cc_coef[2] + 1.96 * cc_se[2])))

#------------------------------------------------------------------------------#
# Fit All 7 Model Combinations
#------------------------------------------------------------------------------#
cat("\n================================================================================\n")
cat("              EBMR-IPW LOGISTIC REGRESSION: ALL MODEL COMBINATIONS\n")
cat("================================================================================\n\n")

results_list <- list()

for (j in seq_along(all_model_sets)) {
  model_set <- all_model_sets[[j]]
  label <- model_set_labels[j]
  name <- model_set_names[j]

  cat(sprintf("Fitting %s (%s)...\n", label, name))

  ps_spec_subset <- list(
    formula.list = full_ps_specifications$formula.list[model_set],
    h_x_names.list = full_ps_specifications$h_x_names.list[model_set],
    outcome = full_ps_specifications$outcome,
    inv_link = full_ps_specifications$inv_link
  )

  ebmr_sub <- EBMRAlgorithmFast2$new("teacher_report", ps_spec_subset, dat, W)

  result_sub <- tryCatch({
    ebmr_sub$EBMR_IPW_regression(
      h_nu = h_nu,
      reg_formula = teacher_report ~ health,
      family = "binomial",
      se.fit = TRUE
    )
  }, error = function(e) {
    cat(sprintf("  Error: %s\n", e$message))
    NULL
  })

  if (!is.null(result_sub)) {
    results_list[[label]] <- list(
      theta = result_sub$theta_hat,
      se = result_sub$se_theta,
      or = exp(result_sub$theta_hat[2]),
      nu = result_sub$nu.hat,
      w = result_sub$w.hat,
      converged = result_sub$convergence
    )
  }
}

#------------------------------------------------------------------------------#
# Summary Table: All Model Combinations
#------------------------------------------------------------------------------#
cat("\n================================================================================\n")
cat("                     SUMMARY: ALL MODEL COMBINATIONS\n")
cat("================================================================================\n\n")

cat("Intercept Estimates:\n")
cat(sprintf("  %-12s %-18s %12s %12s %20s\n",
            "Label", "Models", "Estimate", "SE", "95% CI"))
cat(sprintf("  %s\n", paste(rep("-", 78), collapse="")))

# Complete case
cat(sprintf("  %-12s %-18s %12.4f %12.4f     [%.4f, %.4f]\n",
            "CC", "Complete Case", cc_coef[1], cc_se[1],
            cc_coef[1] - 1.96 * cc_se[1], cc_coef[1] + 1.96 * cc_se[1]))

for (j in seq_along(all_model_sets)) {
  label <- model_set_labels[j]
  name <- model_set_names[j]
  if (!is.null(results_list[[label]])) {
    res <- results_list[[label]]
    ci <- c(res$theta[1] - 1.96 * res$se[1], res$theta[1] + 1.96 * res$se[1])
    cat(sprintf("  %-12s %-18s %12.4f %12.4f     [%.4f, %.4f]\n",
                label, name, res$theta[1], res$se[1], ci[1], ci[2]))
  }
}

cat("\n\nHealth Coefficient Estimates:\n")
cat(sprintf("  %-12s %-18s %12s %12s %20s\n",
            "Label", "Models", "Estimate", "SE", "95% CI"))
cat(sprintf("  %s\n", paste(rep("-", 78), collapse="")))

# Complete case
cat(sprintf("  %-12s %-18s %12.4f %12.4f     [%.4f, %.4f]\n",
            "CC", "Complete Case", cc_coef[2], cc_se[2],
            cc_coef[2] - 1.96 * cc_se[2], cc_coef[2] + 1.96 * cc_se[2]))

for (j in seq_along(all_model_sets)) {
  label <- model_set_labels[j]
  name <- model_set_names[j]
  if (!is.null(results_list[[label]])) {
    res <- results_list[[label]]
    ci <- c(res$theta[2] - 1.96 * res$se[2], res$theta[2] + 1.96 * res$se[2])
    cat(sprintf("  %-12s %-18s %12.4f %12.4f     [%.4f, %.4f]\n",
                label, name, res$theta[2], res$se[2], ci[1], ci[2]))
  }
}

cat("\n\nOdds Ratio for Health Effect:\n")
cat(sprintf("  %-12s %-18s %12s %20s\n",
            "Label", "Models", "OR", "95% CI"))
cat(sprintf("  %s\n", paste(rep("-", 66), collapse="")))

# Complete case
cat(sprintf("  %-12s %-18s %12.4f     [%.4f, %.4f]\n",
            "CC", "Complete Case", cc_or,
            exp(cc_coef[2] - 1.96 * cc_se[2]), exp(cc_coef[2] + 1.96 * cc_se[2])))

for (j in seq_along(all_model_sets)) {
  label <- model_set_labels[j]
  name <- model_set_names[j]
  if (!is.null(results_list[[label]])) {
    res <- results_list[[label]]
    or_ci <- c(exp(res$theta[2] - 1.96 * res$se[2]), exp(res$theta[2] + 1.96 * res$se[2]))
    cat(sprintf("  %-12s %-18s %12.4f     [%.4f, %.4f]\n",
                label, name, res$or, or_ci[1], or_ci[2]))
  }
}

cat("\n\nModel Weights (w.hat) for Each Combination:\n")
cat(sprintf("  %-12s %-18s %40s\n", "Label", "Models", "w.hat"))
cat(sprintf("  %s\n", paste(rep("-", 74), collapse="")))

for (j in seq_along(all_model_sets)) {
  label <- model_set_labels[j]
  name <- model_set_names[j]
  if (!is.null(results_list[[label]])) {
    res <- results_list[[label]]
    w_str <- paste(sprintf("%.4f", res$w), collapse = ", ")
    cat(sprintf("  %-12s %-18s %40s\n", label, name, w_str))
  }
}

#------------------------------------------------------------------------------#
# Bootstrap Standard Errors (Parallel)
#------------------------------------------------------------------------------#
library(parallel)
library(foreach)
library(doSNOW)

B <- 1000
n_cores <- detectCores() - 1

# Output directory (relative path - script runs from Mental_Health_Data_Application)
output_dir <- "MHD_results"

# Check if bootstrap results already exist
boot_file <- file.path(output_dir, "logistic_reg_boot_B1000.RDS")
boot_file_exists <- file.exists(boot_file)

if (boot_file_exists) {
  cat("\n================================================================================\n")
  cat("                 STANDARD BOOTSTRAP (FROM FILE)\n")
  cat("================================================================================\n\n")
  cat("  Bootstrap results already exist. Loading from file...\n")
  boot_results_all <- readRDS(boot_file)
  cat("  Loaded:", boot_file, "\n")
} else {
  cat("\n================================================================================\n")
  cat("                     STANDARD BOOTSTRAP (Parallel)\n")
  cat("================================================================================\n\n")

  cat("Number of bootstrap samples:", B, "\n")
  cat("Number of cores:", n_cores, "\n\n")

  # Set up parallel backend
  cl <- makeCluster(n_cores)
  registerDoSNOW(cl)

  # Progress bar
  pb <- txtProgressBar(max = B, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  cat("Running parallel standard bootstrap:\n")

  boot_results_all <- foreach(
    b = 1:B,
    .options.snow = opts,
    .packages = c("EBMRalgorithmFast2", "Matrix"),
    .export = c("dat", "n", "full_ps_specifications", "W", "h_nu",
                "all_model_sets", "model_set_labels")
  ) %dopar% {
    tryCatch({
      # Bootstrap sample
      set.seed(20200525 + b)
      idx <- sample(1:n, n, replace = TRUE)
      dat_b <- dat[idx, ]

      # Complete case estimate
      cc_fit_b <- glm(teacher_report ~ health, data = dat_b[dat_b$r == 1, ], family = binomial)

      result_vec <- c(cc_intercept = coef(cc_fit_b)[1], cc_health = coef(cc_fit_b)[2])

      # Each model combination
      for (j in seq_along(all_model_sets)) {
        model_set <- all_model_sets[[j]]
        label <- model_set_labels[j]

        ps_spec_b <- list(
          formula.list = full_ps_specifications$formula.list[model_set],
          h_x_names.list = full_ps_specifications$h_x_names.list[model_set],
          outcome = full_ps_specifications$outcome,
          inv_link = full_ps_specifications$inv_link
        )

        res_b <- tryCatch({
          ebmr_b <- EBMRAlgorithmFast2$new("teacher_report", ps_spec_b, dat_b, W)
          result_b <- ebmr_b$EBMR_IPW_regression(
            h_nu = h_nu,
            reg_formula = teacher_report ~ health,
            family = "binomial",
            se.fit = FALSE
          )
          c(result_b$theta_hat[1], result_b$theta_hat[2])
        }, error = function(e) {
          c(NA, NA)
        })

        names(res_b) <- paste0(label, c("_intercept", "_health"))
        result_vec <- c(result_vec, res_b)
      }

      result_vec

    }, error = function(e) {
      # Return NA vector
      result_vec <- rep(NA, 2 + 2 * length(all_model_sets))
      names(result_vec) <- c("cc_intercept", "cc_health",
                              paste0(rep(model_set_labels, each = 2), c("_intercept", "_health")))
      result_vec
    })
  }

  close(pb)
  stopCluster(cl)

  # Save standard bootstrap results
  saveRDS(boot_results_all, boot_file)
  cat("\n\n  Saved:", boot_file, "\n")
}

# Convert to matrix and set column names
boot_matrix <- do.call(rbind, boot_results_all)
expected_colnames <- c("cc_intercept", "cc_health",
                       paste0(rep(model_set_labels, each = 2), c("_intercept", "_health")))
colnames(boot_matrix) <- expected_colnames

cat("\nStandard Bootstrap Results:\n")
cat(sprintf("  Valid bootstrap samples: %d / %d\n", sum(!is.na(boot_matrix[, 1])), B))

#------------------------------------------------------------------------------#
# Perturbation Bootstrap Standard Errors (Parallel)
#------------------------------------------------------------------------------#

# Check if perturbation bootstrap results already exist
perturb_file <- file.path(output_dir, "logistic_reg_perturb_B1000.RDS")
perturb_file_exists <- file.exists(perturb_file)

if (perturb_file_exists) {
  cat("\n================================================================================\n")
  cat("               PERTURBATION BOOTSTRAP (FROM FILE)\n")
  cat("================================================================================\n\n")
  cat("  Perturbation bootstrap results already exist. Loading from file...\n")
  perturb_results_all <- readRDS(perturb_file)
  cat("  Loaded:", perturb_file, "\n")
} else {
  cat("\n================================================================================\n")
  cat("                   PERTURBATION BOOTSTRAP (Parallel)\n")
  cat("================================================================================\n\n")

  cat("Number of perturbation bootstrap samples:", B, "\n")
  cat("Number of cores:", n_cores, "\n\n")

  # Set up parallel backend
  cl_perturb <- makeCluster(n_cores)
  registerDoSNOW(cl_perturb)

  # Progress bar
  pb_perturb <- txtProgressBar(max = B, style = 3)
  progress_perturb <- function(n) setTxtProgressBar(pb_perturb, n)
  opts_perturb <- list(progress = progress_perturb)

  cat("Running parallel perturbation bootstrap:\n")

  perturb_results_all <- foreach(
    b = 1:B,
    .options.snow = opts_perturb,
    .packages = c("EBMRalgorithmFast2", "Matrix"),
    .export = c("dat", "n", "full_ps_specifications", "W", "h_nu",
                "all_model_sets", "model_set_labels")
  ) %dopar% {
    tryCatch({
      # Generate exponential weights
      set.seed(20200525 + b)
      wt <- rexp(n, rate = 1)

      # Weighted complete case estimate
      r_idx <- dat$r == 1
      cc_fit_b <- glm(teacher_report ~ health, data = dat[r_idx, ],
                      family = binomial, weights = wt[r_idx])

      result_vec <- c(cc_intercept = coef(cc_fit_b)[1], cc_health = coef(cc_fit_b)[2])

      # Each model combination
      for (j in seq_along(all_model_sets)) {
        model_set <- all_model_sets[[j]]
        label <- model_set_labels[j]

        ps_spec_b <- list(
          formula.list = full_ps_specifications$formula.list[model_set],
          h_x_names.list = full_ps_specifications$h_x_names.list[model_set],
          outcome = full_ps_specifications$outcome,
          inv_link = full_ps_specifications$inv_link
        )

        res_b <- tryCatch({
          ebmr_b <- EBMRAlgorithmFast2$new("teacher_report", ps_spec_b, dat, W, wt = wt)
          result_b <- ebmr_b$EBMR_IPW_regression(
            h_nu = h_nu,
            reg_formula = teacher_report ~ health,
            family = "binomial",
            se.fit = FALSE,
            wt = wt
          )
          c(result_b$theta_hat[1], result_b$theta_hat[2])
        }, error = function(e) {
          c(NA, NA)
        })

        names(res_b) <- paste0(label, c("_intercept", "_health"))
        result_vec <- c(result_vec, res_b)
      }

      result_vec

    }, error = function(e) {
      # Return NA vector
      result_vec <- rep(NA, 2 + 2 * length(all_model_sets))
      names(result_vec) <- c("cc_intercept", "cc_health",
                              paste0(rep(model_set_labels, each = 2), c("_intercept", "_health")))
      result_vec
    })
  }

  close(pb_perturb)
  stopCluster(cl_perturb)

  # Save perturbation bootstrap results
  saveRDS(perturb_results_all, perturb_file)
  cat("\n\n  Saved:", perturb_file, "\n")
}

# Convert to matrix and set column names
perturb_matrix <- do.call(rbind, perturb_results_all)
colnames(perturb_matrix) <- expected_colnames

cat("\nPerturbation Bootstrap Results:\n")
cat(sprintf("  Valid bootstrap samples: %d / %d\n", sum(!is.na(perturb_matrix[, 1])), B))

#------------------------------------------------------------------------------#
# Bootstrap to Evaluate Variation of Analytical SE Estimates
#------------------------------------------------------------------------------#

# Check if analytical SE bootstrap results already exist
ase_boot_file <- file.path(output_dir, "logistic_reg_ase_boot_B1000.RDS")
ase_boot_file_exists <- file.exists(ase_boot_file)

if (ase_boot_file_exists) {
  cat("\n================================================================================\n")
  cat("           ANALYTICAL SE VARIATION BOOTSTRAP (FROM FILE)\n")
  cat("================================================================================\n\n")
  cat("  Analytical SE bootstrap results already exist. Loading from file...\n")
  ase_boot_results_all <- readRDS(ase_boot_file)
  cat("  Loaded:", ase_boot_file, "\n")
} else {
  cat("\n================================================================================\n")
  cat("           ANALYTICAL SE VARIATION BOOTSTRAP (Parallel)\n")
  cat("================================================================================\n\n")

  cat("This bootstrap evaluates how the analytical SE estimate varies across samples.\n")
  cat("Each bootstrap sample computes both theta_hat AND its analytical SE.\n\n")

  cat("Number of bootstrap samples:", B, "\n")
  cat("Number of cores:", n_cores, "\n\n")

  # Set up parallel backend
  cl_ase <- makeCluster(n_cores)
  registerDoSNOW(cl_ase)

  # Progress bar
  pb_ase <- txtProgressBar(max = B, style = 3)
  progress_ase <- function(n) setTxtProgressBar(pb_ase, n)
  opts_ase <- list(progress = progress_ase)

  cat("Running parallel analytical SE variation bootstrap:\n")

  ase_boot_results_all <- foreach(
    b = 1:B,
    .options.snow = opts_ase,
    .packages = c("EBMRalgorithmFast2", "Matrix"),
    .export = c("dat", "n", "full_ps_specifications", "W", "h_nu",
                "all_model_sets", "model_set_labels")
  ) %dopar% {
    tryCatch({
      # Bootstrap sample
      set.seed(20200525 + b)
      idx <- sample(1:n, n, replace = TRUE)
      dat_b <- dat[idx, ]

      # Complete case estimate with SE
      cc_fit_b <- glm(teacher_report ~ health, data = dat_b[dat_b$r == 1, ], family = binomial)
      cc_coef_b <- coef(cc_fit_b)
      cc_se_b <- summary(cc_fit_b)$coefficients[, "Std. Error"]

      result_vec <- c(
        cc_intercept = cc_coef_b[1], cc_health = cc_coef_b[2],
        cc_se_intercept = cc_se_b[1], cc_se_health = cc_se_b[2]
      )

      # Each model combination
      for (j in seq_along(all_model_sets)) {
        model_set <- all_model_sets[[j]]
        label <- model_set_labels[j]

        ps_spec_b <- list(
          formula.list = full_ps_specifications$formula.list[model_set],
          h_x_names.list = full_ps_specifications$h_x_names.list[model_set],
          outcome = full_ps_specifications$outcome,
          inv_link = full_ps_specifications$inv_link
        )

        res_b <- tryCatch({
          ebmr_b <- EBMRAlgorithmFast2$new("teacher_report", ps_spec_b, dat_b, W)
          result_b <- ebmr_b$EBMR_IPW_regression(
            h_nu = h_nu,
            reg_formula = teacher_report ~ health,
            family = "binomial",
            se.fit = TRUE  # Compute analytical SE in each bootstrap sample
          )
          c(result_b$theta_hat[1], result_b$theta_hat[2],
            result_b$se_theta[1], result_b$se_theta[2])
        }, error = function(e) {
          c(NA, NA, NA, NA)
        })

        names(res_b) <- paste0(label, c("_intercept", "_health", "_se_intercept", "_se_health"))
        result_vec <- c(result_vec, res_b)
      }

      result_vec

    }, error = function(e) {
      # Return NA vector
      result_vec <- rep(NA, 4 + 4 * length(all_model_sets))
      names(result_vec) <- c("cc_intercept", "cc_health", "cc_se_intercept", "cc_se_health",
                              paste0(rep(model_set_labels, each = 4),
                                     c("_intercept", "_health", "_se_intercept", "_se_health")))
      result_vec
    })
  }

  close(pb_ase)
  stopCluster(cl_ase)

  # Save analytical SE bootstrap results
  saveRDS(ase_boot_results_all, ase_boot_file)
  cat("\n\n  Saved:", ase_boot_file, "\n")
}

# Convert to matrix
ase_boot_matrix <- do.call(rbind, ase_boot_results_all)
ase_expected_colnames <- c("cc_intercept", "cc_health", "cc_se_intercept", "cc_se_health",
                            paste0(rep(model_set_labels, each = 4),
                                   c("_intercept", "_health", "_se_intercept", "_se_health")))
colnames(ase_boot_matrix) <- ase_expected_colnames

cat("\nAnalytical SE Variation Bootstrap Results:\n")
cat(sprintf("  Valid bootstrap samples: %d / %d\n", sum(!is.na(ase_boot_matrix[, 1])), B))

#------------------------------------------------------------------------------#
# Check for Outliers in ASE Bootstrap Results
#------------------------------------------------------------------------------#

cat("\n================================================================================\n")
cat("                OUTLIER DIAGNOSTICS FOR BOOTSTRAP ASE\n")
cat("================================================================================\n\n")

# Function to detect outliers using IQR method
detect_outliers_iqr <- function(x, multiplier = 1.5) {
  x <- x[!is.na(x)]
  q1 <- quantile(x, 0.25)
  q3 <- quantile(x, 0.75)
  iqr <- q3 - q1
  lower <- q1 - multiplier * iqr
  upper <- q3 + multiplier * iqr
  outliers <- x < lower | x > upper
  list(
    n_outliers = sum(outliers),
    pct_outliers = 100 * mean(outliers),
    lower = lower,
    upper = upper,
    min = min(x),
    max = max(x),
    median = median(x)
  )
}

cat("Health Coefficient ASE - Outlier Check (IQR x 1.5):\n\n")
cat(sprintf("  %-12s %-18s %8s %8s %10s %10s %10s %10s\n",
            "Label", "Models", "N_out", "%_out", "Min", "Median", "Max", "IQR_upper"))
cat(sprintf("  %s\n", paste(rep("-", 96), collapse="")))

# Complete case
ase_cc_diag <- detect_outliers_iqr(ase_boot_matrix[, "cc_se_health"])
cat(sprintf("  %-12s %-18s %8d %8.1f %10.4f %10.4f %10.4f %10.4f\n",
            "CC", "Complete Case", ase_cc_diag$n_outliers, ase_cc_diag$pct_outliers,
            ase_cc_diag$min, ase_cc_diag$median, ase_cc_diag$max, ase_cc_diag$upper))

for (j in seq_along(all_model_sets)) {
  label <- model_set_labels[j]
  name <- model_set_names[j]
  col_name <- paste0(label, "_se_health")
  ase_diag <- detect_outliers_iqr(ase_boot_matrix[, col_name])
  cat(sprintf("  %-12s %-18s %8d %8.1f %10.4f %10.4f %10.4f %10.4f\n",
              label, name, ase_diag$n_outliers, ase_diag$pct_outliers,
              ase_diag$min, ase_diag$median, ase_diag$max, ase_diag$upper))
}

# Function to compute trimmed statistics (remove outliers)
trimmed_stats <- function(x, multiplier = 1.5) {
  x <- x[!is.na(x)]
  q1 <- quantile(x, 0.25)
  q3 <- quantile(x, 0.75)
  iqr <- q3 - q1
  lower <- q1 - multiplier * iqr
  upper <- q3 + multiplier * iqr
  x_trimmed <- x[x >= lower & x <= upper]
  list(
    mean_raw = mean(x),
    sd_raw = sd(x),
    mean_trimmed = mean(x_trimmed),
    sd_trimmed = sd(x_trimmed),
    n_removed = length(x) - length(x_trimmed)
  )
}

#------------------------------------------------------------------------------#
# Summary: Variation of Analytical SE Estimates
#------------------------------------------------------------------------------#

cat("\n================================================================================\n")
cat("                VARIATION OF ANALYTICAL SE ESTIMATES\n")
cat("================================================================================\n\n")

cat("This shows how the analytical SE estimate varies across bootstrap samples.\n")
cat("SD(ASE) = standard deviation of analytical SE across bootstrap samples.\n")
cat("SD_trim = SD after removing outliers (IQR x 1.5 method).\n")
cat("CV(ASE) = coefficient of variation (SD/Mean) of analytical SE.\n\n")

cat("Health Coefficient (log OR) - Analytical SE Variation:\n\n")
cat(sprintf("  %-12s %-18s %10s %10s %10s %10s %10s\n",
            "Label", "Models", "ASE", "Mean(ASE)", "SD(ASE)", "SD_trim", "CV(ASE)"))
cat(sprintf("  %s\n", paste(rep("-", 88), collapse="")))

# Complete case
ase_original_cc <- cc_se[2]
ase_boot_cc <- ase_boot_matrix[, "cc_se_health"]
ase_stats_cc <- trimmed_stats(ase_boot_cc)
ase_mean_cc <- ase_stats_cc$mean_raw
ase_sd_cc <- ase_stats_cc$sd_raw
ase_sd_trim_cc <- ase_stats_cc$sd_trimmed
ase_cv_cc <- ase_sd_cc / ase_mean_cc
cat(sprintf("  %-12s %-18s %10.4f %10.4f %10.4f %10.4f %10.4f\n",
            "CC", "Complete Case", ase_original_cc, ase_mean_cc, ase_sd_cc, ase_sd_trim_cc, ase_cv_cc))

for (j in seq_along(all_model_sets)) {
  label <- model_set_labels[j]
  name <- model_set_names[j]
  if (!is.null(results_list[[label]])) {
    res <- results_list[[label]]
    col_name <- paste0(label, "_se_health")
    ase_boot <- ase_boot_matrix[, col_name]
    ase_stats <- trimmed_stats(ase_boot)
    ase_mean <- ase_stats$mean_raw
    ase_sd <- ase_stats$sd_raw
    ase_sd_trim <- ase_stats$sd_trimmed
    ase_cv <- ase_sd / ase_mean
    cat(sprintf("  %-12s %-18s %10.4f %10.4f %10.4f %10.4f %10.4f\n",
                label, name, res$se[2], ase_mean, ase_sd, ase_sd_trim, ase_cv))
  }
}
cat(sprintf("  %s\n", paste(rep("-", 88), collapse="")))

cat("\n\nIntercept - Analytical SE Variation:\n\n")
cat(sprintf("  %-12s %-18s %10s %10s %10s %10s %10s\n",
            "Label", "Models", "ASE", "Mean(ASE)", "SD(ASE)", "SD_trim", "CV(ASE)"))
cat(sprintf("  %s\n", paste(rep("-", 88), collapse="")))

# Complete case
ase_original_cc_int <- cc_se[1]
ase_boot_cc_int <- ase_boot_matrix[, "cc_se_intercept"]
ase_stats_cc_int <- trimmed_stats(ase_boot_cc_int)
ase_mean_cc_int <- ase_stats_cc_int$mean_raw
ase_sd_cc_int <- ase_stats_cc_int$sd_raw
ase_sd_trim_cc_int <- ase_stats_cc_int$sd_trimmed
ase_cv_cc_int <- ase_sd_cc_int / ase_mean_cc_int
cat(sprintf("  %-12s %-18s %10.4f %10.4f %10.4f %10.4f %10.4f\n",
            "CC", "Complete Case", ase_original_cc_int, ase_mean_cc_int, ase_sd_cc_int, ase_sd_trim_cc_int, ase_cv_cc_int))

for (j in seq_along(all_model_sets)) {
  label <- model_set_labels[j]
  name <- model_set_names[j]
  if (!is.null(results_list[[label]])) {
    res <- results_list[[label]]
    col_name <- paste0(label, "_se_intercept")
    ase_boot <- ase_boot_matrix[, col_name]
    ase_stats <- trimmed_stats(ase_boot)
    ase_mean <- ase_stats$mean_raw
    ase_sd <- ase_stats$sd_raw
    ase_sd_trim <- ase_stats$sd_trimmed
    ase_cv <- ase_sd / ase_mean
    cat(sprintf("  %-12s %-18s %10.4f %10.4f %10.4f %10.4f %10.4f\n",
                label, name, res$se[1], ase_mean, ase_sd, ase_sd_trim, ase_cv))
  }
}
cat(sprintf("  %s\n", paste(rep("-", 88), collapse="")))

#------------------------------------------------------------------------------#
# Final Summary with All SE Estimates
#------------------------------------------------------------------------------#

cat("\n================================================================================\n")
cat("                    FINAL SUMMARY: SE COMPARISON\n")
cat("================================================================================\n\n")

cat("Health Coefficient (log OR) - All Methods:\n\n")
cat(sprintf("  %-12s %-18s %10s %10s %10s %10s %10s\n",
            "Label", "Models", "Estimate", "ASE", "Boot.SE", "Perturb", "SD_trim"))
cat(sprintf("  %s\n", paste(rep("-", 92), collapse="")))

# Complete case
boot_se_cc <- sd(boot_matrix[, "cc_health"], na.rm = TRUE)
perturb_se_cc <- sd(perturb_matrix[, "cc_health"], na.rm = TRUE)
cat(sprintf("  %-12s %-18s %10.4f %10.4f %10.4f %10.4f %10.4f\n",
            "CC", "Complete Case", cc_coef[2], cc_se[2], boot_se_cc, perturb_se_cc, ase_sd_trim_cc))

for (j in seq_along(all_model_sets)) {
  label <- model_set_labels[j]
  name <- model_set_names[j]
  if (!is.null(results_list[[label]])) {
    res <- results_list[[label]]
    col_name <- paste0(label, "_health")
    boot_se <- sd(boot_matrix[, col_name], na.rm = TRUE)
    perturb_se <- sd(perturb_matrix[, col_name], na.rm = TRUE)
    ase_stats_j <- trimmed_stats(ase_boot_matrix[, paste0(label, "_se_health")])
    ase_sd_trim <- ase_stats_j$sd_trimmed
    cat(sprintf("  %-12s %-18s %10.4f %10.4f %10.4f %10.4f %10.4f\n",
                label, name, res$theta[2], res$se[2], boot_se, perturb_se, ase_sd_trim))
  }
}
cat(sprintf("  %s\n", paste(rep("-", 92), collapse="")))

cat("\nNote: ASE = Analytical SE, Boot.SE = SD of bootstrap estimates,\n")
cat("      Perturb = SD of perturbation bootstrap estimates, SD_trim = trimmed SD of analytical SE (outliers removed).\n")

cat("\n\nOdds Ratio for Health Effect - Bootstrap 95% CI:\n\n")
cat(sprintf("  %-12s %-18s %10s %25s %25s\n",
            "Label", "Models", "OR", "Std.Boot 95% CI", "Perturb 95% CI"))
cat(sprintf("  %s\n", paste(rep("-", 96), collapse="")))

# Complete case
boot_or_cc <- exp(boot_matrix[, "cc_health"])
perturb_or_cc <- exp(perturb_matrix[, "cc_health"])
boot_or_ci_cc <- quantile(boot_or_cc, c(0.025, 0.975), na.rm = TRUE)
perturb_or_ci_cc <- quantile(perturb_or_cc, c(0.025, 0.975), na.rm = TRUE)
cat(sprintf("  %-12s %-18s %10.4f     [%.4f, %.4f]     [%.4f, %.4f]\n",
            "CC", "Complete Case", cc_or,
            boot_or_ci_cc[1], boot_or_ci_cc[2],
            perturb_or_ci_cc[1], perturb_or_ci_cc[2]))

for (j in seq_along(all_model_sets)) {
  label <- model_set_labels[j]
  name <- model_set_names[j]
  if (!is.null(results_list[[label]])) {
    res <- results_list[[label]]
    col_name <- paste0(label, "_health")
    boot_or <- exp(boot_matrix[, col_name])
    perturb_or <- exp(perturb_matrix[, col_name])
    boot_or_ci <- quantile(boot_or, c(0.025, 0.975), na.rm = TRUE)
    perturb_or_ci <- quantile(perturb_or, c(0.025, 0.975), na.rm = TRUE)
    cat(sprintf("  %-12s %-18s %10.4f     [%.4f, %.4f]     [%.4f, %.4f]\n",
                label, name, res$or,
                boot_or_ci[1], boot_or_ci[2],
                perturb_or_ci[1], perturb_or_ci[2]))
  }
}
cat(sprintf("  %s\n", paste(rep("-", 96), collapse="")))

#------------------------------------------------------------------------------#
# Save Final Summary
#------------------------------------------------------------------------------#

# Build summary object
logistic_reg_summary <- list(
  # Complete case
  cc = list(
    theta = cc_coef,
    se = cc_se,
    or = cc_or,
    boot_se = c(intercept = sd(boot_matrix[, "cc_intercept"], na.rm = TRUE),
                health = boot_se_cc),
    perturb_se = c(intercept = sd(perturb_matrix[, "cc_intercept"], na.rm = TRUE),
                   health = perturb_se_cc),
    # Analytical SE variation
    ase_mean = c(intercept = ase_mean_cc_int, health = ase_mean_cc),
    ase_sd = c(intercept = ase_sd_cc_int, health = ase_sd_cc),
    ase_sd_trim = c(intercept = ase_sd_trim_cc_int, health = ase_sd_trim_cc),
    ase_cv = c(intercept = ase_cv_cc_int, health = ase_cv_cc)
  ),
  # All model combinations
  results = results_list,
  # Bootstrap SEs for each combination
  boot_se = sapply(model_set_labels, function(label) {
    c(intercept = sd(boot_matrix[, paste0(label, "_intercept")], na.rm = TRUE),
      health = sd(boot_matrix[, paste0(label, "_health")], na.rm = TRUE))
  }),
  perturb_se = sapply(model_set_labels, function(label) {
    c(intercept = sd(perturb_matrix[, paste0(label, "_intercept")], na.rm = TRUE),
      health = sd(perturb_matrix[, paste0(label, "_health")], na.rm = TRUE))
  }),
  # Analytical SE variation (mean, SD, CV across bootstrap samples)
  ase_mean = sapply(model_set_labels, function(label) {
    c(intercept = mean(ase_boot_matrix[, paste0(label, "_se_intercept")], na.rm = TRUE),
      health = mean(ase_boot_matrix[, paste0(label, "_se_health")], na.rm = TRUE))
  }),
  ase_sd = sapply(model_set_labels, function(label) {
    c(intercept = sd(ase_boot_matrix[, paste0(label, "_se_intercept")], na.rm = TRUE),
      health = sd(ase_boot_matrix[, paste0(label, "_se_health")], na.rm = TRUE))
  }),
  ase_sd_trim = sapply(model_set_labels, function(label) {
    stats_int <- trimmed_stats(ase_boot_matrix[, paste0(label, "_se_intercept")])
    stats_h <- trimmed_stats(ase_boot_matrix[, paste0(label, "_se_health")])
    c(intercept = stats_int$sd_trimmed, health = stats_h$sd_trimmed)
  }),
  ase_cv = sapply(model_set_labels, function(label) {
    ase_int <- ase_boot_matrix[, paste0(label, "_se_intercept")]
    ase_h <- ase_boot_matrix[, paste0(label, "_se_health")]
    c(intercept = sd(ase_int, na.rm = TRUE) / mean(ase_int, na.rm = TRUE),
      health = sd(ase_h, na.rm = TRUE) / mean(ase_h, na.rm = TRUE))
  }),
  # Model labels
  model_set_labels = model_set_labels,
  model_set_names = model_set_names
)

summary_file <- file.path(output_dir, "logistic_reg_summary.RDS")
saveRDS(logistic_reg_summary, summary_file)
cat("\n  Saved:", summary_file, "\n")

cat("\n================================================================================\n")
cat("                              ANALYSIS COMPLETE\n")
cat("================================================================================\n")
