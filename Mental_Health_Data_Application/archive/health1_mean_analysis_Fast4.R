#------------------------------------------------------------------------------#
# Health=1 Subgroup Mean Estimation: E[teacher_report | health=1]
# Using EBMRalgorithmFast4
# Models: r~teacher_report+father, r~teacher_report+parent_report, r~father+parent_report
# h_alpha: (father, parent_report), h_nu: (father, parent_report, fp)
#------------------------------------------------------------------------------#

setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Mental_Health_Data_Application")
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)
library(tidyverse)
library(numDeriv)
library(Matrix)

source("MHD_functions.R")

#------------------------------------------------------------------------------#
# Read and prepare data — subset to health=1
#------------------------------------------------------------------------------#
original_data <- read.csv("data_application.csv")
percent <- original_data$Percentage
n_full <- 2486
class_n <- round(n_full * percent / 100)

dat_full <- gen_data(original_data, class_n, n_full)
dat_full$y <- dat_full$teacher_report

# Subset to health=1
dat <- dat_full[dat_full$health == 1, ]
n <- nrow(dat)
missing_rate <- 1 - mean(dat$r)
mu_cc <- mean(dat$teacher_report[dat$r == 1])
se_cc <- sd(dat$teacher_report[dat$r == 1]) / sqrt(sum(dat$r))

cat("================================================================================\n")
cat("     HEALTH=1 SUBGROUP MEAN ESTIMATION: E[teacher_report | health=1]           \n")
cat("     Using EBMRalgorithmFast4                                                  \n")
cat("================================================================================\n\n")

cat("Sample size (health=1):", n, "\n")
cat("Missing rate:", round(missing_rate * 100, 1), "%\n")
cat("Complete case mean:", round(mu_cc, 4), "\n\n")

#------------------------------------------------------------------------------#
# MAR Analysis: IPW with saturated propensity score model
#------------------------------------------------------------------------------#
cat("================================================================================\n")
cat("                     MAR ANALYSIS: IPW ESTIMATOR                               \n")
cat("================================================================================\n\n")

ps_mar_fit <- glm(r ~ father * parent_report, data = dat, family = binomial)
ps_mar <- fitted(ps_mar_fit)
mu_mar <- sum(dat$r * dat$teacher_report / ps_mar) / n

cat(sprintf("MAR IPW estimate: %.4f\n\n", mu_mar))

#------------------------------------------------------------------------------#
# EBMR Setup
#------------------------------------------------------------------------------#
W <- function(g.matrix) {
  solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
}

ps_specifications <- list(
  formula.list = list(
    r ~ teacher_report + father,
    r ~ teacher_report + parent_report,
    r ~ father + parent_report
  ),
  h_alpha.list = list(
    c("father", "parent_report"),
    c("father", "parent_report"),
    c("father", "parent_report")
  ),
  outcome = "teacher_report",
  inv_link = function(eta) 1 / (1 + exp(eta))
)

h_nu <- function(data) {
  cbind(
    father = data$father,
    parent_report = data$parent_report,
    fp = data$father * data$parent_report
  )
}

#------------------------------------------------------------------------------#
# Fit EBMR
#------------------------------------------------------------------------------#
ebmr <- EBMRAlgorithmFast4$new("teacher_report", ps_specifications, dat, W)
result <- ebmr$EBMR_IPW(h_nu = h_nu, true_ps = NULL)

#------------------------------------------------------------------------------#
# Model Summaries
#------------------------------------------------------------------------------#
cat("================================================================================\n")
cat("                       PROPENSITY SCORE MODEL SUMMARIES                        \n")
cat("================================================================================\n\n")

for (i in seq_along(ebmr$ps_fit.list)) {
  ps_fit <- ebmr$ps_fit.list[[i]]
  formula_str <- deparse(ps_specifications$formula.list[[i]])
  cat(sprintf("Model %d: %s\n", i, formula_str))

  coef <- ps_fit$coefficients
  se <- ps_fit$se
  z_val <- coef / se
  p_val <- 2 * pnorm(-abs(z_val))
  var_names <- colnames(ps_fit$design_matrix)
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
  fitted_ps <- ps_fit$fitted.values
  cat(sprintf("  Fitted PS range: [%.4f, %.4f], >0.99: %d, <0.01: %d\n\n",
              min(fitted_ps), max(fitted_ps), sum(fitted_ps > 0.99), sum(fitted_ps < 0.01)))
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
# All Model Combinations
#------------------------------------------------------------------------------#
cat("\n================================================================================\n")
cat("                     ALL MODEL COMBINATIONS                                    \n")
cat("================================================================================\n\n")

all_model_sets <- list(c(1), c(2), c(3), c(1,2), c(1,3), c(2,3), c(1,2,3))
model_set_labels <- c("100", "010", "001", "110", "101", "011", "111")
model_set_names <- c(
  "pi1 only", "pi2 only", "pi3 only",
  "pi1 + pi2", "pi1 + pi3", "pi2 + pi3", "pi1 + pi2 + pi3"
)

results_list <- list()

for (j in seq_along(all_model_sets)) {
  model_set <- all_model_sets[[j]]
  label <- model_set_labels[j]
  result_sub <- ebmr$EBMR_IPW(h_nu = h_nu, model_indices = model_set, true_ps = NULL)
  results_list[[label]] <- list(
    mu = result_sub$mu_ipw,
    se = result_sub$se_ipw,
    nu = result_sub$nu.hat,
    w = result_sub$w.hat
  )
}

cat(sprintf("  %-20s %-18s %10s %10s %20s\n",
            "Label", "Models", "Estimate", "SE", "95% CI"))
cat("  ", paste(rep("-", 82), collapse = ""), "\n")

cat(sprintf("  %-20s %-18s %10.4f %10.4f     [%.4f, %.4f]\n",
            "CC", "Complete Case", mu_cc, se_cc,
            mu_cc - 1.96 * se_cc, mu_cc + 1.96 * se_cc))

cat(sprintf("  %-20s %-18s %10.4f %10s     %s\n",
            "MAR", "IPW (MAR)", mu_mar, "-", "-"))

for (j in seq_along(all_model_sets)) {
  label <- model_set_labels[j]
  name <- model_set_names[j]
  res <- results_list[[label]]
  ci <- c(res$mu - 1.96 * res$se, res$mu + 1.96 * res$se)
  cat(sprintf("  %-20s %-18s %10.4f %10.4f     [%.4f, %.4f]\n",
              label, name, res$mu, res$se, ci[1], ci[2]))
}

cat("  ", paste(rep("-", 82), collapse = ""), "\n")

# Weights
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
# Bootstrap Standard Errors (Parallel)
#------------------------------------------------------------------------------#
library(parallel)
library(foreach)
library(doSNOW)

B <- 1000
n_cores <- detectCores() - 1

boot_file <- "MHD_results/health1_boot_Fast4_111_B1000.RDS"

if (file.exists(boot_file)) {
  cat("\n================================================================================\n")
  cat("                     BOOTSTRAP STANDARD ERRORS (FROM FILE)                     \n")
  cat("================================================================================\n\n")

  boot_cc <- readRDS("MHD_results/health1_boot_Fast4_CC_B1000.RDS")
  boot_se_cc <- sd(boot_cc, na.rm = TRUE)
  boot_se <- numeric(length(model_set_labels))
  names(boot_se) <- model_set_labels
  boot_results <- list()
  for (label in model_set_labels) {
    bf <- sprintf("MHD_results/health1_boot_Fast4_%s_B1000.RDS", label)
    if (file.exists(bf)) {
      boot_results[[label]] <- readRDS(bf)
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

  cl <- makeCluster(n_cores)
  registerDoSNOW(cl)

  pb <- txtProgressBar(max = B, style = 3)
  progress <- function(nn) setTxtProgressBar(pb, nn)
  opts <- list(progress = progress)

  clusterExport(cl, c("dat_full"), envir = environment())
  clusterEvalQ(cl, devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE))

  boot_results_raw <- foreach(
    b = 1:B,
    .combine = 'rbind',
    .options.snow = opts,
    .packages = c("stringr", "Matrix", "numDeriv")
  ) %dopar% {
    set.seed(12345 + b)
    # Resample from the full dataset, then subset to health=1
    idx <- sample(1:nrow(dat_full), nrow(dat_full), replace = TRUE)
    dat_b_full <- dat_full[idx, ]
    dat_b <- dat_b_full[dat_b_full$health == 1, ]

    W <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
    h_nu <- function(data) cbind(father = data$father, parent_report = data$parent_report,
                                  fp = data$father * data$parent_report)

    ps_spec <- list(
      formula.list = list(
        r ~ teacher_report + father,
        r ~ teacher_report + parent_report,
        r ~ father + parent_report
      ),
      h_alpha.list = list(
        c("father", "parent_report"),
        c("father", "parent_report"),
        c("father", "parent_report")
      ),
      outcome = "teacher_report",
      inv_link = function(eta) 1 / (1 + exp(eta))
    )

    mu_cc_b <- mean(dat_b$teacher_report[dat_b$r == 1])

    all_sets <- list(c(1), c(2), c(3), c(1,2), c(1,3), c(2,3), c(1,2,3))
    labels <- c("100", "010", "001", "110", "101", "011", "111")
    mu_list <- setNames(rep(NA, 7), labels)

    tryCatch({
      ebmr_b <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat_b, W)
      for (j in 1:7) {
        tryCatch({
          res_b <- ebmr_b$EBMR_IPW(h_nu = h_nu, model_indices = all_sets[[j]],
                                     true_ps = NULL, se.fit = FALSE)
          mu_list[j] <- res_b$mu_ipw
        }, error = function(e) {})
      }
    }, error = function(e) {})

    c(CC = mu_cc_b, mu_list)
  }

  close(pb)
  stopCluster(cl)
  cat("\n")

  boot_cc <- boot_results_raw[, "CC"]
  boot_results <- list()
  for (label in model_set_labels) {
    boot_results[[label]] <- boot_results_raw[, label]
  }

  boot_se <- sapply(boot_results, function(x) sd(x, na.rm = TRUE))
  boot_se_cc <- sd(boot_cc, na.rm = TRUE)

  # Save
  saveRDS(boot_cc, "MHD_results/health1_boot_Fast4_CC_B1000.RDS")
  for (label in model_set_labels) {
    saveRDS(boot_results[[label]],
            sprintf("MHD_results/health1_boot_Fast4_%s_B1000.RDS", label))
  }
  cat("  Bootstrap results saved.\n")
}

#------------------------------------------------------------------------------#
# Outlier-trimmed Bootstrap SE (using clean_sim_result logic: IQR x 3, 1% cap)
#------------------------------------------------------------------------------#
compute_trimmed_boot_se <- function(x, label) {
  x <- x[!is.na(x)]
  n_total <- length(x)
  if (n_total < 10) return(list(se_raw = NA, se_trim = NA, n_na = 0, n_outlier = NA))

  Q1 <- quantile(x, 0.25); Q3 <- quantile(x, 0.75)
  IQR_val <- Q3 - Q1
  is_outlier <- x < (Q1 - 3 * IQR_val) | x > (Q3 + 3 * IQR_val)
  max_remove <- floor(0.01 * n_total)
  if (sum(is_outlier) > max_remove && max_remove > 0) {
    dist_med <- abs(x - median(x))
    oi <- which(is_outlier)
    keep <- oi[order(dist_med[oi], decreasing = TRUE)[(max_remove + 1):length(oi)]]
    is_outlier[keep] <- FALSE
  }
  n_outlier <- sum(is_outlier)
  list(se_raw = sd(x), se_trim = sd(x[!is_outlier]), n_na = 0, n_outlier = n_outlier)
}

#------------------------------------------------------------------------------#
# Final Summary
#------------------------------------------------------------------------------#
cat("\n================================================================================\n")
cat("                              FINAL SUMMARY                                     \n")
cat("================================================================================\n\n")

cat("Bootstrap SE with outlier removal (IQR x 3, 1%% cap):\n")
cat(sprintf("  %-10s %12s %12s %12s\n", "Label", "Boot(raw)", "Boot(trim)", "Outliers"))
cat("  ", paste(rep("-", 50), collapse = ""), "\n")

cc_trim <- compute_trimmed_boot_se(boot_cc, "CC")
boot_se_trim_cc <- cc_trim$se_trim
cat(sprintf("  %-10s %12.4f %12.4f %12d\n",
            "CC", cc_trim$se_raw, cc_trim$se_trim, cc_trim$n_outlier))

boot_se_trim <- numeric(length(model_set_labels))
names(boot_se_trim) <- model_set_labels
for (label in model_set_labels) {
  trim_result <- compute_trimmed_boot_se(boot_results[[label]], label)
  boot_se_trim[label] <- trim_result$se_trim
  cat(sprintf("  %-10s %12.4f %12.4f %12d\n",
              label, trim_result$se_raw, trim_result$se_trim, trim_result$n_outlier))
}
cat("  ", paste(rep("-", 50), collapse = ""), "\n")

cat(sprintf("\n  %-20s %-18s %10s %12s %12s %12s\n",
            "Label", "Models", "Estimate", "Analytical", "Boot (raw)", "Boot (trim)"))
cat("  ", paste(rep("-", 88), collapse = ""), "\n")

cat(sprintf("  %-20s %-18s %10.4f %12.4f %12.4f %12.4f\n",
            "CC", "Complete Case", mu_cc, se_cc, boot_se_cc, boot_se_trim_cc))

cat(sprintf("  %-20s %-18s %10.4f %12s %12s %12s\n",
            "MAR", "IPW (MAR)", mu_mar, "-", "-", "-"))

for (j in seq_along(all_model_sets)) {
  label <- model_set_labels[j]
  name <- model_set_names[j]
  res <- results_list[[label]]
  cat(sprintf("  %-20s %-18s %10.4f %12.4f %12.4f %12.4f\n",
              label, name, res$mu, res$se, boot_se[label], boot_se_trim[label]))
}

cat("  ", paste(rep("-", 88), collapse = ""), "\n")

#------------------------------------------------------------------------------#
# Sensitivity Analysis: Perturb Model 1 with Exponential Tilt
#------------------------------------------------------------------------------#
cat("\n================================================================================\n")
cat("                         SENSITIVITY ANALYSIS                                  \n")
cat("================================================================================\n\n")

n_xi <- 30
xi_max <- 10 * n^(-1/2)
xi_values <- seq(0, xi_max, length.out = n_xi)

cat("Sensitivity analysis parameters:\n")
cat("  xi range: [0, ", round(xi_max, 6), "]\n", sep = "")
cat("  Number of xi values:", n_xi, "\n")
cat("  Perturbing model: 1\n")
cat("  Exponential tilt variables: teacher_report, father\n\n")

baseline_ps_matrix <- do.call(cbind, lapply(ebmr$ps_fit.list, function(pf) pf$fitted.values))
exp_tilt_x_names <- c("father")

sensitivity_results <- data.frame(xi = xi_values, mu_ipw = NA)

cat("Running sensitivity analysis:\n")
pb_sa <- txtProgressBar(min = 0, max = n_xi, style = 3)

for (i in 1:n_xi) {
  xi <- xi_values[i]

  exp_tilt <- function(y, x) {
    exp(xi * as.matrix(cbind(y, x)) %*% c(1, 1))
  }

  sa_result <- tryCatch({
    ebmr$EBMR_IPW_with_locally_misspecified_model(
      ps.matrix = baseline_ps_matrix,
      perturb_ps = 1,
      exp_tilt = exp_tilt,
      exp_tilt_x_names = exp_tilt_x_names,
      h_nu = h_nu,
      nu_init = result$nu.hat,
      se.fit = FALSE
    )
  }, error = function(e) NULL)

  if (!is.null(sa_result)) {
    sensitivity_results$mu_ipw[i] <- sa_result$mu_ipw
  }

  setTxtProgressBar(pb_sa, i)
}
close(pb_sa)

cat("\n\nSensitivity Analysis Results:\n")
cat(sprintf("  %-10s %12s\n", "xi", "mu_ipw"))
cat("  ", paste(rep("-", 24), collapse = ""), "\n")
for (i in seq(1, n_xi, by = 5)) {
  cat(sprintf("  %-10.6f %12.4f\n", sensitivity_results$xi[i], sensitivity_results$mu_ipw[i]))
}
cat("  ... (showing every 5th value)\n\n")

cat("Summary:\n")
cat("  mu_ipw range: [", round(min(sensitivity_results$mu_ipw, na.rm = TRUE), 4),
    ", ", round(max(sensitivity_results$mu_ipw, na.rm = TRUE), 4), "]\n", sep = "")
cat("  Baseline (xi=0):", round(sensitivity_results$mu_ipw[1], 4), "\n")
cat("  At xi_max:", round(sensitivity_results$mu_ipw[n_xi], 4), "\n")
cat("  Change from baseline:", round(sensitivity_results$mu_ipw[n_xi] - sensitivity_results$mu_ipw[1], 4), "\n")

png("MHD_results/health1_sensitivity_analysis_Fast4.png", width = 800, height = 600)
plot(sensitivity_results$xi, sensitivity_results$mu_ipw,
     type = "b", pch = 19, col = "blue",
     xlab = expression(xi),
     ylab = expression(hat(mu)[IPW]),
     main = "Sensitivity Analysis: Health=1 Subgroup Mean (Fast4)\nPerturbing Model 1 with Exponential Tilt",
     ylim = range(c(sensitivity_results$mu_ipw, mu_cc), na.rm = TRUE))
abline(h = mu_cc, col = "red", lty = 2, lwd = 2)
abline(h = result$mu_ipw, col = "darkgreen", lty = 3, lwd = 2)
legend("topright",
       legend = c("Sensitivity Analysis", "Complete Case", "EBMR (xi=0)"),
       col = c("blue", "red", "darkgreen"),
       lty = c(1, 2, 3), pch = c(19, NA, NA), lwd = 2)
dev.off()
cat("  Saved: MHD_results/health1_sensitivity_analysis_Fast4.png\n")

cat("\n================================================================================\n")
