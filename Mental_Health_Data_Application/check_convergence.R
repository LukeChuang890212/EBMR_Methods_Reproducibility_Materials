#------------------------------------------------------------------------------#
# Check Convergence of PS Models for Each Health Subgroup
#------------------------------------------------------------------------------#

library(EBMRalgorithmFast)
library(tidyverse)

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
# PS specifications
#------------------------------------------------------------------------------#
full_ps_specifications <- list(
  formula.list = list(
    r ~ o(teacher_report) + father,
    r ~ o(teacher_report) + parent_report,
    r ~ father + parent_report  # MAR model (no o(teacher_report))
  ),
  h_x_names.list = list(
    c("father", "parent_report"),
    c("father", "parent_report"),
    c("father", "parent_report")
  ),
  inv_link = function(eta) 1 / (1 + exp(eta))
)

W <- function(g.matrix) {
  solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
}

init.list <- list(
  alpha_init.list0 = list(
    c(-0.3834341, 1.07945126, -0.1318883),
    c(-4.2670736, 5.9582123, -0.5972129),
    c(0.5, -0.2, -0.5)  # MAR model: intercept, father, parent_report
  ),
  alpha_init.list1 = list(
    c(0.09019492, -1.511548, -0.3946616),
    c(-1.969056, 3.5823611, -1.37232734),
    c(0.5, -0.2, -0.5)  # MAR model: intercept, father, parent_report
  )
)

health <- 0:1
subset_names <- c("health0", "health1")
model_names <- c("Model 1 (MNAR: father)", "Model 2 (MNAR: parent_report)", "Model 3 (MAR)")

#------------------------------------------------------------------------------#
# Fit models and check convergence
#------------------------------------------------------------------------------#
cat("================================================================================\n")
cat("           CONVERGENCE CHECK FOR PS MODELS BY HEALTH SUBGROUP                   \n")
cat("================================================================================\n\n")

for (k in 1:length(health)) {
  subdat <- dat[dat$health == health[k], ]
  subdat$y <- subdat$teacher_report

  cat("--------------------------------------------------------------------------------\n")
  cat("Health =", health[k], "(", ifelse(k == 1, "Unexposed", "Exposed"), ")\n")
  cat("Sample size: n =", nrow(subdat), "\n")
  cat("Missing rate:", round(100 * (1 - mean(subdat$r)), 1), "%\n")
  cat("--------------------------------------------------------------------------------\n\n")

  ps_specifications <- list(
    formula.list = full_ps_specifications$formula.list[1:3],
    h_x_names.list = full_ps_specifications$h_x_names.list[1:3],
    alpha_init.list = init.list[[k]][1:3],
    inv_link = full_ps_specifications$inv_link
  )

  ebmr <- EBMRAlgorithmFast$new("teacher_report", ps_specifications, subdat, W)

  for (j in 1:3) {
    ps_fit <- ebmr$ps_fit.list[[j]]
    gmm_fit <- ps_fit$gmm_fit

    cat("  ", model_names[j], "\n")
    cat("  ", paste(rep("-", 60), collapse = ""), "\n")

    # Check if gmm_fit has opt field (MNAR models)
    if (!is.null(gmm_fit$opt)) {
      opt <- gmm_fit$opt

      # optim convergence codes:
      # 0: successful convergence
      # 1: iteration limit reached
      # 10: degeneracy of Nelder-Mead simplex
      # 51: warning from L-BFGS-B
      # 52: error from L-BFGS-B
      conv_code <- opt$convergence
      conv_status <- switch(as.character(conv_code),
        "0" = "CONVERGED (successful)",
        "1" = "NOT CONVERGED (iteration limit reached)",
        "10" = "DEGENERATE (Nelder-Mead simplex)",
        "51" = "WARNING from L-BFGS-B",
        "52" = "ERROR from L-BFGS-B",
        paste("Unknown code:", conv_code)
      )

      cat("    Convergence code:", conv_code, "-", conv_status, "\n")
      cat("    Final objective value:", format(opt$value, scientific = FALSE, digits = 8), "\n")

      if (!is.null(opt$message) && nchar(opt$message) > 0) {
        cat("    Message:", opt$message, "\n")
      }

      cat("    Coefficients:\n")
      coef_names <- c("Intercept", "o(y)", "father/parent_report")
      if (j == 1) coef_names <- c("Intercept", "o(y)", "father")
      if (j == 2) coef_names <- c("Intercept", "o(y)", "parent_report")

      for (i in 1:length(ps_fit$coefficients)) {
        cat(sprintf("      %-20s: %10.6f", coef_names[i], ps_fit$coefficients[i]))
        if (!is.null(ps_fit$se) && !is.na(ps_fit$se[i])) {
          cat(sprintf("  (SE: %8.6f)", ps_fit$se[i]))
        }
        cat("\n")
      }

      # Check propensity score range
      ps_range <- range(ps_fit$fitted.values)
      cat("    PS range: [", round(ps_range[1], 4), ",", round(ps_range[2], 4), "]\n")

      # Check for extreme propensity scores
      extreme_low <- sum(ps_fit$fitted.values < 0.01)
      extreme_high <- sum(ps_fit$fitted.values > 0.99)
      if (extreme_low > 0 || extreme_high > 0) {
        cat("    WARNING: Extreme PS values - ", extreme_low, " below 0.01, ",
            extreme_high, " above 0.99\n")
      }

    } else {
      # MAR model fitted with glm
      cat("    (Fitted using glm - standard convergence)\n")
      cat("    Coefficients:\n")
      coef_names <- c("Intercept", "father", "parent_report")
      for (i in 1:length(ps_fit$coefficients)) {
        cat(sprintf("      %-20s: %10.6f", coef_names[i], ps_fit$coefficients[i]))
        if (!is.null(ps_fit$se) && !is.na(ps_fit$se[i])) {
          cat(sprintf("  (SE: %8.6f)", ps_fit$se[i]))
        }
        cat("\n")
      }

      # Check propensity score range
      ps_range <- range(ps_fit$fitted.values)
      cat("    PS range: [", round(ps_range[1], 4), ",", round(ps_range[2], 4), "]\n")
    }

    cat("\n")
  }

  # Now run ensemble step and check its convergence
  h_nu <- function(data) cbind(f = data$father, p = data$parent_report, fp = data$father * data$parent_report)
  result <- ebmr$EBMR_IPW(h_nu = h_nu, true_ps = NULL)

  cat("  Ensemble Step\n")
  cat("  ", paste(rep("-", 60), collapse = ""), "\n")

  ensemble_gmm <- result$ensemble_fit$gmm_fit
  if (!is.null(ensemble_gmm$opt)) {
    opt <- ensemble_gmm$opt
    conv_code <- opt$convergence
    conv_status <- switch(as.character(conv_code),
      "0" = "CONVERGED (successful)",
      "1" = "NOT CONVERGED (iteration limit reached)",
      "10" = "DEGENERATE (Nelder-Mead simplex)",
      "51" = "WARNING from L-BFGS-B",
      "52" = "ERROR from L-BFGS-B",
      paste("Unknown code:", conv_code)
    )

    cat("    Convergence code:", conv_code, "-", conv_status, "\n")
    cat("    Final objective value:", format(opt$value, scientific = FALSE, digits = 8), "\n")

    if (!is.null(opt$message) && nchar(opt$message) > 0) {
      cat("    Message:", opt$message, "\n")
    }
  }

  cat("    nu.hat:", round(result$nu.hat, 4), "\n")
  cat("    w.hat:", round(result$w.hat, 4), "(", round(result$w.hat * 100, 1), "%)\n")
  cat("    mu_ipw:", round(result$mu_ipw, 4), "(SE:", round(result$se_ipw, 4), ")\n")

  cat("\n")
}

cat("================================================================================\n")
cat("                         CONVERGENCE CHECK COMPLETE                             \n")
cat("================================================================================\n")
