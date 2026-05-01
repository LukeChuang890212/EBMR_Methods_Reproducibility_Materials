#------------------------------------------------------------------------------#
# EBMR Regression Simulation Functions
#
# This file contains functions for:
#   1. Running regression simulations (simulate_regression, run_scenario_regression)
#   2. Summarizing results (summarize_regression_results, print_regression_summary)
#   3. Utility functions (clean_regression_result, get_theta_true)
#
# For estimating regression coefficients under MNAR using EBMR_IPW_regression.
#------------------------------------------------------------------------------#

# Load required packages for summarization
library(knitr)
library(kableExtra)

#------------------------------------------------------------------------------#
# True Regression Coefficients for Each Setting
#------------------------------------------------------------------------------#
THETA_TRUE <- list(
  # Setting 13: Continuous y (Gaussian family)
  # y ~ N(0.2 + 0.6*z1 + 0.6*z2 + 0.3*u1 + 0.3*u2, 1)
  setting13 = c(intercept = 0.2, z1 = 0.6, z2 = 0.6, u1 = 0.3, u2 = 0.3),

  # Setting 14: Binary y (Binomial family)
  # y ~ Bernoulli(expit(0.2 + 0.8*z1 + 0.8*z2 + 0.4*u1 + 0.4*u2))
  setting14 = c(intercept = 0.2, z1 = 0.8, z2 = 0.8, u1 = 0.4, u2 = 0.4)
)

#------------------------------------------------------------------------------#
# Configuration for Regression Simulations
#------------------------------------------------------------------------------#
CONFIG_REGRESSION <- list(
  replicate_num = 1000,
  n_default = 2000
)

#------------------------------------------------------------------------------#
# Sample Size Vectors for Regression
#------------------------------------------------------------------------------#
n.vector.list.regression <- list(
  correct_model = list(c(2000, 500)),
  misspecified_model = list(c(2000), c(500))
)

#------------------------------------------------------------------------------#
# Data File Lists for Regression (Settings 13/14)
#------------------------------------------------------------------------------#
data_root_regression <- "Simulation_Data/"

correct_model_all_data_file.list.regression <- list(
  setting13 = list(
    miss50 = list(
      paste0(data_root_regression, "setting13.A1_n2000_replicate1000.RDS")
    ),
    miss30 = list(
      paste0(data_root_regression, "setting13.A2_n2000_replicate1000.RDS")
    )
  ),
  setting14 = list(
    miss50 = list(
      paste0(data_root_regression, "setting14.A1_n2000_replicate1000.RDS")
    ),
    miss30 = list(
      paste0(data_root_regression, "setting14.A2_n2000_replicate1000.RDS")
    )
  )
)

misspecified_model_all_data_file.list.regression <- list(
  setting13 = list(
    miss50 = list(
      paste0(data_root_regression, "setting13.B1_n2000_replicate1000.RDS"),
      paste0(data_root_regression, "setting13.B1_n500_replicate1000.RDS")
    ),
    miss30 = list(
      paste0(data_root_regression, "setting13.B2_n2000_replicate1000.RDS"),
      paste0(data_root_regression, "setting13.B2_n500_replicate1000.RDS")
    )
  ),
  setting14 = list(
    miss50 = list(
      paste0(data_root_regression, "setting14.B1_n2000_replicate1000.RDS"),
      paste0(data_root_regression, "setting14.B1_n500_replicate1000.RDS")
    ),
    miss30 = list(
      paste0(data_root_regression, "setting14.B2_n2000_replicate1000.RDS"),
      paste0(data_root_regression, "setting14.B2_n500_replicate1000.RDS")
    )
  )
)

#------------------------------------------------------------------------------#
# Alpha True Lists for Regression (Same as population mean estimation)
#------------------------------------------------------------------------------#
correct_model_alpha.true.list.regression <- list(
  setting13 = list(
    miss50 = list(c(0.1169, 0.2, -0.4, -0.4)),
    miss30 = list(c(-0.7772, 0.2, -0.4, -0.4))
  ),
  setting14 = list(
    miss50 = list(c(0.1154, 0.2, -0.4, -0.4)),
    miss30 = list(c(-0.7696, 0.2, -0.4, -0.4))
  )
)

misspecified_model_alpha.true.list.regression <- list(
  setting13 = list(
    miss50 = list(
      c(0.1169, 0.2, -0.4, -0.4),
      c(0.1169, 0.2, -0.4, -0.4)
    ),
    miss30 = list(
      c(-0.7772, 0.2, -0.4, -0.4),
      c(-0.7772, 0.2, -0.4, -0.4)
    )
  ),
  setting14 = list(
    miss50 = list(
      c(0.1154, 0.2, -0.4, -0.4),
      c(0.1154, 0.2, -0.4, -0.4)
    ),
    miss30 = list(
      c(-0.7696, 0.2, -0.4, -0.4),
      c(-0.7696, 0.2, -0.4, -0.4)
    )
  )
)

#------------------------------------------------------------------------------#
# PS Specifications for Regression
# Matching scenarios from population mean estimation (7-1, 8-1, 8-2, 9-1)
#------------------------------------------------------------------------------#

# Inverse link function (complement logistic)
INV_LINK_REGRESSION <- function(eta) 1 / (1 + exp(eta))

# PS Specification factory
create_ps_spec_regression <- function(formulas, h_alpha_list, outcome = "y") {
  list(
    formula.list = formulas,
    h_alpha.list = h_alpha_list,
    inv_link = INV_LINK_REGRESSION,
    outcome = outcome
  )
}

# PS Specifications registry (matching population mean scenarios)
PS_SPECS_REGRESSION <- list(
  # Scenario 7: full, u1_z1, u2_z1
  `7` = create_ps_spec_regression(
    formulas = list(
      r ~ y + u1 + u2,    # Full model
      r ~ y + u1 + z1,    # u1 + z1
      r ~ y + u2 + z1     # u2 + z1
    ),
    h_alpha_list = list(
      c("u1", "u2", "z1", "z2"),
      c("u1", "u2", "z1", "z2"),
      c("u1", "u2", "z1", "z2")
    )
  ),

  # Scenario 8: full, u1_z1, u2_z2
  `8` = create_ps_spec_regression(
    formulas = list(
      r ~ y + u1 + u2,    # Full model
      r ~ y + u1 + z1,    # u1 + z1
      r ~ y + u2 + z2     # u2 + z2
    ),
    h_alpha_list = list(
      c("u1", "u2", "z1", "z2"),
      c("u1", "u2", "z1", "z2"),
      c("u1", "u2", "z1", "z2")
    )
  ),

  # Scenario 9: full, u1_z2, u2_z2
  `9` = create_ps_spec_regression(
    formulas = list(
      r ~ y + u1 + u2,    # Full model
      r ~ y + u1 + z2,    # u1 + z2
      r ~ y + u2 + z2     # u2 + z2
    ),
    h_alpha_list = list(
      c("u1", "u2", "z1", "z2"),
      c("u1", "u2", "z1", "z2"),
      c("u1", "u2", "z1", "z2")
    )
  )
)

# Default PS specification (for backwards compatibility)
PS_SPEC_REGRESSION <- PS_SPECS_REGRESSION[["7"]]

#------------------------------------------------------------------------------#
# Utility Functions
#------------------------------------------------------------------------------#

#' Get true regression coefficients for a setting
#'
#' @param setting Setting name (e.g., "setting13", "setting14")
#' @return Named vector of true regression coefficients
get_theta_true <- function(setting) {
  if (!setting %in% names(THETA_TRUE)) {
    stop(paste("Unknown setting:", setting,
               "\nAvailable:", paste(names(THETA_TRUE), collapse = ", ")))
  }
  THETA_TRUE[[setting]]
}

#' Get family for a setting
#'
#' @param setting Setting name
#' @return "gaussian" or "binomial"
get_family <- function(setting) {
  if (setting == "setting13") {
    "gaussian"
  } else if (setting == "setting14") {
    "binomial"
  } else {
    stop(paste("Unknown setting:", setting))
  }
}

#' h_nu function for ensemble step
#' @param data Data frame
#' @return Matrix of auxiliary variables
h_nu_regression <- function(data) {
  cbind(u1 = data$u1, u2 = data$u2, z1 = data$z1, z2 = data$z2)
}

#' Weight matrix function for GMM
#' @param g.matrix Matrix of moment conditions
#' @return Weight matrix
W_regression <- function(g.matrix) {
  solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
}

#' Clean regression simulation results by removing NAs and outliers
#'
#' @param sim_result Simulation result matrix (rows are statistics, columns are replicates)
#' @param n_params Number of regression parameters (default: 5)
#' @param multiplier IQR multiplier for outlier detection (default: 3)
#' @param verbose Whether to print summary messages (default: TRUE)
#' @return List with cleaned result and counts
clean_regression_result <- function(sim_result, n_params = 5, multiplier = 3, verbose = TRUE) {
  n_total <- ncol(sim_result)
  n_na <- 0
  n_outliers <- 0

  # Row indices:
  # EBMR: theta (1:5), se (6:10)
  # resp_rate (11), w1 (12)
  # True PS: theta (13:17), se (18:22)
  theta_rows <- 1:n_params
  se_rows <- (n_params + 1):(2 * n_params)

  # Remove NAs (check theta and se rows)
  check_rows <- c(theta_rows, se_rows)
  na_mask <- apply(sim_result[check_rows, , drop = FALSE], 2, function(x) !any(is.na(x)))
  n_na <- sum(!na_mask)
  sim_result <- sim_result[, na_mask, drop = FALSE]

  # Remove outliers (check theta rows using IQR method)
  if (ncol(sim_result) > 0) {
    is_outlier <- rep(FALSE, ncol(sim_result))
    for (row_idx in theta_rows) {
      Q1 <- quantile(sim_result[row_idx, ], 0.25)
      Q3 <- quantile(sim_result[row_idx, ], 0.75)
      IQR_value <- Q3 - Q1
      lower_bound <- Q1 - multiplier * IQR_value
      upper_bound <- Q3 + multiplier * IQR_value
      is_outlier <- is_outlier | sim_result[row_idx, ] > upper_bound |
        sim_result[row_idx, ] < lower_bound
    }
    n_outliers <- sum(is_outlier)
    sim_result <- sim_result[, !is_outlier, drop = FALSE]
  }

  n_successful <- ncol(sim_result)

  if (verbose && (n_na > 0 || n_outliers > 0)) {
    cat("Replicates:", n_total, "-> Used:", n_successful,
        "(NA:", n_na, ", Outliers:", n_outliers, ")\n")
  }

  list(
    result = sim_result,
    n_total = n_total,
    n_na = n_na,
    n_outliers = n_outliers,
    n_successful = n_successful
  )
}

#------------------------------------------------------------------------------#
# Core Simulation Function
#------------------------------------------------------------------------------#

#' Run regression simulation for a single configuration
#'
#' @param all_data Full dataset for all replicates
#' @param ps_model.true Function to compute true propensity scores
#' @param alpha.true True alpha parameter values for PS model
#' @param ps_specifications List of PS model specifications
#' @param n Sample size per replicate
#' @param replicate_num Number of replicates to run
#' @param setting Setting name for determining family
#' @param save_file Path to save results (optional)
#' @return Simulation results matrix (invisibly)
simulate_regression <- function(all_data, ps_model.true, alpha.true, ps_specifications,
                                 n, replicate_num, setting, save_file = NULL) {

  library(EBMRalgorithmFast2)
  library(parallel)
  library(foreach)
  library(doSNOW)

  theta.true <- get_theta_true(setting)
  family <- get_family(setting)
  n_params <- length(theta.true)

  # Setup parallel cluster
  cores <- detectCores()
  cl <- makeCluster(cores - 2)
  registerDoSNOW(cl)

  # Export local variables to workers
  clusterExport(cl, c("ps_specifications", "n", "all_data", "theta.true", "family",
                      "h_nu_regression", "W_regression", "ps_model.true", "alpha.true"),
                envir = environment())

  # Progress bar
  pb <- txtProgressBar(max = replicate_num, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  parallel_packages <- c("EBMRalgorithmFast2", "stringr", "Matrix", "numDeriv")

  start <- Sys.time()

  sim_result <- foreach(
    i = 1:replicate_num,
    .combine = 'cbind',
    .options.snow = opts,
    .packages = parallel_packages
  ) %dopar% {
    tryCatch({
      # Extract data for this replicate
      dat <- all_data[((i - 1) * n + 1):(i * n), ]

      # Initialize result vector
      # EBMR: theta (1:5) + se (6:10) = 10
      # True PS: theta (13:17) + se (18:22) = 10
      # response_rate (11) + w.hat[1] (12)
      result <- rep(NA, 22)

      result[11] <- mean(dat$r)  # Response rate

      #-----------------------------------------------------------------------#
      # 1. True PS IPW Regression (Oracle)
      #-----------------------------------------------------------------------#
      tryCatch({
        # Compute true propensity scores
        pi_true <- ps_model.true(dat, alpha.true)

        # Compute Hajek weights for responders
        resp_idx <- dat$r == 1
        w_true <- rep(0, nrow(dat))
        w_true[resp_idx] <- (1 / pi_true[resp_idx]) / sum(1 / pi_true[resp_idx])

        # Fit weighted regression
        dat_resp <- dat[resp_idx, ]
        w_resp <- w_true[resp_idx] * sum(resp_idx)  # Scale weights for glm

        if (family == "binomial") {
          fit_true <- glm(y ~ z1 + z2 + u1 + u2, data = dat_resp,
                          family = binomial, weights = w_resp)
        } else {
          fit_true <- glm(y ~ z1 + z2 + u1 + u2, data = dat_resp,
                          family = gaussian, weights = w_resp)
        }

        theta_true_ps <- coef(fit_true)
        result[13:17] <- theta_true_ps

        # Sandwich SE for Hajek estimator
        X <- model.matrix(fit_true)
        n_resp <- nrow(dat_resp)
        p <- ncol(X)

        if (family == "binomial") {
          mu_hat <- fitted(fit_true)
          V <- diag(as.vector(mu_hat * (1 - mu_hat)))
          residuals_raw <- dat_resp$y - mu_hat
        } else {
          sigma2 <- summary(fit_true)$dispersion
          V <- diag(rep(1, n_resp))
          residuals_raw <- dat_resp$y - fitted(fit_true)
        }

        # Bread: (X'WVX)^{-1}
        W_diag <- diag(w_resp)
        bread <- solve(t(X) %*% W_diag %*% V %*% X)

        # Meat: sum of weighted squared score contributions
        # Score contribution for each observation
        score_i <- X * as.vector(w_resp * residuals_raw)
        meat <- t(score_i) %*% score_i

        # Sandwich variance
        var_sandwich <- bread %*% meat %*% bread
        se_true_ps <- sqrt(diag(var_sandwich))
        result[18:22] <- se_true_ps

      }, error = function(e) {})

      #-----------------------------------------------------------------------#
      # 2. EBMR-IPW Regression
      #-----------------------------------------------------------------------#
      tryCatch({
        ebmr <- EBMRAlgorithmFast2$new("y", ps_specifications, dat, W_regression)

        ebmr_result <- ebmr$EBMR_IPW_regression(
          h_nu = h_nu_regression,
          reg_formula = y ~ z1 + z2 + u1 + u2,
          family = family,
          se.fit = TRUE
        )

        result[1:5] <- ebmr_result$theta_hat
        result[6:10] <- ebmr_result$se_theta
        result[12] <- ebmr_result$w.hat[1]

      }, error = function(e) {})

      result

    }, error = function(e) {
      cat("ERROR in replicate", i, ":", conditionMessage(e), "\n")
      return(rep(NA, 22))
    })
  }

  close(pb)
  stopCluster(cl)
  cat("Elapsed time:", round(difftime(Sys.time(), start, units = "mins"), 2), "minutes\n")

  # Add row names
  param_names <- c("intercept", "z1", "z2", "u1", "u2")
  rownames(sim_result) <- c(
    paste0("theta_", param_names),
    paste0("se_", param_names),
    "resp_rate", "w1",
    paste0("true_theta_", param_names),
    paste0("true_se_", param_names)
  )

  # Save results
  if (!is.null(save_file)) {
    saveRDS(sim_result, save_file)
    cat("Results saved to:", save_file, "\n")
  }

  gc()
  invisible(sim_result)
}

#------------------------------------------------------------------------------#
# Result Summarization Functions
#------------------------------------------------------------------------------#

#' Summarize regression results for a single parameter
#'
#' @param theta_vals Vector of parameter estimates
#' @param se_vals Vector of standard error estimates
#' @param theta_true True parameter value
#' @return Named vector with Bias, ESD, ESE, SE.Ratio, CP
summarize_param_results <- function(theta_vals, se_vals, theta_true) {
  # Remove NAs
  valid <- !is.na(theta_vals) & !is.na(se_vals) & se_vals > 0
  theta_vals <- theta_vals[valid]
  se_vals <- se_vals[valid]
  n_valid <- length(theta_vals)

  if (n_valid < 10) {
    return(c(n_valid = n_valid, bias = NA, esd = NA, ese = NA,
             se_ratio = NA, cp = NA))
  }

  # Bias
  bias <- mean(theta_vals) - theta_true

  # Empirical SE (ESD)
  esd <- sd(theta_vals)

  # Average estimated SE (ESE)
  ese <- mean(se_vals)

  # SE ratio
  se_ratio <- ese / esd

  # Coverage (95% CI)
  lower <- theta_vals - 1.96 * se_vals
  upper <- theta_vals + 1.96 * se_vals
  cp <- mean(lower <= theta_true & theta_true <= upper)

  c(n_valid = n_valid, bias = bias, esd = esd, ese = ese,
    se_ratio = se_ratio, cp = cp)
}

#' Summarize regression simulation results
#'
#' @param sim_result Simulation result matrix
#' @param setting Setting name
#' @param method "ebmr" or "true_ps"
#' @return Summary matrix
summarize_regression_results <- function(sim_result, setting, method = "ebmr") {
  theta.true <- get_theta_true(setting)
  n_params <- length(theta.true)
  param_names <- c("Intercept", "z1", "z2", "u1", "u2")

  if (method == "ebmr") {
    theta_rows <- 1:n_params
    se_rows <- (n_params + 1):(2 * n_params)
  } else if (method == "true_ps") {
    theta_rows <- 13:17
    se_rows <- 18:22
  } else {
    stop(paste("Unknown method:", method))
  }

  summary_mat <- t(sapply(1:n_params, function(j) {
    summarize_param_results(
      sim_result[theta_rows[j], ],
      sim_result[se_rows[j], ],
      theta.true[j]
    )
  }))
  rownames(summary_mat) <- param_names

  # Ensure consistent column names (sapply may mangle them when combining named vectors)
  colnames(summary_mat) <- c("n_valid", "bias", "esd", "ese", "se_ratio", "cp")
  summary_mat
}

#' Print formatted regression summary
#'
#' @param sim_result Simulation result matrix
#' @param setting Setting name
#' @param n Sample size
#' @param version Version string
print_regression_summary <- function(sim_result, setting, n, version = NULL) {
  theta.true <- get_theta_true(setting)
  family <- get_family(setting)

  cat("\n")
  cat(strrep("=", 80), "\n")
  cat("REGRESSION SIMULATION SUMMARY\n")
  cat(strrep("=", 80), "\n")
  cat("Setting:", setting, "| Family:", family, "| n =", n, "\n")
  if (!is.null(version)) cat("Version:", version, "\n")
  cat("True theta:", paste(round(theta.true, 4), collapse = ", "), "\n")
  cat("Replicates:", ncol(sim_result), "\n")
  cat("Response rate (mean):", round(mean(sim_result["resp_rate", ], na.rm = TRUE), 3), "\n")
  cat(strrep("-", 80), "\n\n")

  # True PS Summary (Oracle)
  cat("True PS IPW Regression (Oracle):\n")
  cat(strrep("-", 80), "\n")
  summary_true <- summarize_regression_results(sim_result, setting, "true_ps")
  print(round(summary_true, 4))

  cat("\n")

  # EBMR Summary
  cat("EBMR-IPW Regression:\n")
  cat(strrep("-", 80), "\n")
  summary_ebmr <- summarize_regression_results(sim_result, setting, "ebmr")
  print(round(summary_ebmr, 4))

  cat(strrep("=", 80), "\n\n")

  invisible(list(true_ps = summary_true, ebmr = summary_ebmr))
}

#------------------------------------------------------------------------------#
# Scenario Management for Regression
#------------------------------------------------------------------------------#

#' Get model parameters for regression based on model type
#'
#' @param model_type Either "correct" or "misspecified"
#' @return List with n_vector, data_files, and alpha_true
get_model_params_regression <- function(model_type) {
  if (model_type == "correct") {
    list(
      n_vector = n.vector.list.regression$correct_model,
      data_files = correct_model_all_data_file.list.regression,
      alpha_true = correct_model_alpha.true.list.regression
    )
  } else if (model_type == "misspecified") {
    list(
      n_vector = n.vector.list.regression$misspecified_model,
      data_files = misspecified_model_all_data_file.list.regression,
      alpha_true = misspecified_model_alpha.true.list.regression
    )
  } else {
    stop(paste("Unknown model type:", model_type))
  }
}

#------------------------------------------------------------------------------#
# Equivalent Scenario Mapping for Regression
#
# When model combinations are identical between scenarios with the same data type
# (correct vs misspecified), results can be copied instead of re-running.
#
# PS Specifications:
#   7: full (r~y+u1+u2), u1_z1 (r~y+u1+z1), u2_z1 (r~y+u2+z1)
#   8: full (r~y+u1+u2), u1_z1 (r~y+u1+z1), u2_z2 (r~y+u2+z2)
#   9: full (r~y+u1+u2), u1_z2 (r~y+u1+z2), u2_z2 (r~y+u2+z2)
#
# Equivalence requires SAME DATA (correct=A data, misspecified=B data)
#------------------------------------------------------------------------------#

# Model 1 equivalence: all have same model 1 (full)
# Correct (A data): 7-1, 8-1, 9-1
# Misspecified (B data): 7-2, 8-2, 9-2
EQUIVALENT_SCENARIOS_MODEL1_REG <- list(
  `7-1` = c("8-1", "9-1"),
  `8-1` = c("7-1", "9-1"),
  `9-1` = c("7-1", "8-1"),
  `7-2` = c("8-2", "9-2"),
  `8-2` = c("7-2", "9-2"),
  `9-2` = c("7-2", "8-2")
)

# Model 2 equivalence: 7 and 8 have same model 2 (u1_z1)
# Correct: 7-1 = 8-1
# Misspecified: 7-2 = 8-2
EQUIVALENT_SCENARIOS_MODEL2_REG <- list(
  `7-1` = c("8-1"),
  `8-1` = c("7-1"),
  `7-2` = c("8-2"),
  `8-2` = c("7-2")
)

# Model 3 equivalence: 8 and 9 have same model 3 (u2_z2)
# Correct: 8-1 = 9-1
# Misspecified: 8-2 = 9-2
EQUIVALENT_SCENARIOS_MODEL3_REG <- list(
  `8-1` = c("9-1"),
  `9-1` = c("8-1"),
  `8-2` = c("9-2"),
  `9-2` = c("8-2")
)

# Model 12 equivalence: same model 1 and model 2
# Correct: 7-1 = 8-1
# Misspecified: 7-2 = 8-2
EQUIVALENT_SCENARIOS_MODEL12_REG <- list(
  `7-1` = c("8-1"),
  `8-1` = c("7-1"),
  `7-2` = c("8-2"),
  `8-2` = c("7-2")
)

# Model 13 equivalence: same model 1 and model 3
# Correct: 8-1 = 9-1
# Misspecified: 8-2 = 9-2
EQUIVALENT_SCENARIOS_MODEL13_REG <- list(
  `8-1` = c("9-1"),
  `9-1` = c("8-1"),
  `8-2` = c("9-2"),
  `9-2` = c("8-2")
)

# Model 23 equivalence: none (all have different model 2 + model 3 combinations)
EQUIVALENT_SCENARIOS_MODEL23_REG <- list()

# Model 123 equivalence: none (all have different full model combinations)
EQUIVALENT_SCENARIOS_MODEL123_REG <- list()

#------------------------------------------------------------------------------#
# Scenario Registry for Regression
# Matching population mean scenarios: 7-1, 8-1, 8-2, 9-1
#------------------------------------------------------------------------------#
SCENARIOS_REGRESSION <- list(
  # Scenario 7-1: Correct model with z1 interactions (full, u1+z1, u2+z1)
  `7-1` = list(
    id = "7-1",
    description = "Correct: full, u1+z1, u2+z1",
    ps_spec_id = "7",
    settings = c("setting13", "setting14"),
    missing_rates = c("miss50", "miss30"),
    model_type = "correct"
  ),

  # Scenario 7-2: Misspecified model with z1 interactions
  `7-2` = list(
    id = "7-2",
    description = "Misspecified: full, u1+z1, u2+z1",
    ps_spec_id = "7",
    settings = c("setting13", "setting14"),
    missing_rates = c("miss50", "miss30"),
    model_type = "misspecified"
  ),

  # Scenario 8-1: Correct model with z1/z2 interactions (full, u1+z1, u2+z2)
  `8-1` = list(
    id = "8-1",
    description = "Correct: full, u1+z1, u2+z2",
    ps_spec_id = "8",
    settings = c("setting13", "setting14"),
    missing_rates = c("miss50", "miss30"),
    model_type = "correct"
  ),

  # Scenario 8-2: Misspecified model with z1/z2 interactions
  `8-2` = list(
    id = "8-2",
    description = "Misspecified: full, u1+z1, u2+z2",
    ps_spec_id = "8",
    settings = c("setting13", "setting14"),
    missing_rates = c("miss50", "miss30"),
    model_type = "misspecified"
  ),

  # Scenario 9-1: Correct model with z2 interactions (full, u1+z2, u2+z2)
  `9-1` = list(
    id = "9-1",
    description = "Correct: full, u1+z2, u2+z2",
    ps_spec_id = "9",
    settings = c("setting13", "setting14"),
    missing_rates = c("miss50", "miss30"),
    model_type = "correct"
  ),

  # Scenario 9-2: Misspecified model with z2 interactions
  `9-2` = list(
    id = "9-2",
    description = "Misspecified: full, u1+z2, u2+z2",
    ps_spec_id = "9",
    settings = c("setting13", "setting14"),
    missing_rates = c("miss50", "miss30"),
    model_type = "misspecified"
  )
)

#' Get PS specification for a scenario
#'
#' @param scenario_id Scenario ID
#' @return PS specification list
get_ps_spec_regression <- function(scenario_id) {
  scenario <- SCENARIOS_REGRESSION[[scenario_id]]
  if (is.null(scenario)) {
    stop(paste("Unknown scenario:", scenario_id))
  }
  ps_spec_id <- scenario$ps_spec_id
  if (is.null(ps_spec_id)) {
    # Default to "7" for backwards compatibility
    ps_spec_id <- "7"
  }
  PS_SPECS_REGRESSION[[ps_spec_id]]
}

#' Run a regression scenario
#'
#' @param scenario_id The scenario ID to run (e.g., "reg-1", "reg-2")
#' @param setting Which setting to run (e.g., "setting13", "setting14")
#' @param n Sample size (optional - if NULL, runs all sizes from n.vector.list)
#' @param replicate_num Number of replicates
#' @param version Version string for output files
#' @param missing_rates Which missing rates to run
run_scenario_regression <- function(scenario_id,
                                     setting = NULL,
                                     n = NULL,
                                     replicate_num = NULL,
                                     version = "v1",
                                     missing_rates = NULL) {

  # Get scenario configuration
  if (!scenario_id %in% names(SCENARIOS_REGRESSION)) {
    stop(paste("Unknown scenario:", scenario_id,
               "\nAvailable:", paste(names(SCENARIOS_REGRESSION), collapse = ", ")))
  }
  scenario <- SCENARIOS_REGRESSION[[scenario_id]]

  # Use defaults if not provided
  if (is.null(setting)) setting <- scenario$settings[1]
  if (is.null(replicate_num)) {
    replicate_num <- if (!is.null(CONFIG_REGRESSION$replicate_num)) {
      CONFIG_REGRESSION$replicate_num
    } else {
      1000
    }
  }
  if (is.null(missing_rates)) missing_rates <- scenario$missing_rates

  # Get model parameters
  params <- get_model_params_regression(scenario$model_type)

  # Get PS specification for this scenario
  ps_spec <- get_ps_spec_regression(scenario_id)

  # Get n.vector from params
  n_vector_list <- params$n_vector

  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("Running Regression Scenario:", scenario_id, "\n")
  cat("Description:", scenario$description, "\n")
  cat("PS Spec:", scenario$ps_spec_id, "\n")
  cat(strrep("=", 70), "\n")
  cat("  Setting:", setting, "\n")
  cat("  Family:", get_family(setting), "\n")
  cat("  True theta:", paste(round(get_theta_true(setting), 4), collapse = ", "), "\n")
  cat("  Replicates:", replicate_num, "\n")
  cat("  Version:", version, "\n")
  cat("  Missing rates:", paste(missing_rates, collapse = ", "), "\n")
  cat("\n")

  # Get data file config for this setting
  data_file_config <- params$data_files[[setting]]
  if (is.null(data_file_config)) {
    stop(paste("No data file configuration for setting:", setting))
  }
  alpha.true_config <- params$alpha_true[[setting]]

  # Create output directory
  if (!dir.exists("Simulation_Results")) {
    dir.create("Simulation_Results")
  }

  # Define true PS model function for setting13/14
  ps_model.true <- function(dat, alpha.true) {
    1 / (1 + exp(cbind(rep(1, nrow(dat)), dat$y, dat$u1, dat$u2) %*% alpha.true))
  }

  # Run for each missing rate
  for (miss_rate in missing_rates) {
    cat("Processing:", miss_rate, "\n")

    # Get data file path(s)
    data_files_for_rate <- data_file_config[[miss_rate]]
    if (is.null(data_files_for_rate)) {
      cat("  WARNING: No data file config for", miss_rate, "\n")
      next
    }

    # Get alpha.true for this missing rate
    alpha.true_for_rate <- alpha.true_config[[miss_rate]]

    # Iterate over data files
    for (file_idx in seq_along(data_files_for_rate)) {
      data_file <- data_files_for_rate[[file_idx]]
      alpha.true <- alpha.true_for_rate[[file_idx]]

      # Get n.vector for this file index
      n_vector <- n_vector_list[[file_idx]]

      # If user specified n, filter to just that n
      if (!is.null(n)) {
        if (n %in% n_vector) {
          n_vector <- n
        } else {
          next
        }
      }

      if (!file.exists(data_file)) {
        cat("  WARNING: Data file not found:", data_file, "\n")
        cat("  Use generate_data() to create it first.\n")
        next
      }

      # Load data
      all_data <- readRDS(data_file)
      cat("  Loaded data:", data_file, "\n")

      # Run for each sample size
      for (current_n in n_vector) {
        cat("  Sample size n =", current_n, "\n")

        # Run simulation for each model combination
        J <- length(ps_spec$formula.list)
        for (model_num in 1:J) {
          model_combinations <- combn(J, model_num)
          for (i in 1:ncol(model_combinations)) {
            model_set <- model_combinations[, i]
            model_str <- paste0(model_set, collapse = "")

            cat("    Running models:", model_str, "\n")

            # Create subset PS specification
            subset_ps_spec <- list(
              formula.list = ps_spec$formula.list[model_set],
              h_alpha.list = ps_spec$h_alpha.list[model_set],
              inv_link = ps_spec$inv_link,
              outcome = ps_spec$outcome
            )

            # Output file
            save_file <- paste0(
              "Simulation_Results/EBMR_IPW_regression_", setting, "-", miss_rate,
              "-scenario", scenario_id, "_", model_str,
              "_n", current_n, "_replicate", replicate_num, "_", version, ".RDS"
            )

            # Check if already exists
            if (file.exists(save_file)) {
              cat("      Already exists, skipping\n")
              next
            }

            # Check for equivalent results from other scenarios
            check_and_copy_equiv_reg <- function(equiv_list) {
              equiv_scenarios <- equiv_list[[scenario_id]]
              if (!is.null(equiv_scenarios)) {
                for (equiv_scen in equiv_scenarios) {
                  equiv_file <- paste0(
                    "Simulation_Results/EBMR_IPW_regression_", setting, "-", miss_rate,
                    "-scenario", equiv_scen, "_", model_str,
                    "_n", current_n, "_replicate", replicate_num, "_", version, ".RDS"
                  )
                  if (file.exists(equiv_file)) {
                    file.copy(equiv_file, save_file)
                    cat("      Copied from equivalent scenario:", equiv_scen, "\n")
                    return(TRUE)
                  }
                }
              }
              return(FALSE)
            }

            # Check equivalence based on model combination
            copied <- FALSE
            if (model_str == "1") {
              copied <- check_and_copy_equiv_reg(EQUIVALENT_SCENARIOS_MODEL1_REG)
            } else if (model_str == "2") {
              copied <- check_and_copy_equiv_reg(EQUIVALENT_SCENARIOS_MODEL2_REG)
            } else if (model_str == "3") {
              copied <- check_and_copy_equiv_reg(EQUIVALENT_SCENARIOS_MODEL3_REG)
            } else if (model_str == "12") {
              copied <- check_and_copy_equiv_reg(EQUIVALENT_SCENARIOS_MODEL12_REG)
            } else if (model_str == "13") {
              copied <- check_and_copy_equiv_reg(EQUIVALENT_SCENARIOS_MODEL13_REG)
            } else if (model_str == "23") {
              copied <- check_and_copy_equiv_reg(EQUIVALENT_SCENARIOS_MODEL23_REG)
            } else if (model_str == "123") {
              copied <- check_and_copy_equiv_reg(EQUIVALENT_SCENARIOS_MODEL123_REG)
            }
            if (copied) next

            # Run simulation
            simulate_regression(
              all_data = all_data,
              ps_model.true = ps_model.true,
              alpha.true = alpha.true,
              ps_specifications = subset_ps_spec,
              n = current_n,
              replicate_num = replicate_num,
              setting = setting,
              save_file = save_file
            )

            # Print detailed summary after simulation completes
            if (file.exists(save_file)) {
              sim_result <- readRDS(save_file)
              cleaned <- clean_regression_result(sim_result, multiplier = 3, verbose = FALSE)
              sim_clean <- cleaned$result

              cat("      Replicates: ", cleaned$n_successful, "/", cleaned$n_total,
                  " (NA:", cleaned$n_na, ", Outliers:", cleaned$n_outliers, ")\n", sep = "")

              theta.true <- get_theta_true(setting)
              param_names <- c("Int", "z1", "z2", "u1", "u2")

              # Print summary for True PS (Oracle)
              cat("\n      True PS IPW (Oracle):\n")
              cat("      ", sprintf("%-6s %8s %8s %8s %8s\n", "Param", "Bias", "ESD", "ASE", "CP"))
              for (j in 1:5) {
                theta_vals <- sim_clean[j + 12, ]
                se_vals <- sim_clean[j + 17, ]
                valid <- !is.na(theta_vals) & !is.na(se_vals) & se_vals > 0
                theta_v <- theta_vals[valid]
                se_v <- se_vals[valid]

                if (length(theta_v) > 10) {
                  bias <- mean(theta_v) - theta.true[j]
                  esd <- sd(theta_v)
                  ase <- mean(se_v)
                  lower <- theta_v - 1.96 * se_v
                  upper <- theta_v + 1.96 * se_v
                  cp <- mean(lower <= theta.true[j] & theta.true[j] <= upper)
                  cat("      ", sprintf("%-6s %8.4f %8.4f %8.4f %8.3f\n",
                      param_names[j], bias, esd, ase, cp))
                } else {
                  cat("      ", sprintf("%-6s %8s %8s %8s %8s\n",
                      param_names[j], "NA", "NA", "NA", "NA"))
                }
              }

              # Print summary for EBMR
              cat("\n      EBMR-IPW Regression:\n")
              cat("      ", sprintf("%-6s %8s %8s %8s %8s\n", "Param", "Bias", "ESD", "ASE", "CP"))
              for (j in 1:5) {
                theta_vals <- sim_clean[j, ]
                se_vals <- sim_clean[j + 5, ]
                valid <- !is.na(theta_vals) & !is.na(se_vals) & se_vals > 0
                theta_v <- theta_vals[valid]
                se_v <- se_vals[valid]

                bias <- mean(theta_v) - theta.true[j]
                esd <- sd(theta_v)
                ase <- mean(se_v)
                lower <- theta_v - 1.96 * se_v
                upper <- theta_v + 1.96 * se_v
                cp <- mean(lower <= theta.true[j] & theta.true[j] <= upper)

                cat("      ", sprintf("%-6s %8.4f %8.4f %8.4f %8.3f\n",
                    param_names[j], bias, esd, ase, cp))
              }
              cat("\n")
            }
          }
        }
      }
    }
  }

  cat("\nScenario", scenario_id, "completed!\n")
}

#------------------------------------------------------------------------------#
# Summary Functions for All Settings
#------------------------------------------------------------------------------#

#' Summarize all settings and missing rates for regression
#'
#' @param settings Vector of settings
#' @param missing_rates Vector of missing rates
#' @param scenario_id Scenario ID
#' @param n.vector Vector of sample sizes (e.g., c(500, 2000))
#' @param replicate_num Number of replicates
#' @param version Version string
summarize_all_regression <- function(settings, missing_rates, scenario_id, n.vector,
                                      replicate_num, version) {
  cat("\n")
  cat(strrep("=", 90), "\n")
  cat("REGRESSION SIMULATION SUMMARY: Scenario", scenario_id, "\n")
  cat(strrep("=", 90), "\n\n")

  for (setting in settings) {
    theta.true <- get_theta_true(setting)
    family <- get_family(setting)

    cat(strrep("-", 90), "\n")
    cat("Setting:", setting, "| Family:", family, "\n")
    cat("True theta:", paste(round(theta.true, 4), collapse = ", "), "\n")
    cat(strrep("-", 90), "\n\n")

    for (miss in missing_rates) {
      cat("  Missing rate:", miss, "\n\n")

      for (current_n in n.vector) {
        cat("    n =", current_n, "\n")

        for (model_str in c("1", "2", "3", "12", "13", "23", "123")) {
          file <- paste0(
            "Simulation_Results/EBMR_IPW_regression_", setting, "-", miss,
            "-scenario", scenario_id, "_", model_str,
            "_n", current_n, "_replicate", replicate_num, "_", version, ".RDS"
          )

          if (file.exists(file)) {
            sim_result <- readRDS(file)

            cat("      Model:", model_str, "\n")

            # True PS summary
            summary_true <- summarize_regression_results(sim_result, setting, "true_ps")
            cat("        True - Bias:", paste(round(summary_true[, "bias"], 4), collapse = ", "),
                "| CP:", paste(round(summary_true[, "cp"], 3), collapse = ", "), "\n")

            # EBMR summary
            summary_ebmr <- summarize_regression_results(sim_result, setting, "ebmr")
            cat("        EBMR - Bias:", paste(round(summary_ebmr[, "bias"], 4), collapse = ", "),
                "| CP:", paste(round(summary_ebmr[, "cp"], 3), collapse = ", "), "\n\n")
          }
        }
      }
    }
  }
}

#------------------------------------------------------------------------------#
# Show Functions (similar to Simulation.r)
#------------------------------------------------------------------------------#

#' Show available regression scenarios
show_scenarios_regression <- function() {
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("Available Regression Scenarios\n")
  cat(strrep("=", 70), "\n\n")

  for (id in names(SCENARIOS_REGRESSION)) {
    s <- SCENARIOS_REGRESSION[[id]]
    cat(sprintf("%-10s %s\n", s$id, s$description))
    cat(sprintf("           PS Spec: %s | Settings: %s | Model type: %s\n",
                s$ps_spec_id, paste(s$settings, collapse = ", "), s$model_type))
    cat("\n")
  }
}

#' Show current regression configuration
show_config_regression <- function() {
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("Regression Simulation Configuration\n")
  cat(strrep("=", 70), "\n\n")
  cat("  Replicate number:", CONFIG_REGRESSION$replicate_num, "\n")
  cat("  Default sample size:", CONFIG_REGRESSION$n_default, "\n")
  cat("\n")
  cat("  True regression coefficients:\n")
  for (setting in names(THETA_TRUE)) {
    cat("    ", setting, ":", paste(round(THETA_TRUE[[setting]], 4), collapse = ", "), "\n")
  }
  cat("\n")
}

#------------------------------------------------------------------------------#
# Comprehensive Summary Table Functions
#------------------------------------------------------------------------------#

#' Generate comprehensive summary table for regression simulation
#'
#' @param setting Setting name (e.g., "setting13", "setting14")
#' @param scenario_id Scenario ID (e.g., "7-1", "8-1")
#' @param n_vector Vector of sample sizes (default: c(2000, 500))
#' @param missing_rates Vector of missing rates (default: c("miss50", "miss30"))
#' @param model_combinations Vector of model combinations (default: c("1", "2", "3", "12", "13", "23", "123"))
#' @param replicate_num Number of replicates
#' @param version Version string
#' @param coef_idx Which coefficient to display (1=Int, 2=z1, 3=z2, 4=u1, 5=u2, or "all")
#' @return Data frame with summary results (invisibly)
generate_summary_table_regression <- function(setting,
                                               scenario_id,
                                               n_vector = c(2000, 500),
                                               missing_rates = c("miss50", "miss30"),
                                               model_combinations = c("1", "2", "3", "12", "13", "23", "123"),
                                               replicate_num = 1000,
                                               version = "v1",
                                               coef_idx = "all") {

  theta.true <- get_theta_true(setting)
  family <- get_family(setting)
  param_names <- c("Int", "z1", "z2", "u1", "u2")

  # Determine which coefficients to show
  if (coef_idx == "all") {
    coef_indices <- 1:5
  } else {
    coef_indices <- as.integer(coef_idx)
  }

  # Build results data frame
  results_list <- list()

  for (miss in missing_rates) {
    for (n in n_vector) {
      for (model_str in model_combinations) {
        file <- paste0(
          "Simulation_Results/EBMR_IPW_regression_", setting, "-", miss,
          "-scenario", scenario_id, "_", model_str,
          "_n", n, "_replicate", replicate_num, "_", version, ".RDS"
        )

        if (file.exists(file)) {
          sim_result <- readRDS(file)

          # Get summaries for both methods
          summary_true <- summarize_regression_results(sim_result, setting, "true_ps")
          summary_ebmr <- summarize_regression_results(sim_result, setting, "ebmr")

          for (j in coef_indices) {
            # True PS row
            results_list[[length(results_list) + 1]] <- data.frame(
              Setting = setting,
              Miss = miss,
              n = n,
              Model = model_str,
              Method = "True",
              Coef = param_names[j],
              Bias = summary_true[j, "bias"],
              ESD = summary_true[j, "esd"],
              ASE = summary_true[j, "ese"],
              CP = summary_true[j, "cp"],
              stringsAsFactors = FALSE
            )

            # EBMR row
            results_list[[length(results_list) + 1]] <- data.frame(
              Setting = setting,
              Miss = miss,
              n = n,
              Model = model_str,
              Method = "EBMR",
              Coef = param_names[j],
              Bias = summary_ebmr[j, "bias"],
              ESD = summary_ebmr[j, "esd"],
              ASE = summary_ebmr[j, "ese"],
              CP = summary_ebmr[j, "cp"],
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }

  if (length(results_list) == 0) {
    cat("No results found for the specified parameters.\n")
    return(invisible(NULL))
  }

  results_df <- do.call(rbind, results_list)

  # Print header
  cat("\n")
  cat(strrep("=", 100), "\n")
  cat("REGRESSION SIMULATION SUMMARY TABLE\n")
  cat(strrep("=", 100), "\n")
  cat("Setting:", setting, "| Family:", family, "| Scenario:", scenario_id, "\n")
  cat("True theta:", paste(round(theta.true, 4), collapse = ", "), "\n")
  cat(strrep("=", 100), "\n\n")

  # Print table for each missing rate

  for (miss in missing_rates) {
    cat(strrep("-", 100), "\n")
    cat("Missing Rate:", miss, "\n")
    cat(strrep("-", 100), "\n\n")

    subset_df <- results_df[results_df$Miss == miss, ]

    # Print header row
    cat(sprintf("%-8s %-6s %-6s %-6s %10s %10s %10s %10s\n",
                "n", "Model", "Method", "Coef", "Bias", "ESD", "ASE", "CP"))
    cat(strrep("-", 80), "\n")

    for (n in n_vector) {
      n_subset <- subset_df[subset_df$n == n, ]

      for (model_str in model_combinations) {
        model_subset <- n_subset[n_subset$Model == model_str, ]

        if (nrow(model_subset) > 0) {
          for (i in 1:nrow(model_subset)) {
            row <- model_subset[i, ]
            cat(sprintf("%-8d %-6s %-6s %-6s %10.4f %10.4f %10.4f %10.3f\n",
                        row$n, row$Model, row$Method, row$Coef,
                        row$Bias, row$ESD, row$ASE, row$CP))
          }
        }
      }
      cat("\n")
    }
  }

  invisible(results_df)
}

#' Generate compact summary table showing only key metrics
#'
#' @param setting Setting name
#' @param scenario_id Scenario ID
#' @param n_vector Vector of sample sizes
#' @param missing_rates Vector of missing rates
#' @param model_combinations Vector of model combinations
#' @param replicate_num Number of replicates
#' @param version Version string
#' @return Data frame with summary results (invisibly)
generate_compact_table_regression <- function(setting,
                                               scenario_id,
                                               n_vector = c(2000, 500),
                                               missing_rates = c("miss50", "miss30"),
                                               model_combinations = c("1", "2", "3", "12", "13", "23", "123"),
                                               replicate_num = 1000,
                                               version = "v1") {

  theta.true <- get_theta_true(setting)
  family <- get_family(setting)
  param_names <- c("Int", "z1", "z2", "u1", "u2")

  # Print header
  cat("\n")
  cat(strrep("=", 120), "\n")
  cat("COMPACT REGRESSION SUMMARY: ", setting, " | Scenario: ", scenario_id, "\n", sep = "")
  cat("Family:", family, "| True theta:", paste(round(theta.true, 4), collapse = ", "), "\n")
  cat(strrep("=", 120), "\n\n")

  # For each missing rate
  for (miss in missing_rates) {
    cat(strrep("-", 120), "\n")
    cat("Missing Rate:", miss, "\n")
    cat(strrep("-", 120), "\n\n")

    # Header: n | Model | Method | Int (Bias/CP) | z1 (Bias/CP) | z2 (Bias/CP) | u1 (Bias/CP) | u2 (Bias/CP)
    cat(sprintf("%-6s %-6s %-6s | %12s | %12s | %12s | %12s | %12s\n",
                "n", "Model", "Method", "Int", "z1", "z2", "u1", "u2"))
    cat(sprintf("%-6s %-6s %-6s | %12s | %12s | %12s | %12s | %12s\n",
                "", "", "", "Bias/CP", "Bias/CP", "Bias/CP", "Bias/CP", "Bias/CP"))
    cat(strrep("-", 120), "\n")

    for (n in n_vector) {
      for (model_str in model_combinations) {
        file <- paste0(
          "Simulation_Results/EBMR_IPW_regression_", setting, "-", miss,
          "-scenario", scenario_id, "_", model_str,
          "_n", n, "_replicate", replicate_num, "_", version, ".RDS"
        )

        if (file.exists(file)) {
          sim_result <- readRDS(file)

          # True PS
          summary_true <- summarize_regression_results(sim_result, setting, "true_ps")
          true_str <- sapply(1:5, function(j) {
            sprintf("%5.3f/%4.2f", summary_true[j, "bias"], summary_true[j, "cp"])
          })
          cat(sprintf("%-6d %-6s %-6s | %12s | %12s | %12s | %12s | %12s\n",
                      n, model_str, "True", true_str[1], true_str[2], true_str[3], true_str[4], true_str[5]))

          # EBMR
          summary_ebmr <- summarize_regression_results(sim_result, setting, "ebmr")
          ebmr_str <- sapply(1:5, function(j) {
            sprintf("%5.3f/%4.2f", summary_ebmr[j, "bias"], summary_ebmr[j, "cp"])
          })
          cat(sprintf("%-6d %-6s %-6s | %12s | %12s | %12s | %12s | %12s\n",
                      n, model_str, "EBMR", ebmr_str[1], ebmr_str[2], ebmr_str[3], ebmr_str[4], ebmr_str[5]))
        }
      }
      cat("\n")
    }
  }

  cat(strrep("=", 120), "\n")
}

#' Generate summary table for a single coefficient across all configurations
#'
#' @param setting Setting name
#' @param scenario_id Scenario ID
#' @param coef_name Coefficient name ("Int", "z1", "z2", "u1", "u2")
#' @param n_vector Vector of sample sizes
#' @param missing_rates Vector of missing rates
#' @param model_combinations Vector of model combinations
#' @param replicate_num Number of replicates
#' @param version Version string
generate_coef_table_regression <- function(setting,
                                            scenario_id,
                                            coef_name = "z1",
                                            n_vector = c(2000, 500),
                                            missing_rates = c("miss50", "miss30"),
                                            model_combinations = c("1", "2", "3", "12", "13", "23", "123"),
                                            replicate_num = 1000,
                                            version = "v1") {

  theta.true <- get_theta_true(setting)
  family <- get_family(setting)
  param_names <- c("Int", "z1", "z2", "u1", "u2")
  coef_idx <- which(param_names == coef_name)

  if (length(coef_idx) == 0) {
    stop(paste("Unknown coefficient:", coef_name,
               "\nAvailable:", paste(param_names, collapse = ", ")))
  }

  # Print header
  cat("\n")
  cat(strrep("=", 100), "\n")
  cat("COEFFICIENT SUMMARY: ", coef_name, " (", setting, " | Scenario: ", scenario_id, ")\n", sep = "")
  cat("Family:", family, "| True value:", round(theta.true[coef_idx], 4), "\n")
  cat(strrep("=", 100), "\n\n")

  # Header
  cat(sprintf("%-8s %-8s %-8s %-8s %10s %10s %10s %10s\n",
              "Miss", "n", "Model", "Method", "Bias", "ESD", "ASE", "CP"))
  cat(strrep("-", 90), "\n")

  for (miss in missing_rates) {
    for (n in n_vector) {
      for (model_str in model_combinations) {
        file <- paste0(
          "Simulation_Results/EBMR_IPW_regression_", setting, "-", miss,
          "-scenario", scenario_id, "_", model_str,
          "_n", n, "_replicate", replicate_num, "_", version, ".RDS"
        )

        if (file.exists(file)) {
          sim_result <- readRDS(file)

          # True PS
          summary_true <- summarize_regression_results(sim_result, setting, "true_ps")
          cat(sprintf("%-8s %-8d %-8s %-8s %10.4f %10.4f %10.4f %10.3f\n",
                      miss, n, model_str, "True",
                      summary_true[coef_idx, "bias"],
                      summary_true[coef_idx, "esd"],
                      summary_true[coef_idx, "ese"],
                      summary_true[coef_idx, "cp"]))

          # EBMR
          summary_ebmr <- summarize_regression_results(sim_result, setting, "ebmr")
          cat(sprintf("%-8s %-8d %-8s %-8s %10.4f %10.4f %10.4f %10.3f\n",
                      miss, n, model_str, "EBMR",
                      summary_ebmr[coef_idx, "bias"],
                      summary_ebmr[coef_idx, "esd"],
                      summary_ebmr[coef_idx, "ese"],
                      summary_ebmr[coef_idx, "cp"]))
        }
      }
    }
    cat("\n")
  }
}

#' Generate wide-format summary table with all coefficients
#'
#' Layout: rows = n (500, 2000) x estimator (IPW, 001, ..., 111) x coef (Int, z1, z2, u1, u2)
#'         columns = miss50 (Bias, ESD, ASE, CP) | miss30 (Bias, ESD, ASE, CP)
#'
#' @param setting Setting name
#' @param scenario_id Scenario ID
#' @param replicate_num Number of replicates
#' @param version Version string
generate_wide_table_regression <- function(setting,
                                            scenario_id,
                                            replicate_num = 1000,
                                            version = "v1") {

  theta.true <- get_theta_true(setting)
  family <- get_family(setting)
  param_names <- c("Int", "z1", "z2", "u1", "u2")

  # Estimator order: IPW first, then single models, pairs, all three
  estimator_order <- c("IPW", "001", "010", "100", "011", "101", "110", "111")
  model_order <- c(NA, "1", "2", "3", "12", "13", "23", "123")  # NA for IPW (uses any file)

  n_vector <- c(500, 2000)
  missing_rates <- c("miss50", "miss30")

  # Print header
  cat("\n")
  cat(strrep("=", 115), "\n")
  cat("REGRESSION SUMMARY: ", setting, " | Scenario: ", scenario_id, "\n", sep = "")
  cat("Family: ", family, " | True theta = (", paste(round(theta.true, 2), collapse = ", "), ")\n", sep = "")
  cat(strrep("=", 115), "\n\n")

  # Column headers
  cat(sprintf("%-4s %-10s %-4s |%s|%s\n",
              "", "", "",
              sprintf("%44s", "miss50"),
              sprintf("%44s", "miss30")))
  cat(sprintf("%-4s %-10s %-4s |%10s %10s %10s %10s |%10s %10s %10s %10s\n",
              "n", "Estimator", "Coef", "Bias", "ESD", "ASE", "CP", "Bias", "ESD", "ASE", "CP"))
  cat(strrep("-", 115), "\n")

  for (n in n_vector) {
    first_n_row <- TRUE

    for (i in seq_along(estimator_order)) {
      est_name <- estimator_order[i]
      model_str <- model_order[i]

      # For display
      if (est_name == "IPW") {
        display_name <- "theta_IPW"
      } else {
        display_name <- paste0("theta_", est_name)
      }

      # Load data for both missing rates
      summary_miss50 <- NULL
      summary_miss30 <- NULL

      for (miss in missing_rates) {
        # For IPW, we can use any model file (true PS is the same)
        if (is.na(model_str)) {
          test_file <- paste0(
            "Simulation_Results/EBMR_IPW_regression_", setting, "-", miss,
            "-scenario", scenario_id, "_1",
            "_n", n, "_replicate", replicate_num, "_", version, ".RDS"
          )
        } else {
          test_file <- paste0(
            "Simulation_Results/EBMR_IPW_regression_", setting, "-", miss,
            "-scenario", scenario_id, "_", model_str,
            "_n", n, "_replicate", replicate_num, "_", version, ".RDS"
          )
        }

        if (file.exists(test_file)) {
          sim_result <- readRDS(test_file)

          if (est_name == "IPW") {
            summary_res <- summarize_regression_results(sim_result, setting, "true_ps")
          } else {
            summary_res <- summarize_regression_results(sim_result, setting, "ebmr")
          }

          if (miss == "miss50") {
            summary_miss50 <- summary_res
          } else {
            summary_miss30 <- summary_res
          }
        }
      }

      # Print rows for each coefficient
      for (j in 1:5) {
        # Format values
        fmt_val <- function(x) if (is.null(x) || is.na(x)) sprintf("%10s", "-") else sprintf("%10.4f", x)
        fmt_cp <- function(x) if (is.null(x) || is.na(x)) sprintf("%10s", "-") else sprintf("%10.3f", x)

        # Get values for miss50
        if (!is.null(summary_miss50)) {
          v50 <- c(summary_miss50[j, "bias"], summary_miss50[j, "esd"],
                   summary_miss50[j, "ese"], summary_miss50[j, "cp"])
        } else {
          v50 <- c(NA, NA, NA, NA)
        }

        # Get values for miss30
        if (!is.null(summary_miss30)) {
          v30 <- c(summary_miss30[j, "bias"], summary_miss30[j, "esd"],
                   summary_miss30[j, "ese"], summary_miss30[j, "cp"])
        } else {
          v30 <- c(NA, NA, NA, NA)
        }

        # Print row
        n_display <- ifelse(first_n_row && j == 1, as.character(n), "")
        est_display <- ifelse(j == 1, display_name, "")

        cat(sprintf("%-4s %-10s %-4s |%s %s %s %s |%s %s %s %s\n",
                    n_display, est_display, param_names[j],
                    fmt_val(v50[1]), fmt_val(v50[2]), fmt_val(v50[3]), fmt_cp(v50[4]),
                    fmt_val(v30[1]), fmt_val(v30[2]), fmt_val(v30[3]), fmt_cp(v30[4])))

        first_n_row <- FALSE
      }
    }
    cat(strrep("-", 115), "\n")
  }
}
