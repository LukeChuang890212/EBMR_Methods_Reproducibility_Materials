#------------------------------------------------------------------------------#
# EBMR Simulation Functions
#
# This file contains functions for:
#   1. Running simulations (simulate, simulate_all_*)
#   2. Summarizing results (summarize_results, summarize_scenario, summarize_all_*)
#   3. Utility functions (clean_sim_result, get_mu_true, generate_data)
#------------------------------------------------------------------------------#

# Load required packages for summarization
library(knitr)
library(kableExtra)

#------------------------------------------------------------------------------#
# Utility Functions
#------------------------------------------------------------------------------#

#' Clean simulation results by removing NAs and outliers
#'
#' This is the common logic used for cleaning simulation results.
#' It checks all first 4 rows (mu_ipw, mu_ipw.true, se_ipw, se_ipw.true)
#' for NAs and outliers using the IQR method with multiplier=30.
#'
#' @param sim_result Simulation result matrix (rows are statistics, columns are replicates)
#' @param multiplier IQR multiplier for outlier detection (default: 30)
#' @param verbose Whether to print summary messages (default: TRUE)
#' @return List with:
#'   - result: cleaned simulation result matrix
#'   - n_total: total number of replicates
#'   - n_na: number of replicates removed due to NA
#'   - n_outliers: number of replicates removed as outliers
#'   - n_successful: number of successful (clean) replicates
clean_sim_result <- function(sim_result, multiplier = 30, verbose = TRUE) {
  n_total <- ncol(sim_result)
  n_na <- 0
  n_outliers <- 0

  # Remove NAs (check all first 4 rows: mu_ipw, mu_ipw.true, se_ipw, se_ipw.true)
  na_mask <- apply(sim_result[1:4, , drop = FALSE], 2, function(x) !any(is.na(x)))
  n_na <- sum(!na_mask)
  sim_result <- sim_result[, na_mask, drop = FALSE]

  # Remove outliers (check all first 4 rows using IQR method)
  if (ncol(sim_result) > 0) {
    is_outlier <- rep(FALSE, ncol(sim_result))
    for (row_idx in 1:4) {
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

#' Get true mean for a setting by generating a large sample
#'
#' @param setting Setting name (e.g., "setting11")
#' @param n_large Sample size for computing true mean (default: 10^7)
#' @return True mean value (population mean of y)
get_mu_true <- function(setting, n_large = 10^7) {
  # Use a cache environment to store computed values
  if (!exists(".mu_true_cache", envir = .GlobalEnv)) {
    assign(".mu_true_cache", new.env(), envir = .GlobalEnv)
  }
  cache <- get(".mu_true_cache", envir = .GlobalEnv)

  # Check if already computed
  cache_key <- paste0(setting, "_n", n_large)
  if (exists(cache_key, envir = cache)) {
    return(get(cache_key, envir = cache))
  }

  # Determine data generation function name (use A1 variant for true mean)
  data_gen_func_name <- paste0(setting, ".A1")

  # Check if function exists
  if (!exists(data_gen_func_name, mode = "function")) {
    warning("Data generation function ", data_gen_func_name, " not found. Using fallback values.")
    # Fallback to hardcoded values if function doesn't exist
    mu_true <- switch(setting,
                      "setting7" = 0.7,
                      "setting8" = 0.617,
                      "setting9" = 0.7,
                      "setting10" = 0.646,
                      "setting11" = 0.7,
                      "setting12" = 0.646,
                      "Cho_RM2" = 2,
                      "Cho_RM3" = 2,
                      "Cho_RM2p" = 0.5,
                      "Cho_RM3p" = 0.5,
                      "Cho_RM2q" = 0.5,
                      "Cho_RM3q" = 0.5,
                      0.7)  # default
    return(mu_true)
  }

  # Generate large sample and compute true mean
  cat("Computing true mean for", setting, "with n =", format(n_large, scientific = FALSE), "...\n")
  data_gen_func <- get(data_gen_func_name)

  # Some functions require response.rate argument
  tryCatch({
    data <- data_gen_func(n_large)
  }, error = function(e) {
    # Try with response.rate if needed
    data <<- data_gen_func(n_large, response.rate = 0.5)
  })

  mu_true <- mean(data$y)
  cat("True mean (mu.true) for", setting, "=", round(mu_true, 6), "\n")

  # Cache the result
  assign(cache_key, mu_true, envir = cache)

  return(mu_true)
}

#' Generate simulation data for a setting
#'
#' @param setting_name Name of the data generation function (e.g., "setting11.A1")
#' @param n Sample size per replicate
#' @param replicate_num Number of replicates to generate
#' @param seed Random seed for reproducibility (default: 2345)
#' @param cores Number of cores to use (default: all available minus 2)
#' @return Invisibly returns the generated data; also saves to Simulation_Data/
generate_data <- function(setting_name, n = 2000, replicate_num = 1000,
                          seed = 2345, cores = NULL) {
  # Check if function exists
  if (!exists(setting_name, mode = "function")) {
    stop(paste("Data generation function not found:", setting_name))
  }

  setting_func <- get(setting_name)

  # Create output directory if needed
  if (!dir.exists("Simulation_Data")) {
    dir.create("Simulation_Data")
  }

  # Output file path
  output_file <- paste0("Simulation_Data/", setting_name, "_n", n,
                        "_replicate", replicate_num, ".RDS")

  # Check if file already exists
  if (file.exists(output_file)) {
    cat("File already exists:", output_file, "\n")
    cat("Loading existing data...\n")
    return(invisible(readRDS(output_file)))
  }

  cat("Generating data for:", setting_name, "\n")
  cat("  Sample size (n):", n, "\n")
  cat("  Replicates:", replicate_num, "\n")

  # Set seed for reproducibility
  set.seed(seed)

  # Setup parallel cluster
  library(parallel)
  library(foreach)
  library(doSNOW)

  if (is.null(cores)) {
    cores <- detectCores() - 2
  }
  cores <- max(1, cores)  # Ensure at least 1 core

  cl <- makeCluster(cores)
  registerDoSNOW(cl)

  # Progress bar
  pb <- txtProgressBar(max = replicate_num, style = 3)
  progress <- function(i) setTxtProgressBar(pb, i)
  opts <- list(progress = progress)

  cat("  Using", cores, "cores\n")
  start <- Sys.time()

  # Generate data in parallel
  dat <- foreach(
    i = 1:replicate_num,
    .combine = 'rbind',
    .options.snow = opts
  ) %dopar% {
    setting_func(n = n)
  }

  close(pb)
  stopCluster(cl)

  elapsed <- Sys.time() - start
  cat("  Completed in", format(elapsed), "\n")

  # Save to file
  saveRDS(dat, output_file)
  cat("  Saved to:", output_file, "\n")

  return(invisible(dat))
}

#------------------------------------------------------------------------------#
# Core Simulation Function
#------------------------------------------------------------------------------#

#' Run simulation for a single configuration
#'
#' @param all_data Full dataset for all replicates
#' @param ps_model.true Function to compute true propensity scores
#' @param alpha.true True alpha parameter values
#' @param ps_specifications List of PS model specifications
#' @param n Sample size per replicate
#' @param replicate_num Number of replicates to run
#' @param save_file Path to save results (optional)
#' @return Simulation results matrix (invisibly)
simulate <- function(all_data, ps_model.true, alpha.true, ps_specifications,
                     n, replicate_num, save_file = NULL) {
  # NOTE: EBMRalgorithm and Cho2025 are temporarily disabled
  # Only EBMRalgorithmOld is used for now
  # library(EBMRalgorithm)  # DISABLED
  library(EBMRalgorithmOld)
  library(parallel)
  library(foreach)
  library(doSNOW)

  mu.true <- mean(all_data$y)

  # Setup parallel cluster

  cores <- detectCores()
  cl <- makeCluster(cores - 6)
  registerDoSNOW(cl)

  # Export local variables to workers
  clusterExport(cl, c("ps_model.true", "alpha.true", "ps_specifications", "n", "all_data"),
                envir = environment())
  # Export global functions - Cho2025 disabled
  # clusterExport(cl, "Cho2025", envir = .GlobalEnv)

  # Progress bar

  pb <- txtProgressBar(max = replicate_num, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  # EBMRalgorithm disabled from packages
  parallel_packages <- c("EBMRalgorithmOld", "stringr", "Matrix", "numDeriv")

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

      # Weight matrix function
      W <- function(g.matrix) {
        return(diag(ncol(g.matrix)))
      }

      # Old EBMR algorithm (the only method currently enabled)
      ebmr_old <- EBMRAlgorithmOld$new("y", ps_specifications, dat, W)
      result_old <- ebmr_old$EBMR_IPW(
        h_x_names = c("u1", "u2", "z1", "z2"),
        true_ps = ps_model.true(dat, alpha.true)
      )
      estimates_old <- unlist(result_old[1:4])

      # # New EBMR algorithm - DISABLED
      # ebmr <- EBMRAlgorithm$new("y", ps_specifications, dat, W)
      # result <- ebmr$EBMR_IPW(
      #   h_x_names = ps_specifications$h_x_names.list[[1]],
      #   true_ps = ps_model.true(dat, alpha.true)
      # )
      # estimates <- unlist(result[1:4])

      # # Cho2025 estimator - DISABLED
      # Cho2025_est <- Cho2025(
      #   dat,
      #   or_x_names = ps_specifications$h_x_names.list[[1]],
      #   s = sqrt(1/3),
      #   ebmr$ps_fit.list,
      #   Aux_names = ps_specifications$h_x_names.list[[1]],
      #   inv_link = ps_specifications$inv_link
      # )

      # Return only EBMRAlgorithmOld results
      # Output structure: mu_ipw, mu_ipw.true, se_ipw, se_ipw.true, alpha.hat, se_alpha.hat, nu.hat, w.hat
      c(estimates_old,
        alpha.hat = unlist(lapply(ebmr_old$ps_fit.list, function(ps_fit) ps_fit$coefficients)),
        se_alpha.hat = unlist(lapply(ebmr_old$ps_fit.list, function(ps_fit) ps_fit$se)),
        nu.hat = result_old$nu.hat,
        w.hat = result_old$w.hat)
    }, error = function(e) {
      cat("ERROR in replicate", i, ":", conditionMessage(e), "\n")
      return(NA)  # Return NA so cbind still works
    })
  }

  close(pb)
  stopCluster(cl)
  cat("Elapsed time:", round(difftime(Sys.time(), start, units = "mins"), 2), "minutes\n")

  # Save results
  if (!is.null(save_file)) {
    saveRDS(sim_result, save_file)
    cat("Results saved to:", save_file, "\n")
  }

  # Summary is printed at the end via summarize_scenario()

  gc()
  invisible(sim_result)
}

#' Print summary of simulation results
#'
#' @param sim_result Simulation result matrix
#' @param mu.true True mean value
print_simulation_summary <- function(sim_result, mu.true) {
  cat("\n")
  cat(strrep("=", 60), "\n")
  cat("Simulation Summary\n")
  cat(strrep("=", 60), "\n")
  cat("Number of replicates:", ncol(sim_result), "\n")
  cat("True mean (mu.true):", mu.true, "\n\n")

  cat("Before removing outliers:\n")
  print(rbind(
    Mean = apply(sim_result, 1, mean, na.rm = TRUE),
    SD = apply(sim_result, 1, sd, na.rm = TRUE)
  ))

  # Remove NAs and outliers using common logic
  cleaned <- clean_sim_result(sim_result, multiplier = 30, verbose = TRUE)
  sim_result <- cleaned$result

  cat("\nAfter removing outliers:\n")
  print(rbind(
    Mean = apply(sim_result, 1, mean, na.rm = TRUE),
    SD = apply(sim_result, 1, sd, na.rm = TRUE)
  ))

  # Coverage probability
  cp <- rep(NA, 2)
  for (i in 1:2) {
    ci <- cbind(
      sim_result[i, ] - 1.96 * sim_result[2 + i, ],
      sim_result[i, ] + 1.96 * sim_result[2 + i, ]
    )
    cp[i] <- mean(na.omit(apply(ci, 1, function(v) {
      ifelse(v[1] < mu.true & v[2] > mu.true, 1, 0)
    })))
  }
  names(cp) <- c("IPW_CP", "IPW_true_CP")

  cat("\nCoverage Probability (95%):\n")
  print(cp)
  cat(strrep("=", 60), "\n\n")
}

#------------------------------------------------------------------------------#
# Orchestration Functions
#------------------------------------------------------------------------------#

#' Run simulation for all model combinations and sample sizes
simulate_all_model_combinations_and_sample_sizes <- function(
    setting, scenario, missing_rate, full_ps_specifications,
    n.vector, all_data_file, alpha.true, replicate_num, version) {

  all_data <- readRDS(all_data_file)
  J <- length(full_ps_specifications$formula.list)

  for (n in n.vector) {
    # Define true PS model based on setting
    ps_model.true <- if (setting %in% c("Cho_RM2", "Cho_RM2p", "Cho_RM2q")) {
      function(dat, alpha.true) {
        eta <- cbind(rep(1, nrow(dat)), dat$x1, dat$y) %*% alpha.true
        exp(eta) / (1 + exp(eta))
      }
    } else if (setting %in% c("Cho_RM3", "Cho_RM3p", "Cho_RM3q")) {
      function(dat, alpha.true) {
        eta <- cbind(rep(1, nrow(dat)), dat$x2, dat$y) %*% alpha.true
        exp(eta) / (1 + exp(eta))
      }
    } else {
      function(dat, alpha.true) {
        1 / (1 + exp(cbind(rep(1, nrow(dat)), dat$y, dat$u1, dat$u2) %*% alpha.true))
      }
    }

    # Iterate over model combinations
    for (model_num in 1:J) {
      model_combinations <- combn(J, model_num)
      for (i in 1:ncol(model_combinations)) {
        model_set <- model_combinations[, i]
        save_file <- paste0(
          "Simulation_Results/EBMR_IPW_", setting, "-", missing_rate,
          "-scenario", scenario, "_", paste0(model_set, collapse = ""),
          "_n", n, "_replicate", replicate_num, "_", version, ".RDS"
        )
        cat("\nModel set:", paste0(model_set, collapse = ""), "| File:", save_file, "\n")

        ps_specifications <- list(
          formula.list = full_ps_specifications$formula.list[model_set],
          h_x_names.list = full_ps_specifications$h_x_names.list[model_set],
          inv_link = full_ps_specifications$inv_link
        )

        simulate(all_data, ps_model.true, alpha.true, ps_specifications,
                 n, replicate_num, save_file)
      }
    }
  }
}

#' Run simulation for all settings and missing rates
simulate_all_settings_with_all_missing_rates <- function(
    settings, missing_rates, scenario, full_ps_specifications,
    n.vector.list, all_data_file.list, alpha.true.list, version) {

  for (setting in settings) {
    for (missing_rate in missing_rates) {
      for (i in 1:length(all_data_file.list[[setting]][[missing_rate]])) {
        simulate_all_model_combinations_and_sample_sizes(
          setting = setting,
          scenario = scenario,
          missing_rate = missing_rate,
          full_ps_specifications = full_ps_specifications,
          n.vector = n.vector.list[[i]],
          all_data_file = all_data_file.list[[setting]][[missing_rate]][[i]],
          alpha.true = alpha.true.list[[setting]][[missing_rate]][[i]],
          replicate_num = replicate_num,
          version = version
        )
      }
    }
  }
}

#------------------------------------------------------------------------------#
# Result Summarization Functions
#------------------------------------------------------------------------------#

#' Summarize results for a single simulation
#'
#' @param sim_result Simulation result matrix
#' @param pe_index Row index for point estimate
#' @param ese_index Row index for estimated standard error
#' @param mu.true True mean value
#' @param is.original If FALSE, remove outliers using clean_sim_result()
#' @return Formatted vector of (Bias, ESD, ESE, CP)
summarize_results <- function(sim_result, pe_index, ese_index, mu.true, is.original) {
  if (!is.original) {
    # Use common cleaning logic
    cleaned <- clean_sim_result(sim_result, multiplier = 30, verbose = TRUE)
    sim_result <- cleaned$result
  } else {
    cat("Replicates:", ncol(sim_result), "\n")
  }

  # Helper function: round to 3 decimal places, if last digit is 0 change to 1
  round_3dp <- function(x) {
    result <- round(x, 3)
    # If the 3rd decimal digit is 0, change to 1
    result <- ifelse(round(result * 1000) %% 10 == 0, result + 0.001, result)
    return(result)
  }

  pe <- mean(sim_result[pe_index, ], na.rm = TRUE)
  esd <- sd(sim_result[pe_index, ], na.rm = TRUE)
  bias <- round_3dp(pe - mu.true)
  esd <- round_3dp(esd)

  ese <- sim_result[ese_index, ]
  ci <- cbind(
    sim_result[pe_index, ] - 1.96 * ese,
    sim_result[pe_index, ] + 1.96 * ese
  )
  coverage <- apply(ci, 1, function(interval) {
    ifelse(mu.true >= interval[1] & mu.true <= interval[2], 1, 0)
  })
  cp <- round_3dp(mean(na.omit(coverage)))
  ese <- round_3dp(mean(ese, na.rm = TRUE))

  return(format(c(bias, esd, ese, cp), nsmall = 3))
}

#' Summarize a single scenario's results (simplified interface)
#'
#' @param scenario Scenario ID (e.g., "7-2")
#' @param setting Setting name (e.g., "setting11")
#' @param missing_rate Missing rate (e.g., "miss50")
#' @param n Sample size
#' @param replicate_num Number of replicates
#' @param version Version string
#' @param model_set Model set (default: "1")
#' @param mu.true True mean (if NULL, auto-determined from setting)
#' @return Summary data frame
summarize_scenario <- function(scenario, setting, missing_rate, n, replicate_num,
                               version, model_set = "1", mu.true = NULL) {
  # Build filename
  result_file <- paste0(
    "Simulation_Results/EBMR_IPW_", setting, "-", missing_rate,
    "-scenario", scenario, "_", model_set,
    "_n", n, "_replicate", replicate_num, "_", version, ".RDS"
  )

  if (!file.exists(result_file)) {
    cat("File not found:", result_file, "\n")
    return(NULL)
  }

  sim_result <- readRDS(result_file)

  # Get true mean
  if (is.null(mu.true)) {
    mu.true <- get_mu_true(setting)
  }

  cat("\n")
  cat(strrep("=", 60), "\n")
  cat("Summary: Scenario", scenario, "| Setting:", setting, "\n")
  cat(strrep("=", 60), "\n")
  cat("File:", basename(result_file), "\n")
  cat("Replicates:", ncol(sim_result), "\n")
  cat("Sample size:", n, "\n")
  cat("True mean:", mu.true, "\n")
  cat(strrep("-", 60), "\n\n")

  # Remove NA columns
  sim_result <- sim_result[, !is.na(sim_result[3, ]), drop = FALSE]

  if (ncol(sim_result) == 0) {
    cat("No valid results.\n")
    return(NULL)
  }

  # IPW estimates (rows 1-4: ipw, ipw.true, se_ipw, se_ipw.true)
  results <- data.frame(
    Estimator = c("IPW", "IPW (true PS)"),
    Mean = c(mean(sim_result[1, ], na.rm = TRUE),
             mean(sim_result[2, ], na.rm = TRUE)),
    Bias = c(mean(sim_result[1, ], na.rm = TRUE) - mu.true,
             mean(sim_result[2, ], na.rm = TRUE) - mu.true),
    ESD = c(sd(sim_result[1, ], na.rm = TRUE),
            sd(sim_result[2, ], na.rm = TRUE)),
    ESE = c(mean(sim_result[3, ], na.rm = TRUE),
            mean(sim_result[4, ], na.rm = TRUE))
  )

  # Coverage probability
  for (i in 1:2) {
    ci_lower <- sim_result[i, ] - 1.96 * sim_result[i + 2, ]
    ci_upper <- sim_result[i, ] + 1.96 * sim_result[i + 2, ]
    results$CP[i] <- mean((ci_lower <= mu.true) & (mu.true <= ci_upper), na.rm = TRUE)
  }

  # Helper function: round to 3 decimal places, if last digit is 0 change to 1
  round_3dp <- function(x) {
    result <- round(x, 3)
    result <- ifelse(round(result * 1000) %% 10 == 0, result + 0.001, result)
    return(result)
  }

  # Format and print
  results[, 2:6] <- sapply(results[, 2:6], round_3dp)
  print(results, row.names = FALSE)
  cat(strrep("=", 60), "\n\n")

  invisible(results)
}

#' Summarize all model combinations and sample sizes
#' NOTE: Updated for EBMRalgorithmOld-only output format
#' Output structure: mu_ipw(1), mu_ipw.true(2), se_ipw(3), se_ipw.true(4), alpha.hat, se_alpha.hat, nu.hat, w.hat
summarize_all_model_combinations_and_sample_sizes <- function(
    setting, scenario, J, missing_rate, n.vector,
    replicate_num, mu.true, version, is.original) {

  # With Cho2025 disabled, we have 8 rows per n:
  # IPW true PS (1), 7 model combinations (1,2,3,12,13,23,123)
  result <- matrix(NA, 8 * length(n.vector), 4)
  j <- 1

  for (n in n.vector) {
    cat("n =", n, "\n")
    for (model_num in 1:J) {
      model_combinations <- combn(J, model_num)
      for (i in 1:ncol(model_combinations)) {
        model_set <- model_combinations[, i]
        read_file <- paste0(
          "Simulation_Results/EBMR_IPW_", setting, "-", missing_rate,
          "-scenario", scenario, "_", paste0(model_set, collapse = ""),
          "_n", n, "_replicate", replicate_num, "_", version, ".RDS"
        )
        sim_result <- readRDS(read_file)

        # For the first single model (model 1), also report IPW with true PS
        if (model_num == 1 & i == 1) {
          # IPW with true PS: row 2 (mu_ipw.true), row 4 (se_ipw.true)
          result[j, ] <- summarize_results(sim_result, pe_index = 2, ese_index = 4,
                                           mu.true, is.original)
          j <- j + 1
        }

        # Ensemble IPW estimate: row 1 (mu_ipw), row 3 (se_ipw)
        result[j, ] <- summarize_results(
          sim_result,
          pe_index = 1,
          ese_index = 3,
          mu.true, is.original
        )
        j <- j + 1

        # For full model (all 3), also report ensemble result
        if (model_num == 3) {
          # This was for ensemble - already reported above, skip duplicate
          # Previously: result[j, ] <- summarize_results(sim_result, pe_index = 1, ese_index = 3, ...)

          # Cho2025 estimator is DISABLED - skip this row
          # sim_result <- sim_result[, which.not.extreme(sim_result[nrow(sim_result), ])]
          # result[j, ] <- c(
          #   format(round(mean(sim_result[nrow(sim_result), ]) - mu.true, 3), nsmall = 3),
          #   format(round(sd(sim_result[nrow(sim_result), ]), 3), nsmall = 3),
          #   "-", "-"
          # )
          # j <- j + 1
        }
      }
    }
  }
  return(result)
}

#' Print console summary table for a setting
#' NOTE: Updated for EBMRalgorithmOld-only output format
#' Output structure: mu_ipw(1), mu_ipw.true(2), se_ipw(3), se_ipw.true(4), ...
print_console_summary <- function(setting, scenario, missing_rates, n.vector,
                                   replicate_num, mu.true, version) {
  # Helper function: round to 3 decimal places, if last digit is 0 change to 1
  round_3dp <- function(x) {
    result <- round(x, 3)
    # If the 3rd decimal digit is 0, change to 1
    result <- ifelse(round(result * 1000) %% 10 == 0, result + 0.001, result)
    return(result)
  }

  cat("\n")
  cat(strrep("=", 90), "\n")
  cat("SIMULATION RESULTS: Scenario", scenario, "| Setting", setting, "| mu.true =", mu.true, "\n")
  cat(strrep("=", 90), "\n\n")

  for (n in n.vector) {
    cat("Sample size n =", n, "\n")
    cat(strrep("-", 90), "\n")
    cat(sprintf("%-20s | %8s %8s %8s %8s | %8s %8s %8s %8s\n",
                "Estimator", "Bias", "ESD", "ESE", "CP", "Bias", "ESD", "ESE", "CP"))
    cat(sprintf("%-20s | %35s | %35s\n", "",
                paste0(missing_rates[1], " missing"),
                paste0(missing_rates[2], " missing")))
    cat(strrep("-", 90), "\n")

    # IPW with true PS (from model "1" file, rows 2 and 4)
    row_ipw <- c()
    for (miss in missing_rates) {
      file <- paste0("Simulation_Results/EBMR_IPW_", setting, "-", miss,
                     "-scenario", scenario, "_1_n", n, "_replicate", replicate_num, "_", version, ".RDS")
      if (file.exists(file)) {
        sim_result <- readRDS(file)
        bias <- round_3dp(mean(sim_result[2, ], na.rm = TRUE) - mu.true)
        esd <- round_3dp(sd(sim_result[2, ], na.rm = TRUE))
        ese <- round_3dp(mean(sim_result[4, ], na.rm = TRUE))
        ci_lower <- sim_result[2, ] - 1.96 * sim_result[4, ]
        ci_upper <- sim_result[2, ] + 1.96 * sim_result[4, ]
        cp <- round_3dp(mean((ci_lower <= mu.true) & (mu.true <= ci_upper), na.rm = TRUE))
        row_ipw <- c(row_ipw, bias, esd, ese, cp)
      } else {
        row_ipw <- c(row_ipw, NA, NA, NA, NA)
      }
    }
    cat(sprintf("%-20s | %8.3f %8.3f %8.3f %8.3f | %8.3f %8.3f %8.3f %8.3f\n",
                "mu_IPW (true PS)", row_ipw[1], row_ipw[2], row_ipw[3], row_ipw[4],
                row_ipw[5], row_ipw[6], row_ipw[7], row_ipw[8]))
    cat(strrep("-", 90), "\n")

    # Model combinations - with EBMRalgorithmOld only:
    # mu_ipw is row 1, se_ipw is row 3
    for (model_set in c("1", "2", "3", "12", "13", "23", "123")) {
      row <- c()
      for (miss in missing_rates) {
        file <- paste0("Simulation_Results/EBMR_IPW_", setting, "-", miss,
                       "-scenario", scenario, "_", model_set, "_n", n,
                       "_replicate", replicate_num, "_", version, ".RDS")
        if (file.exists(file)) {
          sim_result <- readRDS(file)
          # Fixed indices for EBMRalgorithmOld-only output
          pe_idx <- 1   # mu_ipw
          ese_idx <- 3  # se_ipw
          bias <- round_3dp(mean(sim_result[pe_idx, ], na.rm = TRUE) - mu.true)
          esd <- round_3dp(sd(sim_result[pe_idx, ], na.rm = TRUE))
          ese <- round_3dp(mean(sim_result[ese_idx, ], na.rm = TRUE))
          ci_lower <- sim_result[pe_idx, ] - 1.96 * sim_result[ese_idx, ]
          ci_upper <- sim_result[pe_idx, ] + 1.96 * sim_result[ese_idx, ]
          cp <- round_3dp(mean((ci_lower <= mu.true) & (mu.true <= ci_upper), na.rm = TRUE))
          row <- c(row, bias, esd, ese, cp)
        } else {
          row <- c(row, NA, NA, NA, NA)
        }
      }
      cat(sprintf("%-20s | %8.3f %8.3f %8.3f %8.3f | %8.3f %8.3f %8.3f %8.3f\n",
                  paste0("mu_", model_set), row[1], row[2], row[3], row[4],
                  row[5], row[6], row[7], row[8]))
    }
    cat("\n")
  }
}

#' Summarize all settings with all missing rates (generates LaTeX tables)
#' NOTE: Updated for EBMRalgorithmOld-only output (Cho2025 disabled)
summarize_all_settings_with_all_missing_rates <- function(
    settings, missing_rates, scenario, J, n.vector,
    all_data_file.list, alpha_true.list, version, is.original = TRUE) {

  summary_tbls <- list()

  for (j in 1:length(settings)) {
    setting <- settings[j]
    cat("Processing:", setting, "\n")

    # 9 rows per n (Cho2025 disabled): IPW true PS + 7 model combinations + ensemble
    # Actually: IPW_true, mu_1, mu_2, mu_3, mu_12, mu_13, mu_23, mu_123 = 8 rows per n
    results_with_all_missing_rates <- matrix(NA, 8 * length(n.vector), 4 * length(missing_rates))
    mu.true <- get_mu_true(setting)

    # Print console summary table
    print_console_summary(setting, scenario, missing_rates, n.vector,
                          replicate_num, mu.true, version)

    for (i in 1:length(missing_rates)) {
      missing_rate <- missing_rates[i]
      results_with_all_missing_rates[, ((i - 1) * 4 + 1):((i - 1) * 4 + 4)] <-
        summarize_all_model_combinations_and_sample_sizes(
          setting = setting,
          scenario = scenario,
          J = J,
          missing_rate = missing_rate,
          n.vector = n.vector,
          replicate_num = replicate_num,
          mu.true = mu.true,
          is.original = is.original,
          version = version
        )
    }

    # Estimator names for LaTeX (Cho2025/MCEL removed)
    estimator_names <- if (substr(scenario, 3, 3) == "1" || scenario %in% c("cho1")) {
      rep(c("$\\hat{\\mu}_{\\text{IPW}}$",
            "$\\hat{\\mu}_{100}$", "$\\hat{\\mu}_{010}$", "$\\hat{\\mu}_{001}$",
            "$\\hat{\\mu}_{110}$", "$\\hat{\\mu}_{101}$", "$\\hat{\\mu}_{011}$",
            "$\\hat{\\mu}_{111}$"),
          length(n.vector))
    } else {
      rep(c("$\\tilde{\\mu}_{\\text{IPW}}$",
            "$\\tilde{\\mu}_{100}$", "$\\tilde{\\mu}_{010}$", "$\\tilde{\\mu}_{001}$",
            "$\\tilde{\\mu}_{110}$", "$\\tilde{\\mu}_{101}$", "$\\tilde{\\mu}_{011}$",
            "$\\tilde{\\mu}_{111}$"),
          length(n.vector))
    }

    # Check if there are any negative numbers in the results
    has_negative <- any(as.numeric(results_with_all_missing_rates) < 0, na.rm = TRUE)

    # If there are negative numbers, add "~" prefix to positive numbers for alignment
    if (has_negative) {
      results_with_all_missing_rates <- apply(results_with_all_missing_rates, c(1, 2), function(x) {
        if (is.na(x)) return(NA)
        num_val <- as.numeric(x)
        if (!is.na(num_val) && num_val >= 0) {
          # Remove leading/trailing whitespace before adding "~"
          return(paste0("~", trimws(x)))
        }
        # Also trim whitespace for negative numbers
        return(trimws(x))
      })
    }

    results_with_all_missing_rates <- cbind(estimator_names,
                                            as.data.frame(results_with_all_missing_rates))
    colnames(results_with_all_missing_rates) <- c("",
                                                  rep(c("Bias", "ESD", "ESE", "CP"), length(missing_rates)))

    summary_tbls[[j]] <- kable(
      results_with_all_missing_rates,
      format = "latex", align = "c", booktabs = TRUE,
      escape = FALSE, linesep = "",
      caption = paste0(
        "Comparison between different estimators under the Scenario ", scenario,
        " of Setting ", substr(setting, 9, 9),
        " with $\\mu_0$ approximately ", round(mu.true, 3), ". ",
        "The $\\bm{\\alpha}_0$ in $\\pi(\\bm{U}, Y; \\bm{\\alpha}_0)$ ",
        "that leads to $50\\%$ of missingness in $Y$ is $(",
        paste(alpha_true.list[[setting]][[1]][[1]], collapse = ", "),
        ")^{\\top}$ and that leads to $30\\%$ of missingness is $(",
        paste(alpha_true.list[[setting]][[2]][[1]], collapse = ", "), ")^{\\top}$."
      )
    ) %>%
      kable_styling(full_width = FALSE, latex_options = c("hold_position", "scale_down")) %>%
      add_header_above(c("", "$50\\%$ missing" = 4, "$30\\%$ missing" = 4))
  }

  return(summary_tbls)
}
