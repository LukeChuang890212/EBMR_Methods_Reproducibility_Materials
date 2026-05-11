#------------------------------------------------------------------------------#
# LaTeX Table Generation Script
#
# This script generates LaTeX tables from simulation results, matching the
# format in simulation_tables.rtf. Tables are generated for each scenario
# with captions that include true alpha values and population mean of y.
#------------------------------------------------------------------------------#

# Load required packages
library(dplyr)
library(knitr)
library(kableExtra)

# Load the simulation framework components
source("Basic_setup.r")
source("Data_Generation.r")
source("Simulation.r")
source("config/scenarios.R")

#------------------------------------------------------------------------------#
# Table Configuration
#------------------------------------------------------------------------------#

# Define table configurations: scenario -> setting mapping
# Each entry specifies which setting and scenario combination produces a table
TABLE_CONFIG <- list(
  # Main paper tables: setting3 (continuous Y) and setting4 (binary Y)
  # Scenarios: 7-1, 8-1, 8-2, 8-3, 8-4, 9-1

  # Scenario 7-1: Correct model with z1 interactions (Model 3: u2+z1)
  list(setting = "setting3", scenario = "7-1",
       description = "correct model with z1 interactions",
       mu_desc = "continuous outcome $Y$, where the set of candidate models includes a correctly specified model"),

  list(setting = "setting4", scenario = "7-1",
       description = "correct model with z1 interactions",
       mu_desc = "binary outcome $Y$, where the set of candidate models includes a correctly specified model"),

  # Scenario 7-2: Locally misspecified model with z1 interactions
  list(setting = "setting3", scenario = "7-2",
       description = "locally misspecified model with z1 interactions",
       mu_desc = "continuous outcome $Y$, where the set of candidate models includes a locally misspecified model"),

  list(setting = "setting4", scenario = "7-2",
       description = "locally misspecified model with z1 interactions",
       mu_desc = "binary outcome $Y$, where the set of candidate models includes a locally misspecified model"),

  # Scenario 7-3: Locally misspecified model with z1 interactions (h_x alt1)
  list(setting = "setting3", scenario = "7-3",
       description = "locally misspecified model with z1 interactions (h_x alt1)",
       mu_desc = "continuous outcome $Y$, where the set of candidate models includes a locally misspecified model"),

  list(setting = "setting4", scenario = "7-3",
       description = "locally misspecified model with z1 interactions (h_x alt1)",
       mu_desc = "binary outcome $Y$, where the set of candidate models includes a locally misspecified model"),

  # Scenario 7-4: Locally misspecified model with z1 interactions (h_x alt2)
  list(setting = "setting3", scenario = "7-4",
       description = "locally misspecified model with z1 interactions (h_x alt2)",
       mu_desc = "continuous outcome $Y$, where the set of candidate models includes a locally misspecified model"),

  list(setting = "setting4", scenario = "7-4",
       description = "locally misspecified model with z1 interactions (h_x alt2)",
       mu_desc = "binary outcome $Y$, where the set of candidate models includes a locally misspecified model"),

  # Scenario 8-1: Correct model with z1/z2 interactions (Model 2: u1+z1, Model 3: u2+z2)
  list(setting = "setting3", scenario = "8-1",
       description = "correct model with z1/z2 interactions",
       mu_desc = "continuous outcome $Y$, where the set of candidate models includes a correctly specified model"),

  list(setting = "setting4", scenario = "8-1",
       description = "correct model with z1/z2 interactions",
       mu_desc = "binary outcome $Y$, where the set of candidate models includes a correctly specified model"),

  # Scenario 8-2: Locally misspecified model with z1/z2 interactions
  list(setting = "setting3", scenario = "8-2",
       description = "locally misspecified model with z1/z2 interactions",
       mu_desc = "continuous outcome $Y$, where the set of candidate models includes a locally misspecified model"),

  list(setting = "setting4", scenario = "8-2",
       description = "locally misspecified model with z1/z2 interactions",
       mu_desc = "binary outcome $Y$, where the set of candidate models includes a locally misspecified model"),

  # Scenario 8-3: Locally misspecified model with z1/z2 interactions (h_x alt1)
  list(setting = "setting3", scenario = "8-3",
       description = "locally misspecified model with z1/z2 interactions (h_x alt1)",
       mu_desc = "continuous outcome $Y$, where the set of candidate models includes a locally misspecified model"),

  list(setting = "setting4", scenario = "8-3",
       description = "locally misspecified model with z1/z2 interactions (h_x alt1)",
       mu_desc = "binary outcome $Y$, where the set of candidate models includes a locally misspecified model"),

  # Scenario 8-4: Locally misspecified model with z1/z2 interactions (h_x alt2)
  list(setting = "setting3", scenario = "8-4",
       description = "locally misspecified model with z1/z2 interactions (h_x alt2)",
       mu_desc = "continuous outcome $Y$, where the set of candidate models includes a locally misspecified model"),

  list(setting = "setting4", scenario = "8-4",
       description = "locally misspecified model with z1/z2 interactions (h_x alt2)",
       mu_desc = "binary outcome $Y$, where the set of candidate models includes a locally misspecified model"),

  # Scenario 9-1: Correct model with z2 interactions (Model 2: u1+z2, Model 3: u2+z2)
  list(setting = "setting3", scenario = "9-1",
       description = "correct model with z2 interactions",
       mu_desc = "continuous outcome $Y$, where the set of candidate models includes a correctly specified model"),

  list(setting = "setting4", scenario = "9-1",
       description = "correct model with z2 interactions",
       mu_desc = "binary outcome $Y$, where the set of candidate models includes a correctly specified model"),

  # Scenario 9-2: Locally misspecified model with z2 interactions
  list(setting = "setting3", scenario = "9-2",
       description = "locally misspecified model with z2 interactions",
       mu_desc = "continuous outcome $Y$, where the set of candidate models includes a locally misspecified model"),

  list(setting = "setting4", scenario = "9-2",
       description = "locally misspecified model with z2 interactions",
       mu_desc = "binary outcome $Y$, where the set of candidate models includes a locally misspecified model"),

  # Scenario 9-3: Locally misspecified model with z2 interactions (h_x alt1)
  list(setting = "setting3", scenario = "9-3",
       description = "locally misspecified model with z2 interactions (h_x alt1)",
       mu_desc = "continuous outcome $Y$, where the set of candidate models includes a locally misspecified model"),

  list(setting = "setting4", scenario = "9-3",
       description = "locally misspecified model with z2 interactions (h_x alt1)",
       mu_desc = "binary outcome $Y$, where the set of candidate models includes a locally misspecified model"),

  # Scenario 9-4: Locally misspecified model with z2 interactions (h_x alt2)
  list(setting = "setting3", scenario = "9-4",
       description = "locally misspecified model with z2 interactions (h_x alt2)",
       mu_desc = "continuous outcome $Y$, where the set of candidate models includes a locally misspecified model"),

  list(setting = "setting4", scenario = "9-4",
       description = "locally misspecified model with z2 interactions (h_x alt2)",
       mu_desc = "binary outcome $Y$, where the set of candidate models includes a locally misspecified model"),

  # Setting 15: Binary Y with reduced m(x) coefficients
  list(setting = "setting15", scenario = "7-1",
       description = "correct model with z1 interactions (reduced m(x))",
       mu_desc = "binary outcome $Y$ with reduced $m(\\bm{x})$ coefficients, where the set of candidate models includes a correctly specified model"),

  list(setting = "setting15", scenario = "8-1",
       description = "correct model with z1/z2 interactions (reduced m(x))",
       mu_desc = "binary outcome $Y$ with reduced $m(\\bm{x})$ coefficients, where the set of candidate models includes a correctly specified model"),

  list(setting = "setting15", scenario = "8-2",
       description = "locally misspecified model with z1/z2 interactions (reduced m(x))",
       mu_desc = "binary outcome $Y$ with reduced $m(\\bm{x})$ coefficients, where the set of candidate models includes a locally misspecified model"),

  list(setting = "setting15", scenario = "8-3",
       description = "locally misspecified model with z1/z2 interactions (h_x alt1, reduced m(x))",
       mu_desc = "binary outcome $Y$ with reduced $m(\\bm{x})$ coefficients, where the set of candidate models includes a locally misspecified model"),

  list(setting = "setting15", scenario = "8-4",
       description = "locally misspecified model with z1/z2 interactions (h_x alt2, reduced m(x))",
       mu_desc = "binary outcome $Y$ with reduced $m(\\bm{x})$ coefficients, where the set of candidate models includes a locally misspecified model"),

  list(setting = "setting15", scenario = "9-1",
       description = "correct model with z2 interactions (reduced m(x))",
       mu_desc = "binary outcome $Y$ with reduced $m(\\bm{x})$ coefficients, where the set of candidate models includes a correctly specified model"),

  # Setting 16: Continuous Y with tilt = y+u1+u2 (no z1, z2 in tilt)
  list(setting = "setting16", scenario = "8-2",
       description = "locally misspecified model (tilt=y+u1+u2)",
       mu_desc = "continuous outcome $Y$, where the set of candidate models includes a locally misspecified model with tilt function $\\exp(n^{-1/2}(Y+U_1+U_2))$"),

  # Setting 17: Binary Y with tilt = y+u1+u2 (no z1, z2 in tilt)
  list(setting = "setting17", scenario = "8-2",
       description = "locally misspecified model (tilt=y+u1+u2, binary Y)",
       mu_desc = "binary outcome $Y$, where the set of candidate models includes a locally misspecified model with tilt function $\\exp(n^{-1/2}(Y+U_1+U_2))$"),

  # Cho et al. scenarios
  list(setting = "Cho_M1_gamma000", scenario = "cho1",
       description = "Cho M1 gamma=0",
       mu_desc = "$\\gamma = 0$ and missing mechanism (M1), where $\\pi_1(\\bm{X}, Y; \\bm{\\alpha}_1)$ is the correctly specified model"),

  list(setting = "Cho_M2_gamma000", scenario = "cho2",
       description = "Cho M2 gamma=0",
       mu_desc = "$\\gamma = 0$ and missing mechanism (M2), where $\\pi_2(\\bm{X}, Y; \\bm{\\alpha}_2)$ is the correctly specified model"),

  list(setting = "Cho_M1_gamma005", scenario = "cho1",
       description = "Cho M1 gamma=0.05",
       mu_desc = "$\\gamma = 0.05$ and missing mechanism (M1), where $\\pi_1(\\bm{X}, Y; \\bm{\\alpha}_1)$ is the correctly specified model"),

  list(setting = "Cho_M2_gamma005", scenario = "cho2",
       description = "Cho M2 gamma=0.05",
       mu_desc = "$\\gamma = 0.05$ and missing mechanism (M2), where $\\pi_2(\\bm{X}, Y; \\bm{\\alpha}_2)$ is the correctly specified model")
)

#------------------------------------------------------------------------------#
# Helper Functions
#------------------------------------------------------------------------------#

#' Get alpha string for caption
#'
#' @param setting Setting name
#' @param missing_rate "miss50" or "miss30"
#' @param scenario Scenario ID (e.g., "7-1", "8-2") - used to determine correct vs misspecified
#' @return Formatted alpha vector string
get_alpha_string <- function(setting, missing_rate, scenario = NULL) {
  # For locally misspecified scenarios, use misspecified alpha list
  is_misspecified <- !is.null(scenario) && scenario %in% c("7-2", "7-3", "7-4", "8-2", "8-3", "8-4", "9-2", "9-3", "9-4")

  if (is_misspecified) {
    alpha_list <- misspecified_model_alpha.true.list[[setting]][[missing_rate]]
  } else {
    alpha_list <- correct_model_alpha.true.list[[setting]][[missing_rate]]
    if (is.null(alpha_list)) {
      alpha_list <- misspecified_model_alpha.true.list[[setting]][[missing_rate]]
    }
  }
  if (is.null(alpha_list)) return(NULL)

  alpha <- alpha_list[[1]]
  # Format each alpha value with consistent decimal places
  alpha_formatted <- sapply(alpha, function(x) sprintf("%.2f", x))
  paste0("(", paste(alpha_formatted, collapse = ", "), ")^{\\top}")
}

#' Format number with proper spacing
#'
#' @param x Numeric value
#' @return Formatted string with ~ for positive numbers
format_num <- function(x) {
  if (is.na(x)) return("-")
  rounded <- round(x, 3)
  # If value rounds to exactly 0, show as 0.001 instead
  if (rounded == 0) rounded <- 0.001
  formatted <- sprintf("%.3f", rounded)
  # Ensure 3 decimal places, add ~ prefix for non-negative numbers for alignment
  # Check rounded value (not original x) to handle small negatives that round to 0
  if (rounded >= 0) {
    formatted <- paste0("~", formatted)
  }
  formatted
}

#' Generate table for a specific setting and scenario
#'
#' @param setting Setting name (e.g., "setting13")
#' @param scenario Scenario ID (e.g., "7-1")
#' @param version Version string
#' @param n.vector Vector of sample sizes
#' @param missing_rates Vector of missing rates
#' @param replicate_num Number of replicates
#' @return LaTeX table string
generate_table <- function(setting, scenario, version = "test15",
                           n.vector = c(500, 2000),
                           missing_rates = c("miss50", "miss30"),
                           replicate_num = 1000,
                           mu.true = NULL) {

  # Get true mean: use provided value or compute from large sample
  if (is.null(mu.true)) {
    mu.true <- get_mu_true(setting)
  }

  # Determine if this is a correct model scenario (use hat) or misspecified (use tilde)
  use_hat <- (substr(scenario, 3, 3) == "1" || scenario %in% c("cho1", "cho2"))
  mu_symbol <- if (use_hat) "\\hat{\\mu}" else "\\tilde{\\mu}"

  # Estimator row names
  estimator_names <- c(
    paste0("$", mu_symbol, "_{\\text{IPW}}$"),
    paste0("$", mu_symbol, "_{100}$"),
    paste0("$", mu_symbol, "_{010}$"),
    paste0("$", mu_symbol, "_{001}$"),
    paste0("$", mu_symbol, "_{110}$"),
    paste0("$", mu_symbol, "_{101}$"),
    paste0("$", mu_symbol, "_{011}$"),
    paste0("$", mu_symbol, "_{111}$")
  )

  # Model sets in order: IPW (from model 1), then 1,2,3,12,13,23,123
  model_sets <- c("1", "2", "3", "12", "13", "23", "123")

  # Results matrix: 8 rows per n (IPW + 7 model combinations)
  n_rows_per_n <- 8
  results <- matrix(NA, n_rows_per_n * length(n.vector), 4 * length(missing_rates))

  for (n_idx in seq_along(n.vector)) {
    n <- n.vector[n_idx]
    cat("  Processing n =", n, "\n")

    for (miss_idx in seq_along(missing_rates)) {
      missing_rate <- missing_rates[miss_idx]
      col_start <- (miss_idx - 1) * 4 + 1

      # Row 1: IPW with true PS (from model "1" file)
      row_base <- (n_idx - 1) * n_rows_per_n

      read_file <- paste0(
        "Simulation_Results/EBMR_IPW_", setting, "-", missing_rate,
        "-scenario", scenario, "_1",
        "_n", n, "_replicate", replicate_num, "_", version, ".RDS"
      )

      if (file.exists(read_file)) {
        sim_result <- readRDS(read_file)
        cleaned <- clean_sim_result(sim_result, multiplier = 3, verbose = FALSE)
        sim_result <- cleaned$result

        # IPW with true PS: row 2 (mu_ipw.true), row 4 (se_ipw.true)
        pe <- mean(sim_result[2, ], na.rm = TRUE)
        esd <- sd(sim_result[2, ], na.rm = TRUE)
        ese <- mean(sim_result[4, ], na.rm = TRUE)
        ci <- cbind(
          sim_result[2, ] - 1.96 * sim_result[4, ],
          sim_result[2, ] + 1.96 * sim_result[4, ]
        )
        cp <- mean(apply(ci, 1, function(v) mu.true >= v[1] & mu.true <= v[2]), na.rm = TRUE)

        results[row_base + 1, col_start:(col_start + 3)] <- c(
          format_num(pe - mu.true), format_num(esd), format_num(ese), format_num(cp)
        )
      } else {
        cat("    WARNING: File not found:", basename(read_file), "\n")
      }

      # Rows 2-8: Model combinations (1, 2, 3, 12, 13, 23, 123)
      for (model_idx in seq_along(model_sets)) {
        model_set <- model_sets[model_idx]
        row_num <- row_base + 1 + model_idx  # +1 for IPW row, +model_idx for this model

        read_file <- paste0(
          "Simulation_Results/EBMR_IPW_", setting, "-", missing_rate,
          "-scenario", scenario, "_", model_set,
          "_n", n, "_replicate", replicate_num, "_", version, ".RDS"
        )

        if (!file.exists(read_file)) {
          cat("    WARNING: File not found:", basename(read_file), "\n")
          next
        }

        sim_result <- readRDS(read_file)
        cleaned <- clean_sim_result(sim_result, multiplier = 3, verbose = FALSE)
        sim_result <- cleaned$result

        # Ensemble IPW: row 1 (mu_ipw), row 3 (se_ipw)
        pe <- mean(sim_result[1, ], na.rm = TRUE)
        esd <- sd(sim_result[1, ], na.rm = TRUE)
        ese <- mean(sim_result[3, ], na.rm = TRUE)
        ci <- cbind(
          sim_result[1, ] - 1.96 * sim_result[3, ],
          sim_result[1, ] + 1.96 * sim_result[3, ]
        )
        cp <- mean(apply(ci, 1, function(v) mu.true >= v[1] & mu.true <= v[2]), na.rm = TRUE)

        results[row_num, col_start:(col_start + 3)] <- c(
          format_num(pe - mu.true), format_num(esd), format_num(ese), format_num(cp)
        )
      }
    }
  }

  # Create data frame with row labels
  results_df <- data.frame(
    Estimator = rep(estimator_names, length(n.vector)),
    results,
    stringsAsFactors = FALSE
  )
  colnames(results_df) <- c("", rep(c("Bias", "ESD", "ESE", "CP"), length(missing_rates)))

  return(list(
    data = results_df,
    mu.true = mu.true,
    setting = setting,
    scenario = scenario
  ))
}

#' Generate LaTeX table code
#'
#' @param table_data Output from generate_table()
#' @param mu_desc Description for caption
#' @param n.vector Sample sizes
#' @return LaTeX code string
generate_latex <- function(table_data, mu_desc, n.vector = c(500, 2000)) {
  setting <- table_data$setting
  scenario <- table_data$scenario
  mu.true <- table_data$mu.true
  df <- table_data$data

  # Build caption - include alpha values for all scenarios
  # Get alpha strings for caption (pass scenario to use correct alpha list)
  alpha_50 <- get_alpha_string(setting, "miss50", scenario)
  alpha_30 <- get_alpha_string(setting, "miss30", scenario)

  if (!is.null(alpha_50) && !is.null(alpha_30)) {
    caption <- paste0(
      "Estimation of $\\mu_0 = ", round(mu.true, 3), "$ is conducted under ", mu_desc, ". ",
      "The parameter vectors $\\bm{\\alpha}_0 = ", alpha_50, "$ and $", alpha_30, "$ ",
      "correspond to missingness rates of 50\\% and 30\\%, respectively."
    )
  } else {
    caption <- paste0(
      "Estimation of $\\mu_0 = ", round(mu.true, 3), "$ is conducted under ", mu_desc, "."
    )
  }

  # Start LaTeX code
  latex <- paste0(
    "\\begin{table}[ht]\n",
    "\\centering\n",
    "\\caption{", caption, "}\n",
    "\\begin{threeparttable}\n",
    "\\begin{tabularx}{\\textwidth}{l *{8}{>{\\centering\\arraybackslash}X}}\n",
    "\\toprule\n",
    "\\multicolumn{1}{c}{} & \\multicolumn{4}{c}{$50\\%$ missing} & \\multicolumn{4}{c}{$30\\%$ missing} \\\\\n",
    "\\cmidrule(lr){2-5} \\cmidrule(lr){6-9}\n",
    " & Bias & ESD & ESE & CP & Bias & ESD & ESE & CP \\\\\n"
  )

  # Add data rows for each sample size
  rows_per_n <- 8
  for (n_idx in seq_along(n.vector)) {
    n <- n.vector[n_idx]
    start_row <- (n_idx - 1) * rows_per_n + 1
    end_row <- n_idx * rows_per_n

    latex <- paste0(latex, "\\midrule\n")
    latex <- paste0(latex, "\\multicolumn{9}{c}{$n = ", n, "$} \\\\\n")
    latex <- paste0(latex, "\\midrule\n")

    for (i in start_row:end_row) {
      row_data <- df[i, ]
      latex <- paste0(latex,
                      row_data[1], " & ",
                      paste(row_data[2:9], collapse = " & "),
                      "\\\\\n")
    }
  }

  # Close table
  latex <- paste0(latex,
                  "\\bottomrule\n",
                  "\\end{tabularx}\n",
                  "\\end{threeparttable}\n",
                  "\\label{tab:", setting, "_", scenario, "}\n",
                  "\\end{table}\n")

  return(latex)
}

#------------------------------------------------------------------------------#
# Main Table Generation
#------------------------------------------------------------------------------#

#' Generate all tables for specified configurations
#'
#' @param configs List of table configurations
#' @param version Version string for result files (can be a named list for setting-specific versions)
#' @param output_file Output file path for LaTeX tables
generate_all_tables <- function(configs = TABLE_CONFIG,
                                version = "test15",
                                output_file = "simulation_tables_generated.tex") {

  cat("Generating simulation tables...\n\n")

  # Pre-compute mu.true for each unique setting using a very large sample
  # This ensures the same true mean is used across all scenarios of the same setting
  unique_settings <- unique(sapply(configs, function(x) x$setting))
  mu_true_map <- list()
  for (s in unique_settings) {
    mu_true_map[[s]] <- get_mu_true(s)
    cat("  mu.true for", s, "=", round(mu_true_map[[s]], 6), "\n")
  }
  cat("\n")

  all_latex <- ""

  for (i in seq_along(configs)) {
    config <- configs[[i]]
    setting <- config$setting
    scenario <- config$scenario

    # Get version for this setting (support named list or single string)
    if (is.list(version)) {
      current_version <- version[[setting]]
      if (is.null(current_version)) current_version <- version[["default"]]
    } else {
      current_version <- version
    }

    cat("Processing:", setting, "- Scenario", scenario, "(version:", current_version, ")\n")

    # Determine n.vector based on setting
    n.vector <- if (grepl("^Cho_", setting)) {
      c(500, 2000)  # Cho settings
    } else if (setting %in% c("setting3", "setting4", "setting13", "setting14", "setting15", "setting16", "setting17")) {
      c(500, 2000)  # Main paper settings
    } else {
      c(300, 1000)  # Legacy settings
    }

    # Generate table data (pass pre-computed mu.true for this setting)
    table_data <- tryCatch({
      generate_table(
        setting = setting,
        scenario = scenario,
        version = current_version,
        n.vector = n.vector,
        mu.true = mu_true_map[[setting]]
      )
    }, error = function(e) {
      cat("  ERROR:", e$message, "\n")
      return(NULL)
    })

    if (is.null(table_data)) next

    # Generate LaTeX
    latex <- generate_latex(
      table_data = table_data,
      mu_desc = config$mu
      _desc,
      n.vector = n.vector
    )

    all_latex <- paste0(all_latex, "\n", latex, "\n")
    cat("  Done.\n")
  }

  # Write to file
  writeLines(all_latex, output_file)
  cat("\nTables saved to:", output_file, "\n")

  # Automatically open the file
  if (Sys.info()["sysname"] == "Darwin") {
    system(paste("open", shQuote(output_file)))
  } else if (Sys.info()["sysname"] == "Windows") {
    shell.exec(output_file)
  } else {
    system(paste("xdg-open", shQuote(output_file)))
  }

  return(invisible(all_latex))
}

#------------------------------------------------------------------------------#
# Run Table Generation
#------------------------------------------------------------------------------#

# Uncomment to generate tables:
# Filter TABLE_CONFIG to only include setting3, setting4, and setting15 entries
# Use test47 for setting3/setting4, test22 for setting15

# setting3, setting4 with scenarios 7-1, 8-1, 9-1, 9-2, 9-3, 9-4 using test49
configs_test49 <- TABLE_CONFIG[sapply(TABLE_CONFIG, function(x) {
  x$setting %in% c("setting3", "setting4") &&
  x$scenario %in% c("7-1", "8-1", "9-1", "9-2", "9-3", "9-4", "7-2", "7-3", "7-4")
})]
generate_all_tables(configs = configs_test49, version = "test53")

# # Second: setting15 (excluding 8-2)
# configs_test22 <- TABLE_CONFIG[sapply(TABLE_CONFIG, function(x) {
#   x$setting %in% c("setting15") &&
#   !(x$scenario == "8-2")
# })]
# generate_all_tables(configs = configs_test22, version = "test22")

# # Third: setting15 with scenarios 8-2, 8-4 using test25
# configs_test25 <- TABLE_CONFIG[sapply(TABLE_CONFIG, function(x) {
#   x$setting %in% c("setting15") && x$scenario %in% c("8-2", "8-4")
# })]
# generate_all_tables(configs = configs_test25, version = "test25",
#                     output_file = "simulation_tables_8-2.tex")

# # Fourth: setting15 with scenario 8-3 using test26
# configs_test26 <- TABLE_CONFIG[sapply(TABLE_CONFIG, function(x) {
#   x$setting %in% c("setting15") && x$scenario %in% c("8-3")
# })]
# generate_all_tables(configs = configs_test26, version = "test26",
#                     output_file = "simulation_tables_8-3.tex")

# # Fifth: setting16, setting17 with scenario 8-2, 8-3, 8-4 using test26
# configs_test26b <- TABLE_CONFIG[sapply(TABLE_CONFIG, function(x) {
#   x$setting %in% c("setting16", "setting17") && x$scenario %in% c("8-2", "8-3", "8-4")
# })]
# generate_all_tables(configs = configs_test26b, version = "test26",
#                     output_file = "simulation_tables_8-2_8-3_8-4.tex")

# Generate Cho et al. tables with test22-6 version
# Cho_M1_gamma005, Cho_M1_gamma000 under cho1; Cho_M2_gamma005, Cho_M2_gamma000 under cho2
cho_configs <- TABLE_CONFIG[sapply(TABLE_CONFIG, function(x) {
  (x$setting %in% c("Cho_M1_gamma005", "Cho_M1_gamma000") && x$scenario == "cho1") ||
  (x$setting %in% c("Cho_M2_gamma005", "Cho_M2_gamma000") && x$scenario == "cho2")
})]
generate_all_tables(configs = cho_configs, version = "test22-6",
                    output_file = "simulation_tables_cho.tex")

# Or generate a single table:
# table_data <- generate_table("setting13", "7-1", version = "test15", n.vector = c(500, 2000))
# cat(generate_latex(table_data, "continuous outcome Y with correctly specified model"))
