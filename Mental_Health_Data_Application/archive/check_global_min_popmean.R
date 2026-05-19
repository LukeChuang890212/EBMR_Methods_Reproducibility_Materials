#------------------------------------------------------------------------------#
# Check if alpha.hat and nu.hat are Global Minimizers for Population Mean Analysis
#------------------------------------------------------------------------------#

library(EBMRalgorithmFast)
library(tidyverse)
library(numDeriv)

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
# PS specifications (same as Population_Mean_Analysis.R)
#------------------------------------------------------------------------------#
full_ps_specifications <- list(
  formula.list = list(
    r ~ o(teacher_report) + health + father,          # Model 1
    r ~ o(teacher_report) + health + parent_report,   # Model 2
    r ~ o(teacher_report) + father + parent_report    # Model 3
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

alpha_init.list <- list(
  c(0, 1, 0, 0),
  c(0, 1, 0, 0),
  c(0, 1, 0, 0)
)

# h_nu function (4 terms as per current Population_Mean_Analysis.R)
h_nu <- function(data) cbind(
  health = data$health,
  father = data$father,
  parent_report = data$parent_report,
  hf = data$health * data$father
)

model_names <- c(
  "Model 1 (o(y) + health + father)",
  "Model 2 (o(y) + health + parent_report)",
  "Model 3 (o(y) + father + parent_report)"
)

#------------------------------------------------------------------------------#
# Fit EBMR with full model (all 3 PS models)
#------------------------------------------------------------------------------#
cat("================================================================================\n")
cat("  CHECK GLOBAL MINIMUM FOR POPULATION MEAN ANALYSIS                             \n")
cat("================================================================================\n\n")

cat("Sample size: n =", nrow(dat), "\n")
cat("Missing rate:", round(100 * (1 - mean(dat$r)), 1), "%\n\n")

ps_specifications <- list(
  formula.list = full_ps_specifications$formula.list[1:3],
  h_x_names.list = full_ps_specifications$h_x_names.list[1:3],
  alpha_init.list = alpha_init.list[1:3],
  inv_link = full_ps_specifications$inv_link
)

ebmr <- EBMRAlgorithmFast$new("teacher_report", ps_specifications, dat, W)
result <- ebmr$EBMR_IPW(h_nu = h_nu, true_ps = NULL)

#------------------------------------------------------------------------------#
# Check alpha estimates for each PS model
#------------------------------------------------------------------------------#
cat("================================================================================\n")
cat("                    CHECKING ALPHA ESTIMATES (PS MODELS)                        \n")
cat("================================================================================\n\n")

set.seed(12345)
n_starts <- 100

for (j in 1:3) {
  cat("--------------------------------------------------------------------------------\n")
  cat(model_names[j], "\n")
  cat("--------------------------------------------------------------------------------\n\n")

  ps_fit <- ebmr$ps_fit.list[[j]]
  gmm_fit <- ps_fit$gmm_fit

  current_alpha <- ps_fit$coefficients
  current_obj <- gmm_fit$opt$value

  cat("  Current estimate:\n")
  cat("    alpha.hat =", round(current_alpha, 4), "\n")
  cat("    Objective value =", format(current_obj, scientific = TRUE, digits = 6), "\n")
  cat("    Convergence code =", gmm_fit$opt$convergence, "\n")

  # Get the actual objective function
  obj_func <- gmm_fit$obj

  # Re-evaluate current solution
  current_obj_check <- obj_func(current_alpha)
  cat("    Objective at current (re-evaluated):", format(current_obj_check, scientific = TRUE, digits = 6), "\n\n")

  # Try many random starting points
  cat("  Trying", n_starts, "random starting points...\n")

  alpha_dim <- length(current_alpha)

  # Generate random starting points
  random_starts <- matrix(
    c(runif(n_starts, -3, 3),     # intercept
      runif(n_starts, -5, 5),     # o(y)
      runif(n_starts, -3, 3),     # x1
      runif(n_starts, -3, 3)),    # x2
    ncol = alpha_dim
  )

  # Add structured starting points
  structured_starts <- rbind(
    rep(0, alpha_dim),
    c(0, 1, 0, 0),
    c(0, -1, 0, 0),
    c(0, 2, 0, 0),
    c(0, -2, 0, 0),
    c(0, 5, 0, 0),
    c(-1, 1, -0.5, -0.5),
    c(1, -1, 0.5, 0.5),
    current_alpha,
    current_alpha * 1.1,
    current_alpha * 0.9,
    current_alpha + c(0.1, 0.1, 0.1, 0.1),
    current_alpha - c(0.1, 0.1, 0.1, 0.1)
  )

  all_starts <- rbind(random_starts, structured_starts)

  results_df <- data.frame(
    start_idx = integer(),
    final_obj = numeric(),
    alpha1 = numeric(),
    alpha2 = numeric(),
    alpha3 = numeric(),
    alpha4 = numeric(),
    convergence = integer()
  )

  for (i in 1:nrow(all_starts)) {
    start <- all_starts[i, ]

    start_obj <- tryCatch(obj_func(start), error = function(e) Inf)
    if (!is.finite(start_obj) || start_obj > 1e6) next

    tryCatch({
      opt_result <- optim(start, obj_func, method = "L-BFGS-B",
                          lower = rep(-10, alpha_dim), upper = rep(10, alpha_dim),
                          control = list(maxit = 2000))

      if (is.finite(opt_result$value)) {
        results_df <- rbind(results_df, data.frame(
          start_idx = i,
          final_obj = opt_result$value,
          alpha1 = opt_result$par[1],
          alpha2 = opt_result$par[2],
          alpha3 = opt_result$par[3],
          alpha4 = opt_result$par[4],
          convergence = opt_result$convergence
        ))
      }
    }, error = function(e) {})
  }

  if (nrow(results_df) > 0) {
    results_df <- results_df[order(results_df$final_obj), ]

    cat("\n  Top 10 solutions found (out of", nrow(results_df), "successful runs):\n")
    cat(sprintf("  %5s  %15s  %10s  %10s  %10s  %10s  %5s\n",
                "Rank", "Obj Value", "alpha[1]", "alpha[2]", "alpha[3]", "alpha[4]", "Conv"))
    cat("  ", paste(rep("-", 80), collapse = ""), "\n")

    for (i in 1:min(10, nrow(results_df))) {
      cat(sprintf("  %5d  %15.6e  %10.4f  %10.4f  %10.4f  %10.4f  %5d\n",
                  i, results_df$final_obj[i],
                  results_df$alpha1[i], results_df$alpha2[i],
                  results_df$alpha3[i], results_df$alpha4[i],
                  results_df$convergence[i]))
    }

    best_obj <- min(results_df$final_obj)
    best_idx <- which.min(results_df$final_obj)

    if (abs(current_obj_check - best_obj) / (abs(best_obj) + 1e-10) < 0.01) {
      cat("\n  >> Current solution IS the global minimum (within 1%)\n")
    } else if (current_obj_check <= best_obj * 1.01) {
      cat("\n  >> Current solution is CLOSE to the best found\n")
    } else {
      cat("\n  >> WARNING: Found BETTER solution than current estimate!\n")
      cat("     Current obj:", format(current_obj_check, scientific = TRUE, digits = 6), "\n")
      cat("     Best obj:   ", format(best_obj, scientific = TRUE, digits = 6), "\n")
      cat("     Ratio:      ", format(current_obj_check / best_obj, digits = 4), "\n")
      cat("     Better alpha: (", round(results_df$alpha1[best_idx], 4), ",",
          round(results_df$alpha2[best_idx], 4), ",",
          round(results_df$alpha3[best_idx], 4), ",",
          round(results_df$alpha4[best_idx], 4), ")\n")
    }

    # Check for distinct local minima
    unique_objs <- unique(round(results_df$final_obj, 6))
    if (length(unique_objs) > 1) {
      cat("\n  Found", length(unique_objs), "distinct objective values (potential local minima)\n")
    }
  } else {
    cat("\n  No successful optimizations from random starts.\n")
  }

  # Check PS range
  ps_range <- range(ps_fit$fitted.values)
  cat("\n  PS range: [", round(ps_range[1], 4), ",", round(ps_range[2], 4), "]\n")

  extreme_low <- sum(ps_fit$fitted.values < 0.01)
  extreme_high <- sum(ps_fit$fitted.values > 0.99)
  if (extreme_low > 0 || extreme_high > 0) {
    cat("  WARNING: Extreme PS values -", extreme_low, "below 0.01,", extreme_high, "above 0.99\n")
  }

  cat("\n")
}

#------------------------------------------------------------------------------#
# Check nu estimates (ensemble step)
#------------------------------------------------------------------------------#
cat("================================================================================\n")
cat("                    CHECKING NU ESTIMATES (ENSEMBLE STEP)                       \n")
cat("================================================================================\n\n")

ensemble_fit <- result$ensemble_fit
ensemble_gmm <- ensemble_fit$gmm_fit

current_nu <- result$nu.hat
current_w <- result$w.hat

cat("Current estimate:\n")
cat("  nu.hat =", round(current_nu, 4), "\n")
cat("  w.hat  =", round(current_w, 4), "(", round(current_w * 100, 1), "%)\n")
cat("  mu_ipw =", round(result$mu_ipw, 4), "(SE:", round(result$se_ipw, 4), ")\n\n")

if (!is.null(ensemble_gmm$opt)) {
  current_obj_nu <- ensemble_gmm$opt$value
  cat("  Objective value =", format(current_obj_nu, scientific = TRUE, digits = 6), "\n")
  cat("  Convergence code =", ensemble_gmm$opt$convergence, "\n")

  # Get ensemble objective function
  obj_func_nu <- ensemble_gmm$obj

  current_obj_nu_check <- obj_func_nu(current_nu)
  cat("  Objective at current (re-evaluated):", format(current_obj_nu_check, scientific = TRUE, digits = 6), "\n\n")

  # Try random starting points for nu
  cat("Trying", n_starts, "random starting points for nu...\n")

  J <- length(current_nu)

  # Generate starting points (nu should sum to ~1 in absolute value for reasonable weights)
  random_starts_nu <- matrix(runif(n_starts * J, -2, 2), ncol = J)

  # Normalize some starts
  for (i in 1:(n_starts/2)) {
    random_starts_nu[i, ] <- random_starts_nu[i, ] / sum(abs(random_starts_nu[i, ]))
  }

  # Add structured starting points
  structured_starts_nu <- rbind(
    c(1, 0, 0),
    c(0, 1, 0),
    c(0, 0, 1),
    c(1/3, 1/3, 1/3),
    c(0.5, 0.25, 0.25),
    c(0.25, 0.5, 0.25),
    c(0.25, 0.25, 0.5),
    c(0.9, 0.05, 0.05),
    c(0.05, 0.9, 0.05),
    c(0.05, 0.05, 0.9),
    current_nu,
    current_nu * 1.1,
    current_nu * 0.9,
    -current_nu
  )

  all_starts_nu <- rbind(random_starts_nu, structured_starts_nu)

  results_nu <- data.frame(
    start_idx = integer(),
    final_obj = numeric(),
    nu1 = numeric(),
    nu2 = numeric(),
    nu3 = numeric(),
    convergence = integer()
  )

  for (i in 1:nrow(all_starts_nu)) {
    start <- all_starts_nu[i, ]

    start_obj <- tryCatch(obj_func_nu(start), error = function(e) Inf)
    if (!is.finite(start_obj) || start_obj > 1e6) next

    tryCatch({
      opt_result <- optim(start, obj_func_nu, method = "L-BFGS-B",
                          lower = rep(-10, J), upper = rep(10, J),
                          control = list(maxit = 2000))

      if (is.finite(opt_result$value)) {
        results_nu <- rbind(results_nu, data.frame(
          start_idx = i,
          final_obj = opt_result$value,
          nu1 = opt_result$par[1],
          nu2 = opt_result$par[2],
          nu3 = opt_result$par[3],
          convergence = opt_result$convergence
        ))
      }
    }, error = function(e) {})
  }

  if (nrow(results_nu) > 0) {
    results_nu <- results_nu[order(results_nu$final_obj), ]

    # Compute w for each solution
    results_nu$w1 <- results_nu$nu1^2 / (results_nu$nu1^2 + results_nu$nu2^2 + results_nu$nu3^2)
    results_nu$w2 <- results_nu$nu2^2 / (results_nu$nu1^2 + results_nu$nu2^2 + results_nu$nu3^2)
    results_nu$w3 <- results_nu$nu3^2 / (results_nu$nu1^2 + results_nu$nu2^2 + results_nu$nu3^2)

    cat("\nTop 10 solutions found (out of", nrow(results_nu), "successful runs):\n")
    cat(sprintf("%5s  %15s  %10s  %10s  %10s  %8s  %8s  %8s\n",
                "Rank", "Obj Value", "nu[1]", "nu[2]", "nu[3]", "w[1]%", "w[2]%", "w[3]%"))
    cat(paste(rep("-", 90), collapse = ""), "\n")

    for (i in 1:min(10, nrow(results_nu))) {
      cat(sprintf("%5d  %15.6e  %10.4f  %10.4f  %10.4f  %8.1f  %8.1f  %8.1f\n",
                  i, results_nu$final_obj[i],
                  results_nu$nu1[i], results_nu$nu2[i], results_nu$nu3[i],
                  results_nu$w1[i] * 100, results_nu$w2[i] * 100, results_nu$w3[i] * 100))
    }

    best_obj_nu <- min(results_nu$final_obj)
    best_idx_nu <- which.min(results_nu$final_obj)

    if (abs(current_obj_nu_check - best_obj_nu) / (abs(best_obj_nu) + 1e-10) < 0.01) {
      cat("\n>> Current nu solution IS the global minimum (within 1%)\n")
    } else if (current_obj_nu_check <= best_obj_nu * 1.01) {
      cat("\n>> Current nu solution is CLOSE to the best found\n")
    } else {
      cat("\n>> WARNING: Found BETTER nu solution than current estimate!\n")
      cat("   Current obj:", format(current_obj_nu_check, scientific = TRUE, digits = 6), "\n")
      cat("   Best obj:   ", format(best_obj_nu, scientific = TRUE, digits = 6), "\n")
      cat("   Ratio:      ", format(current_obj_nu_check / best_obj_nu, digits = 4), "\n")
      cat("   Better nu: (", round(results_nu$nu1[best_idx_nu], 4), ",",
          round(results_nu$nu2[best_idx_nu], 4), ",",
          round(results_nu$nu3[best_idx_nu], 4), ")\n")
      cat("   Better w:  (", round(results_nu$w1[best_idx_nu] * 100, 1), "%,",
          round(results_nu$w2[best_idx_nu] * 100, 1), "%,",
          round(results_nu$w3[best_idx_nu] * 100, 1), "%)\n")
    }

    # Check for distinct local minima
    unique_objs_nu <- unique(round(results_nu$final_obj, 6))
    if (length(unique_objs_nu) > 1) {
      cat("\nFound", length(unique_objs_nu), "distinct objective values for nu:\n")
      for (obj_val in head(unique_objs_nu, 5)) {
        subset_res <- results_nu[round(results_nu$final_obj, 6) == obj_val, ]
        cat(sprintf("  Obj = %.6e: %d solutions, avg w = (%.1f%%, %.1f%%, %.1f%%)\n",
                    obj_val, nrow(subset_res),
                    mean(subset_res$w1) * 100, mean(subset_res$w2) * 100, mean(subset_res$w3) * 100))
      }
    }
  } else {
    cat("\nNo successful optimizations for nu from random starts.\n")
  }
}

cat("\n================================================================================\n")
cat("                        GLOBAL MINIMUM CHECK COMPLETE                           \n")
cat("================================================================================\n")
