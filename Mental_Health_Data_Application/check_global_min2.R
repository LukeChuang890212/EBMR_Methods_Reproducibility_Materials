#------------------------------------------------------------------------------#
# Check if alpha.hat is Global Minimizer for Each PS Model
# Version 2: Use the actual objective function from the gmm_fit
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
# PS specifications
#------------------------------------------------------------------------------#
full_ps_specifications <- list(
  formula.list = list(
    r ~ o(teacher_report) + father,
    r ~ o(teacher_report) + parent_report,
    r ~ father + parent_report  # MAR model
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

# Current initial values
init.list <- list(
  alpha_init.list0 = list(
    c(-0.3834341, 1.07945126, -0.1318883),
    c(-4.2670736, 5.9582123, -0.5972129),
    c(0.5, -0.2, -0.5)
  ),
  alpha_init.list1 = list(
    c(0.09019492, -1.511548, -0.3946616),
    c(-1.969056, 3.5823611, -1.37232734),
    c(0.5, -0.2, -0.5)
  )
)

health <- 0:1
subset_names <- c("health0", "health1")
model_names <- c("Model 1 (MNAR: father)", "Model 2 (MNAR: parent_report)", "Model 3 (MAR)")

#------------------------------------------------------------------------------#
# Check global minimum using the actual objective function from gmm_fit
#------------------------------------------------------------------------------#
cat("================================================================================\n")
cat("        CHECK GLOBAL MINIMUM FOR PS MODEL ALPHA ESTIMATES                       \n")
cat("================================================================================\n\n")

set.seed(12345)
n_starts <- 100

for (k in 1:length(health)) {
  subdat <- dat[dat$health == health[k], ]
  subdat$y <- subdat$teacher_report

  cat("--------------------------------------------------------------------------------\n")
  cat("Health =", health[k], "(", ifelse(k == 1, "Unexposed", "Exposed"), ")\n")
  cat("Sample size: n =", nrow(subdat), "\n")
  cat("--------------------------------------------------------------------------------\n\n")

  # Fit models with current initial values
  ps_specifications <- list(
    formula.list = full_ps_specifications$formula.list[1:3],
    h_x_names.list = full_ps_specifications$h_x_names.list[1:3],
    alpha_init.list = init.list[[k]][1:3],
    inv_link = full_ps_specifications$inv_link
  )

  ebmr <- EBMRAlgorithmFast$new("teacher_report", ps_specifications, subdat, W)

  for (j in 1:2) {  # Only check MNAR models (1 and 2)
    cat("  ", model_names[j], "\n")
    cat("  ", paste(rep("-", 70), collapse = ""), "\n")

    ps_fit <- ebmr$ps_fit.list[[j]]
    gmm_fit <- ps_fit$gmm_fit

    current_alpha <- ps_fit$coefficients
    current_obj <- gmm_fit$opt$value

    cat("    Current estimate:\n")
    cat("      alpha.hat =", round(current_alpha, 4), "\n")
    cat("      Objective value =", format(current_obj, scientific = TRUE, digits = 6), "\n")
    cat("      Convergence code =", gmm_fit$opt$convergence, "\n\n")

    # Get the actual objective function from gmm_fit
    obj_func <- gmm_fit$obj

    # Evaluate current solution
    current_obj_check <- obj_func(current_alpha)
    cat("    Objective at current alpha (re-evaluated):", format(current_obj_check, scientific = TRUE, digits = 6), "\n\n")

    # Try many random starting points
    cat("    Trying", n_starts, "random starting points...\n")

    # Generate random starting points with various ranges
    random_starts <- rbind(
      matrix(c(runif(n_starts/2, -3, 3), runif(n_starts/2, -5, 5), runif(n_starts/2, -3, 3)), ncol = 3),
      matrix(c(runif(n_starts/2, -5, 5), runif(n_starts/2, -10, 10), runif(n_starts/2, -5, 5)), ncol = 3)
    )

    # Add some structured starting points around the current solution
    structured_starts <- rbind(
      c(0, 0, 0),
      c(0, 1, 0),
      c(0, -1, 0),
      c(0, 2, 0),
      c(0, -2, 0),
      c(0, 5, 0),
      c(0, -5, 0),
      c(-1, 1, -0.5),
      c(1, -1, 0.5),
      c(-2, 3, -1),
      c(-3, 5, -1),
      c(-5, 8, -1),
      current_alpha,  # Current solution
      current_alpha + c(0.1, 0.1, 0.1),  # Perturbations around current
      current_alpha - c(0.1, 0.1, 0.1),
      current_alpha + c(0.5, 0.5, 0.5),
      current_alpha - c(0.5, 0.5, 0.5),
      current_alpha * 1.1,
      current_alpha * 0.9
    )

    all_starts <- rbind(random_starts, structured_starts)

    results <- data.frame(
      start_idx = integer(),
      final_obj = numeric(),
      alpha1 = numeric(),
      alpha2 = numeric(),
      alpha3 = numeric(),
      convergence = integer()
    )

    for (i in 1:nrow(all_starts)) {
      start <- all_starts[i, ]

      # Skip if starting point gives infinite objective
      start_obj <- tryCatch(obj_func(start), error = function(e) Inf)
      if (!is.finite(start_obj) || start_obj > 1e6) next

      tryCatch({
        opt_result <- optim(start, obj_func, method = "L-BFGS-B",
                            lower = rep(-10, 3), upper = rep(10, 3),
                            control = list(maxit = 2000))

        if (is.finite(opt_result$value)) {
          results <- rbind(results, data.frame(
            start_idx = i,
            final_obj = opt_result$value,
            alpha1 = opt_result$par[1],
            alpha2 = opt_result$par[2],
            alpha3 = opt_result$par[3],
            convergence = opt_result$convergence
          ))
        }
      }, error = function(e) {
        # Skip failed optimizations
      })
    }

    if (nrow(results) > 0) {
      # Sort by objective value
      results <- results[order(results$final_obj), ]

      # Show top 10 solutions
      cat("\n    Top 10 solutions found (out of", nrow(results), "successful runs):\n")
      cat(sprintf("    %5s  %15s  %12s  %12s  %12s  %5s\n",
                  "Rank", "Obj Value", "alpha[1]", "alpha[2]", "alpha[3]", "Conv"))
      cat("    ", paste(rep("-", 70), collapse = ""), "\n")

      for (i in 1:min(10, nrow(results))) {
        cat(sprintf("    %5d  %15.6e  %12.4f  %12.4f  %12.4f  %5d\n",
                    i, results$final_obj[i],
                    results$alpha1[i], results$alpha2[i], results$alpha3[i],
                    results$convergence[i]))
      }

      # Check if current solution is the best
      best_obj <- min(results$final_obj)
      best_idx <- which.min(results$final_obj)

      if (abs(current_obj_check - best_obj) / (abs(best_obj) + 1e-10) < 0.01) {
        cat("\n    >> Current solution IS the global minimum (within 1%)\n")
      } else if (current_obj_check < best_obj * 1.01) {
        cat("\n    >> Current solution is CLOSE to the best found\n")
      } else {
        cat("\n    >> WARNING: Found BETTER solution than current estimate!\n")
        cat("       Current obj:", format(current_obj_check, scientific = TRUE, digits = 6), "\n")
        cat("       Best obj:   ", format(best_obj, scientific = TRUE, digits = 6), "\n")
        cat("       Ratio:      ", format(current_obj_check / best_obj, digits = 4), "\n")

        # Show the better solution
        cat("       Better alpha: (", round(results$alpha1[best_idx], 4), ",",
            round(results$alpha2[best_idx], 4), ",",
            round(results$alpha3[best_idx], 4), ")\n")

        # Compare propensity scores
        pi_current <- 1 / (1 + exp(cbind(1, subdat$teacher_report, subdat[[ifelse(j==1, "father", "parent_report")]]) %*% current_alpha))
        pi_best <- 1 / (1 + exp(cbind(1, subdat$teacher_report, subdat[[ifelse(j==1, "father", "parent_report")]]) %*% c(results$alpha1[best_idx], results$alpha2[best_idx], results$alpha3[best_idx])))

        cat("\n       PS comparison:\n")
        cat("         Current PS range: [", round(min(pi_current), 4), ",", round(max(pi_current), 4), "]\n")
        cat("         Best PS range:    [", round(min(pi_best), 4), ",", round(max(pi_best), 4), "]\n")
        cat("         Correlation:      ", round(cor(as.vector(pi_current), as.vector(pi_best)), 4), "\n")
      }

      # Check for distinct local minima
      cat("\n    Checking for distinct local minima...\n")
      # Cluster solutions by objective value
      unique_objs <- unique(round(results$final_obj, 6))
      if (length(unique_objs) > 1) {
        cat("    Found", length(unique_objs), "distinct objective values:\n")
        for (obj_val in head(unique_objs, 5)) {
          subset_res <- results[round(results$final_obj, 6) == obj_val, ]
          cat(sprintf("      Obj = %.6e: %d solutions, avg alpha = (%.3f, %.3f, %.3f)\n",
                      obj_val, nrow(subset_res),
                      mean(subset_res$alpha1), mean(subset_res$alpha2), mean(subset_res$alpha3)))
        }
      }

    } else {
      cat("\n    No successful optimizations from random starts.\n")
    }

    cat("\n")
  }

  cat("\n")
}

cat("================================================================================\n")
cat("                        GLOBAL MINIMUM CHECK COMPLETE                           \n")
cat("================================================================================\n")
