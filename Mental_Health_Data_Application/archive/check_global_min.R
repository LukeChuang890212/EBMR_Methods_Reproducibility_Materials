#------------------------------------------------------------------------------#
# Check if alpha.hat is Global Minimizer for Each PS Model
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
# Function to create objective function for a given model
#------------------------------------------------------------------------------#
create_obj_function <- function(subdat, model_idx, inv_link) {
  r <- subdat$r
  y <- subdat$teacher_report
  n <- nrow(subdat)

  if (model_idx == 1) {
    # Model 1: r ~ o(y) + father
    x <- subdat$father
    h_x <- cbind(1, subdat$father, subdat$parent_report)
  } else if (model_idx == 2) {
    # Model 2: r ~ o(y) + parent_report
    x <- subdat$parent_report
    h_x <- cbind(1, subdat$father, subdat$parent_report)
  } else {
    # Model 3: MAR - use glm, not GMM
    return(NULL)
  }

  # Design matrix: intercept, o(y), x
  design_matrix <- cbind(1, y, x)

  # Objective function: quadratic form of moment conditions
  obj <- function(alpha) {
    eta <- design_matrix %*% alpha
    pi_vec <- inv_link(eta)

    # Check for invalid propensity scores
    if (any(pi_vec <= 0) || any(pi_vec >= 1) || any(!is.finite(pi_vec))) {
      return(1e10)
    }

    g_matrix <- (r / pi_vec - 1) * h_x
    G <- colMeans(g_matrix)

    # First stage: identity weight
    # Later stages use optimal weight
    tryCatch({
      W_mat <- solve(t(g_matrix) %*% g_matrix / n)
      value <- as.numeric(t(G) %*% W_mat %*% G)
      if (!is.finite(value)) return(1e10)
      return(value)
    }, error = function(e) {
      return(1e10)
    })
  }

  return(obj)
}

#------------------------------------------------------------------------------#
# Check global minimum with multiple starting points
#------------------------------------------------------------------------------#
cat("================================================================================\n")
cat("        CHECK GLOBAL MINIMUM FOR PS MODEL ALPHA ESTIMATES                       \n")
cat("================================================================================\n\n")

# Generate random starting points
set.seed(12345)
n_starts <- 50

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
    current_alpha <- ps_fit$coefficients
    current_obj <- ps_fit$gmm_fit$opt$value

    cat("    Current estimate:\n")
    cat("      alpha.hat =", round(current_alpha, 4), "\n")
    cat("      Objective value =", format(current_obj, scientific = TRUE, digits = 6), "\n")
    cat("      Convergence code =", ps_fit$gmm_fit$opt$convergence, "\n\n")

    # Create objective function for this model
    obj_func <- create_obj_function(subdat, j, full_ps_specifications$inv_link)

    if (!is.null(obj_func)) {
      # Try many random starting points
      cat("    Trying", n_starts, "random starting points...\n")

      # Generate random starting points
      random_starts <- matrix(
        c(runif(n_starts, -5, 5),    # intercept
          runif(n_starts, -10, 10),  # coefficient on y
          runif(n_starts, -5, 5)),   # coefficient on x
        ncol = 3
      )

      # Add some structured starting points
      structured_starts <- rbind(
        c(0, 0, 0),
        c(0, 1, 0),
        c(0, -1, 0),
        c(0, 2, 0),
        c(0, -2, 0),
        c(-1, 1, -0.5),
        c(1, -1, 0.5),
        c(-2, 3, -1),
        c(-3, 5, -1),
        c(-5, 8, -1),
        init.list[[k]][[j]]  # Current init
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

        tryCatch({
          opt_result <- optim(start, obj_func, method = "L-BFGS-B",
                              lower = rep(-10, 3), upper = rep(10, 3),
                              control = list(maxit = 1000))

          results <- rbind(results, data.frame(
            start_idx = i,
            final_obj = opt_result$value,
            alpha1 = opt_result$par[1],
            alpha2 = opt_result$par[2],
            alpha3 = opt_result$par[3],
            convergence = opt_result$convergence
          ))
        }, error = function(e) {
          # Skip failed optimizations
        })
      }

      # Sort by objective value
      results <- results[order(results$final_obj), ]

      # Show top 5 solutions
      cat("\n    Top 5 solutions found:\n")
      cat(sprintf("    %5s  %15s  %12s  %12s  %12s  %5s\n",
                  "Rank", "Obj Value", "alpha[1]", "alpha[2]", "alpha[3]", "Conv"))
      cat("    ", paste(rep("-", 70), collapse = ""), "\n")

      for (i in 1:min(5, nrow(results))) {
        cat(sprintf("    %5d  %15.6e  %12.4f  %12.4f  %12.4f  %5d\n",
                    i, results$final_obj[i],
                    results$alpha1[i], results$alpha2[i], results$alpha3[i],
                    results$convergence[i]))
      }

      # Check if current solution is the best
      best_obj <- min(results$final_obj)
      if (abs(current_obj - best_obj) < 1e-6) {
        cat("\n    >> Current solution IS the global minimum (or very close)\n")
      } else if (current_obj < best_obj) {
        cat("\n    >> Current solution is BETTER than all random starts\n")
      } else {
        cat("\n    >> WARNING: Found BETTER solution than current estimate!\n")
        cat("       Current obj:", format(current_obj, scientific = TRUE, digits = 6), "\n")
        cat("       Best obj:   ", format(best_obj, scientific = TRUE, digits = 6), "\n")
        cat("       Improvement:", format(current_obj - best_obj, scientific = TRUE, digits = 6), "\n")

        # Show the better solution
        best_idx <- which.min(results$final_obj)
        cat("       Better alpha: (", round(results$alpha1[best_idx], 4), ",",
            round(results$alpha2[best_idx], 4), ",",
            round(results$alpha3[best_idx], 4), ")\n")
      }

      # Check for multiple local minima
      unique_solutions <- results %>%
        mutate(obj_rounded = round(final_obj, 8)) %>%
        group_by(obj_rounded) %>%
        summarize(
          count = n(),
          avg_alpha1 = mean(alpha1),
          avg_alpha2 = mean(alpha2),
          avg_alpha3 = mean(alpha3),
          .groups = "drop"
        ) %>%
        filter(count >= 2)  # Only show solutions found multiple times

      if (nrow(unique_solutions) > 1) {
        cat("\n    Multiple local minima detected:\n")
        print(unique_solutions, n = 10)
      }
    }

    cat("\n")
  }

  cat("\n")
}

cat("================================================================================\n")
cat("                        GLOBAL MINIMUM CHECK COMPLETE                           \n")
cat("================================================================================\n")
