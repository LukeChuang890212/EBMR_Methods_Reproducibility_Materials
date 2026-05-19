#------------------------------------------------------------------------------#
# Check if Case 2 solution is really the global minimum
#------------------------------------------------------------------------------#

library(EBMRalgorithmFast)
library(Matrix)

source("MHD_functions.R")

# Read and prepare data
original_data <- read.csv("data_application.csv")
percent <- original_data$Percentage
n <- 2486
class_n <- round(n * percent / 100)
dat <- gen_data(original_data, class_n, n)

cat("================================================================================\n")
cat("  CHECKING GLOBAL MINIMUM FOR CASE 2 (OVER-IDENTIFIED)\n")
cat("================================================================================\n\n")

# Case 2: h_x_names = c("health", "father", "parent_report")
# This gives h(X) = (1, health, father, parent_report) -> 4 moment conditions
# Model has 3 parameters

R <- dat$r
Y <- dat$teacher_report

# Construct h(X) matrix for Case 2
h_X <- cbind(
  intercept = 1,
  health = dat$health,
  father = dat$father,
  parent_report = dat$parent_report
)

# Define GMM objective function for Case 2
# g(alpha) = (R/pi - 1) * h(X)
# Objective = g_bar' * W * g_bar where g_bar = mean(g)

gmm_objective <- function(alpha, W_matrix = NULL) {
  eta <- alpha[1] + alpha[2] * Y + alpha[3] * dat$health
  pi_vec <- 1 / (1 + exp(eta))
  
  # Moment conditions
  g <- (R / pi_vec - 1) * h_X
  g_bar <- colMeans(g)
  
  # If W not provided, use identity
  if (is.null(W_matrix)) {
    W_matrix <- diag(4)
  }
  
  obj <- as.numeric(t(g_bar) %*% W_matrix %*% g_bar)
  return(obj)
}

# Compute optimal W matrix at a given alpha
compute_W <- function(alpha) {
  eta <- alpha[1] + alpha[2] * Y + alpha[3] * dat$health
  pi_vec <- 1 / (1 + exp(eta))
  g <- (R / pi_vec - 1) * h_X
  Omega <- t(g) %*% g / n
  tryCatch(solve(Omega), error = function(e) diag(4))
}

#------------------------------------------------------------------------------#
# Evaluate objective at different points
#------------------------------------------------------------------------------#

# Case 1 solution (from just-identified)
alpha_case1 <- c(-1.0529, 2.2007, -0.3050)

# Case 2 solution (from over-identified with default init)
alpha_case2 <- c(-0.2415, 0.1524, -0.1924)

cat("Candidate solutions:\n")
cat("  Case 1 solution: (", paste(round(alpha_case1, 4), collapse=", "), ")\n")
cat("  Case 2 solution: (", paste(round(alpha_case2, 4), collapse=", "), ")\n\n")

#------------------------------------------------------------------------------#
# Evaluate with identity W
#------------------------------------------------------------------------------#
cat("--------------------------------------------------------------------------------\n")
cat("Using Identity Weighting Matrix (W = I)\n")
cat("--------------------------------------------------------------------------------\n\n")

obj_case1_I <- gmm_objective(alpha_case1, diag(4))
obj_case2_I <- gmm_objective(alpha_case2, diag(4))

cat(sprintf("  Objective at Case 1 solution: %.6e\n", obj_case1_I))
cat(sprintf("  Objective at Case 2 solution: %.6e\n", obj_case2_I))
cat(sprintf("  Winner: %s\n\n", ifelse(obj_case1_I < obj_case2_I, "Case 1", "Case 2")))

#------------------------------------------------------------------------------#
# Evaluate with optimal W computed at Case 2 solution
#------------------------------------------------------------------------------#
cat("--------------------------------------------------------------------------------\n")
cat("Using Optimal W (computed at Case 2 solution)\n")
cat("--------------------------------------------------------------------------------\n\n")

W_opt <- compute_W(alpha_case2)

obj_case1_W <- gmm_objective(alpha_case1, W_opt)
obj_case2_W <- gmm_objective(alpha_case2, W_opt)

cat(sprintf("  Objective at Case 1 solution: %.6e\n", obj_case1_W))
cat(sprintf("  Objective at Case 2 solution: %.6e\n", obj_case2_W))
cat(sprintf("  Winner: %s\n\n", ifelse(obj_case1_W < obj_case2_W, "Case 1", "Case 2")))

#------------------------------------------------------------------------------#
# Try optimizing from Case 1 solution as initial value
#------------------------------------------------------------------------------#
cat("--------------------------------------------------------------------------------\n")
cat("Re-optimizing Case 2 starting from Case 1 solution\n")
cat("--------------------------------------------------------------------------------\n\n")

# Two-step GMM: first with identity W, then with optimal W
result_step1 <- optim(
  alpha_case1, 
  function(a) gmm_objective(a, diag(4)),
  method = "L-BFGS-B",
  lower = c(-20, -20, -20),
  upper = c(20, 20, 20)
)

cat("Step 1 (Identity W):\n")
cat(sprintf("  Starting from: (%s)\n", paste(round(alpha_case1, 4), collapse=", ")))
cat(sprintf("  Converged to:  (%s)\n", paste(round(result_step1$par, 4), collapse=", ")))
cat(sprintf("  Objective:     %.6e\n", result_step1$value))
cat(sprintf("  Convergence:   %d\n\n", result_step1$convergence))

# Step 2 with optimal W
W_step1 <- compute_W(result_step1$par)
result_step2 <- optim(
  result_step1$par,
  function(a) gmm_objective(a, W_step1),
  method = "L-BFGS-B",
  lower = c(-20, -20, -20),
  upper = c(20, 20, 20)
)

cat("Step 2 (Optimal W):\n")
cat(sprintf("  Starting from: (%s)\n", paste(round(result_step1$par, 4), collapse=", ")))
cat(sprintf("  Converged to:  (%s)\n", paste(round(result_step2$par, 4), collapse=", ")))
cat(sprintf("  Objective:     %.6e\n", result_step2$value))
cat(sprintf("  Convergence:   %d\n\n", result_step2$convergence))

#------------------------------------------------------------------------------#
# Grid search to find global minimum
#------------------------------------------------------------------------------#
cat("--------------------------------------------------------------------------------\n")
cat("Grid Search for Global Minimum\n")
cat("--------------------------------------------------------------------------------\n\n")

# Search over a grid of starting points
intercept_grid <- seq(-3, 1, by = 1)
oy_grid <- seq(-1, 3, by = 1)
health_grid <- seq(-1, 1, by = 0.5)

best_obj <- Inf
best_alpha <- NULL
all_results <- list()

cat("Searching over", length(intercept_grid) * length(oy_grid) * length(health_grid), "starting points...\n\n")

for (int in intercept_grid) {
  for (oy in oy_grid) {
    for (h in health_grid) {
      init <- c(int, oy, h)
      
      # Two-step GMM
      res1 <- tryCatch({
        optim(init, function(a) gmm_objective(a, diag(4)),
              method = "L-BFGS-B",
              lower = c(-20, -20, -20),
              upper = c(20, 20, 20))
      }, error = function(e) list(par = init, value = Inf))
      
      W_temp <- tryCatch(compute_W(res1$par), error = function(e) diag(4))
      
      res2 <- tryCatch({
        optim(res1$par, function(a) gmm_objective(a, W_temp),
              method = "L-BFGS-B",
              lower = c(-20, -20, -20),
              upper = c(20, 20, 20))
      }, error = function(e) list(par = res1$par, value = Inf))
      
      if (res2$value < best_obj) {
        best_obj <- res2$value
        best_alpha <- res2$par
      }
    }
  }
}

cat(sprintf("Best solution found:\n"))
cat(sprintf("  Alpha:     (%s)\n", paste(round(best_alpha, 4), collapse=", ")))
cat(sprintf("  Objective: %.6e\n\n", best_obj))

#------------------------------------------------------------------------------#
# Compare all solutions
#------------------------------------------------------------------------------#
cat("================================================================================\n")
cat("  SUMMARY COMPARISON\n")
cat("================================================================================\n\n")

# Recompute objectives with same W for fair comparison
W_final <- compute_W(best_alpha)

obj_case1_final <- gmm_objective(alpha_case1, W_final)
obj_case2_final <- gmm_objective(alpha_case2, W_final)
obj_best_final <- gmm_objective(best_alpha, W_final)

cat(sprintf("  %-40s %12s %12s\n", "Solution", "Alpha", "Objective"))
cat("  ", paste(rep("-", 70), collapse = ""), "\n")
cat(sprintf("  %-40s (%s) %12.6e\n", 
            "Case 1 (just-identified)", 
            paste(round(alpha_case1, 4), collapse=", "),
            obj_case1_final))
cat(sprintf("  %-40s (%s) %12.6e\n", 
            "Case 2 (package solution)", 
            paste(round(alpha_case2, 4), collapse=", "),
            obj_case2_final))
cat(sprintf("  %-40s (%s) %12.6e\n", 
            "Grid search best", 
            paste(round(best_alpha, 4), collapse=", "),
            obj_best_final))
cat("  ", paste(rep("-", 70), collapse = ""), "\n\n")

# Show moment values for each solution
cat("Moment conditions at each solution:\n\n")

show_moments <- function(alpha, label) {
  eta <- alpha[1] + alpha[2] * Y + alpha[3] * dat$health
  pi_vec <- 1 / (1 + exp(eta))
  g <- (R / pi_vec - 1) * h_X
  g_bar <- colMeans(g)
  cat(sprintf("  %s:\n", label))
  cat(sprintf("    E[(R/pi-1)*1]:            %12.6e\n", g_bar[1]))
  cat(sprintf("    E[(R/pi-1)*health]:       %12.6e\n", g_bar[2]))
  cat(sprintf("    E[(R/pi-1)*father]:       %12.6e\n", g_bar[3]))
  cat(sprintf("    E[(R/pi-1)*parent_report]:%12.6e\n\n", g_bar[4]))
}

show_moments(alpha_case1, "Case 1 solution")
show_moments(alpha_case2, "Case 2 solution")
show_moments(best_alpha, "Grid search best")

cat("================================================================================\n")
