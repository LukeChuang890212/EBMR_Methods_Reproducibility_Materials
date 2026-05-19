#------------------------------------------------------------------------------#
# Compare Analytical SE vs Bootstrap SE for Both Approaches
#------------------------------------------------------------------------------#

library(EBMRalgorithmFast)
library(tidyverse)
library(numDeriv)

source("MHD_functions.R")

original_data <- read.csv("data_application.csv")
percent <- original_data$Percentage
n <- 2486
class_n <- round(n * percent / 100)
dat <- gen_data(original_data, class_n, n)

cat("================================================================================\n")
cat("     COMPARE ANALYTICAL SE VS BOOTSTRAP SE                                     \n")
cat("================================================================================\n\n")

# Focus on Health = 1, Model 1
health_val <- 1
subdat <- dat[dat$health == health_val, ]
n_sub <- nrow(subdat)

#------------------------------------------------------------------------------#
# CURRENT APPROACH (using EBMRalgorithmFast)
#------------------------------------------------------------------------------#
cat("================================================================================\n")
cat("CURRENT APPROACH: r ~ o(y) + father\n")
cat("================================================================================\n\n")

ps_specifications <- list(
  formula.list = list(r ~ o(teacher_report) + father),
  h_x_names.list = list(c("father", "parent_report")),
  alpha_init.list = list(c(0.09, -1.5, -0.4)),
  inv_link = function(eta) 1 / (1 + exp(eta))
)

W <- function(g.matrix) {
  solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
}

ebmr <- EBMRAlgorithmFast$new("teacher_report", ps_specifications, subdat, W)
ps_fit <- ebmr$ps_fit.list[[1]]

alpha_hat <- ps_fit$coefficients
se_analytical <- ps_fit$se

cat("Point estimates:", round(alpha_hat, 4), "\n")
cat("Analytical SE:  ", round(se_analytical, 4), "\n\n")

# Compute Jacobian column for alpha_y to check
gmm_fit <- ps_fit$gmm_fit
Gamma_hat <- gmm_fit$Gamma.hat

cat("Jacobian Gamma.hat (column 2 = d(E[g])/d(alpha_y)):\n")
print(round(Gamma_hat, 6))
cat("\n")

cat("Note: Column 2 has small values because derivative involves y,\n")
cat("and y=0 for most observed cases.\n\n")

# Bootstrap for comparison
r_vec <- subdat$r
y_vec <- subdat$teacher_report
x_vec <- subdat$father
h_x <- cbind(1, subdat$father, subdat$parent_report)
design_matrix <- cbind(1, y_vec, x_vec)
inv_link <- function(eta) 1 / (1 + exp(eta))

set.seed(12345)
n_boot <- 500
boot_curr <- matrix(NA, n_boot, 3)

for (b in 1:n_boot) {
  idx <- sample(1:n_sub, replace = TRUE)
  r_b <- r_vec[idx]; design_b <- design_matrix[idx,]; h_x_b <- h_x[idx,]

  obj_b <- function(alpha) {
    eta <- design_b %*% alpha
    pi_vec <- inv_link(eta)
    if (any(pi_vec <= 0.001) || any(pi_vec >= 0.999)) return(1e10)
    rw <- r_b / as.vector(pi_vec)
    g.mat <- (rw - 1) * h_x_b
    G <- colMeans(g.mat)
    W_mat <- tryCatch(solve(crossprod(g.mat) / n_sub), error = function(e) diag(3))
    as.numeric(t(G) %*% W_mat %*% G)
  }

  tryCatch({
    opt <- optim(alpha_hat, obj_b, method = "L-BFGS-B",
                 lower = rep(-10, 3), upper = rep(10, 3), control = list(maxit = 1000))
    if (opt$convergence == 0 && all(abs(opt$par) < 9.9)) boot_curr[b,] <- opt$par
  }, error = function(e) {})
}

se_boot_curr <- apply(boot_curr[complete.cases(boot_curr),], 2, sd)

cat("SE Comparison:\n")
cat(sprintf("  %-12s  Analytical  Bootstrap  Ratio\n", ""))
cat(sprintf("  %-12s  %9.4f  %9.4f  %5.3f\n", "Intercept", se_analytical[1], se_boot_curr[1], se_analytical[1]/se_boot_curr[1]))
cat(sprintf("  %-12s  %9.4f  %9.4f  %5.3f\n", "o(y)", se_analytical[2], se_boot_curr[2], se_analytical[2]/se_boot_curr[2]))
cat(sprintf("  %-12s  %9.4f  %9.4f  %5.3f\n", "father", se_analytical[3], se_boot_curr[3], se_analytical[3]/se_boot_curr[3]))

cat("\n")

#------------------------------------------------------------------------------#
# CENTERED APPROACH (manual implementation)
#------------------------------------------------------------------------------#
cat("================================================================================\n")
cat("CENTERED APPROACH: r ~ o(y - y_bar) + father\n")
cat("================================================================================\n\n")

y_bar <- mean(y_vec[r_vec == 1])
y_centered <- y_vec - y_bar
design_centered <- cbind(1, y_centered, x_vec)

cat("y_bar =", round(y_bar, 4), "\n\n")

# Fit manually with GMM
g_centered <- function(alpha) {
  eta <- design_centered %*% alpha
  pi_vec <- inv_link(eta)
  rw <- r_vec / as.vector(pi_vec)
  g.matrix <- (rw - 1) * h_x
  return(g.matrix)
}

obj_centered <- function(alpha) {
  g.mat <- g_centered(alpha)
  G <- colMeans(g.mat)
  W_mat <- tryCatch(solve(crossprod(g.mat) / n_sub), error = function(e) diag(3))
  as.numeric(t(G) %*% W_mat %*% G)
}

# Find solution
opt_cent <- optim(c(-0.2, -1.5, -0.4), obj_centered, method = "L-BFGS-B",
                  lower = rep(-10, 3), upper = rep(10, 3), control = list(maxit = 2000))
alpha_cent <- opt_cent$par

cat("Point estimates:", round(alpha_cent, 4), "\n\n")

# Compute analytical SE manually for centered approach
# SE = sqrt(diag((Gamma' W Gamma)^{-1} / n))
g_mat_cent <- g_centered(alpha_cent)
W_cent <- solve(crossprod(g_mat_cent) / n_sub)

# Compute Gamma numerically
Gamma_cent <- jacobian(function(a) colMeans(g_centered(a)), alpha_cent)

cat("Jacobian Gamma (column 2 = d(E[g])/d(alpha_y)):\n")
print(round(Gamma_cent, 6))
cat("\n")

cat("Note: Column 2 now has larger values because derivative involves (y-y_bar),\n")
cat("which is non-zero for y=0 cases.\n\n")

# Analytical SE (simple sandwich)
GtW <- t(Gamma_cent) %*% W_cent
bread <- solve(GtW %*% Gamma_cent)
se_analytical_cent <- sqrt(diag(bread / n_sub))

cat("Analytical SE (simple):", round(se_analytical_cent, 4), "\n")

# Bootstrap SE
boot_cent <- matrix(NA, n_boot, 3)

for (b in 1:n_boot) {
  idx <- sample(1:n_sub, replace = TRUE)
  r_b <- r_vec[idx]; design_b <- design_centered[idx,]; h_x_b <- h_x[idx,]

  obj_b <- function(alpha) {
    eta <- design_b %*% alpha
    pi_vec <- inv_link(eta)
    if (any(pi_vec <= 0.001) || any(pi_vec >= 0.999)) return(1e10)
    rw <- r_b / as.vector(pi_vec)
    g.mat <- (rw - 1) * h_x_b
    G <- colMeans(g.mat)
    W_mat <- tryCatch(solve(crossprod(g.mat) / n_sub), error = function(e) diag(3))
    as.numeric(t(G) %*% W_mat %*% G)
  }

  tryCatch({
    opt <- optim(alpha_cent, obj_b, method = "L-BFGS-B",
                 lower = rep(-10, 3), upper = rep(10, 3), control = list(maxit = 1000))
    if (opt$convergence == 0 && all(abs(opt$par) < 9.9)) boot_cent[b,] <- opt$par
  }, error = function(e) {})
}

se_boot_cent <- apply(boot_cent[complete.cases(boot_cent),], 2, sd)

cat("Bootstrap SE:         ", round(se_boot_cent, 4), "\n\n")

cat("SE Comparison (Centered):\n")
cat(sprintf("  %-12s  Analytical  Bootstrap  Ratio\n", ""))
cat(sprintf("  %-12s  %9.4f  %9.4f  %5.3f\n", "Intercept", se_analytical_cent[1], se_boot_cent[1], se_analytical_cent[1]/se_boot_cent[1]))
cat(sprintf("  %-12s  %9.4f  %9.4f  %5.3f\n", "o(y-y_bar)", se_analytical_cent[2], se_boot_cent[2], se_analytical_cent[2]/se_boot_cent[2]))
cat(sprintf("  %-12s  %9.4f  %9.4f  %5.3f\n", "father", se_analytical_cent[3], se_boot_cent[3], se_analytical_cent[3]/se_boot_cent[3]))

#------------------------------------------------------------------------------#
# SUMMARY
#------------------------------------------------------------------------------#
cat("\n================================================================================\n")
cat("SUMMARY\n")
cat("================================================================================\n\n")

cat("                          CURRENT              CENTERED\n")
cat("                       (1, y, x)           (1, y-y_bar, x)\n")
cat("--------------------------------------------------------------------\n")
cat(sprintf("Coef on y:            %8.4f               %8.4f\n", alpha_hat[2], alpha_cent[2]))
cat(sprintf("Analytical SE:        %8.4f               %8.4f\n", se_analytical[2], se_analytical_cent[2]))
cat(sprintf("Bootstrap SE:         %8.4f               %8.4f\n", se_boot_curr[2], se_boot_cent[2]))
cat(sprintf("Ratio (Ana/Boot):     %8.3f               %8.3f\n",
            se_analytical[2]/se_boot_curr[2], se_analytical_cent[2]/se_boot_cent[2]))
cat("--------------------------------------------------------------------\n")

cat("\nKey insight:\n")
cat("- Current approach: Analytical SE is", round(se_analytical[2]/se_boot_curr[2] * 100), "% of bootstrap SE\n")
cat("- Centered approach: Analytical SE is", round(se_analytical_cent[2]/se_boot_cent[2] * 100), "% of bootstrap SE\n")

cat("\n================================================================================\n")
