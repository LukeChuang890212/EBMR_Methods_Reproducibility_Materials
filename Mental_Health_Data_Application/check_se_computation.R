#------------------------------------------------------------------------------#
# Check Standard Error Computation in GMM
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

cat("================================================================================\n")
cat("                    CHECK STANDARD ERROR COMPUTATION                            \n")
cat("================================================================================\n\n")

#------------------------------------------------------------------------------#
# Setup
#------------------------------------------------------------------------------#
full_ps_specifications <- list(
  formula.list = list(
    r ~ o(teacher_report) + father,
    r ~ o(teacher_report) + parent_report,
    r ~ father + parent_report
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
    c(0.5, -0.2, -0.5)
  ),
  alpha_init.list1 = list(
    c(0.09019492, -1.511548, -0.3946616),
    c(-1.969056, 3.5823611, -1.37232734),
    c(0.5, -0.2, -0.5)
  )
)

#------------------------------------------------------------------------------#
# Focus on Health = 1, Model 1 for detailed check
#------------------------------------------------------------------------------#
health_val <- 1
subdat <- dat[dat$health == health_val, ]
n_sub <- nrow(subdat)

cat("Checking Health =", health_val, ", Model 1 (o(y) + father)\n")
cat("Sample size: n =", n_sub, "\n\n")

ps_specifications <- list(
  formula.list = full_ps_specifications$formula.list[1],
  h_x_names.list = full_ps_specifications$h_x_names.list[1],
  alpha_init.list = list(init.list[[health_val + 1]][[1]]),
  inv_link = full_ps_specifications$inv_link
)

ebmr <- EBMRAlgorithmFast$new("teacher_report", ps_specifications, subdat, W)
ps_fit <- ebmr$ps_fit.list[[1]]
gmm_fit <- ps_fit$gmm_fit

alpha_hat <- ps_fit$coefficients
se_reported <- ps_fit$se

cat("--------------------------------------------------------------------------------\n")
cat("REPORTED ESTIMATES AND SE\n")
cat("--------------------------------------------------------------------------------\n")
cat("alpha.hat:", round(alpha_hat, 6), "\n")
cat("SE:       ", round(se_reported, 6), "\n")
cat("z-values: ", round(alpha_hat / se_reported, 3), "\n\n")

#------------------------------------------------------------------------------#
# Manually reconstruct the SE computation
#------------------------------------------------------------------------------#
cat("--------------------------------------------------------------------------------\n")
cat("MANUAL SE COMPUTATION CHECK\n")
cat("--------------------------------------------------------------------------------\n\n")

# Get components from gmm_fit
Gamma.hat <- gmm_fit$Gamma.hat
W.hat <- gmm_fit$W.hat
g.matrix <- gmm_fit$g.matrix
eta_s <- gmm_fit$eta_s
Q <- gmm_fit$Q
psi <- gmm_fit$psi

cat("Dimensions:\n")
cat("  Gamma.hat:", dim(Gamma.hat), "\n")
cat("  W.hat:", dim(W.hat), "\n")
cat("  g.matrix:", dim(g.matrix), "\n")
cat("  eta_s:", length(eta_s), "\n")
cat("  Q:", dim(Q), "\n")
cat("  psi:", dim(psi), "\n\n")

# Check eta_s (should be close to 0 at solution)
cat("eta_s (moment conditions at solution, should be ~0):\n")
cat("  ", round(eta_s, 8), "\n\n")

# Check psi (influence function)
cat("psi summary (influence function, param_dim x n):\n")
cat("  Mean of psi rows:", round(rowMeans(psi), 8), "\n")
cat("  Var of psi rows: ", round(apply(psi, 1, var), 6), "\n\n")

# Recompute covariance
cov_manual <- var(t(psi)) / n_sub
se_manual <- sqrt(diag(cov_manual))

cat("SE comparison:\n")
cat("  Reported SE: ", round(se_reported, 6), "\n")
cat("  Manual SE:   ", round(se_manual, 6), "\n")
cat("  Difference:  ", round(se_reported - se_manual, 10), "\n\n")

#------------------------------------------------------------------------------#
# Check the sandwich formula components
#------------------------------------------------------------------------------#
cat("--------------------------------------------------------------------------------\n")
cat("SANDWICH FORMULA COMPONENTS\n")
cat("--------------------------------------------------------------------------------\n\n")

# Standard GMM SE formula (simplified, when eta_s = 0):
# Var(alpha.hat) = (Gamma' W Gamma)^{-1} Gamma' W Var(g) W Gamma (Gamma' W Gamma)^{-1} / n
# where Var(g) = E[g g'] - E[g] E[g]' ≈ E[g g'] when E[g] = 0

# Bread matrix: (Gamma' W Gamma)^{-1}
GtW <- t(Gamma.hat) %*% W.hat
bread <- solve(GtW %*% Gamma.hat)

cat("Bread matrix (Gamma' W Gamma)^{-1}:\n")
print(round(bread, 6))
cat("\n")

# Meat matrix: Gamma' W Var(g) W Gamma
Var_g <- var(g.matrix)  # This is the sample variance
meat <- GtW %*% Var_g %*% W.hat %*% Gamma.hat

cat("Meat matrix (Gamma' W Var(g) W Gamma):\n")
print(round(meat, 6))
cat("\n")

# Sandwich covariance
cov_sandwich <- bread %*% meat %*% bread / n_sub
se_sandwich <- sqrt(diag(cov_sandwich))

cat("SE from sandwich formula:\n")
cat("  ", round(se_sandwich, 6), "\n\n")

#------------------------------------------------------------------------------#
# Alternative: Simple asymptotic variance (when W is optimal)
#------------------------------------------------------------------------------#
cat("--------------------------------------------------------------------------------\n")
cat("ALTERNATIVE SE FORMULAS\n")
cat("--------------------------------------------------------------------------------\n\n")

# When W = (E[g g'])^{-1} is the optimal weight, the asymptotic variance simplifies to:
# Var(alpha.hat) = (Gamma' W Gamma)^{-1} / n

cov_simple <- bread / n_sub
se_simple <- sqrt(diag(cov_simple))

cat("SE from simple formula (assuming optimal W):\n")
cat("  ", round(se_simple, 6), "\n\n")

# Check if W.hat is close to optimal
W_optimal <- solve(t(g.matrix) %*% g.matrix / n_sub)
cat("Is W.hat close to optimal W?\n")
cat("  Max abs diff:", max(abs(W.hat - W_optimal)), "\n\n")

#------------------------------------------------------------------------------#
# Bootstrap check
#------------------------------------------------------------------------------#
cat("--------------------------------------------------------------------------------\n")
cat("BOOTSTRAP SE CHECK\n")
cat("--------------------------------------------------------------------------------\n\n")

set.seed(12345)
n_boot <- 200

# Create the moment function for this model
r_vec <- subdat$r
y_vec <- subdat$teacher_report
x_vec <- subdat$father
h_x <- cbind(1, subdat$father, subdat$parent_report)
design_matrix <- cbind(1, y_vec, x_vec)
inv_link <- function(eta) 1 / (1 + exp(eta))

Phi_alpha <- function(param) {
  eta <- design_matrix %*% param
  pi_vec <- inv_link(eta)
  rw <- r_vec / as.vector(pi_vec)
  g.matrix <- (rw - 1) * h_x
  return(g.matrix)
}

obj_func <- function(param) {
  g.mat <- Phi_alpha(param)
  G <- colMeans(g.mat)
  W_mat <- tryCatch(solve(crossprod(g.mat) / n_sub), error = function(e) diag(3))
  return(as.numeric(t(G) %*% W_mat %*% G))
}

boot_estimates <- matrix(NA, n_boot, 3)

cat("Running", n_boot, "bootstrap replications...\n")

for (b in 1:n_boot) {
  idx <- sample(1:n_sub, replace = TRUE)
  boot_dat <- subdat[idx, ]

  # Recreate for bootstrap sample
  r_boot <- boot_dat$r
  y_boot <- boot_dat$teacher_report
  x_boot <- boot_dat$father
  h_x_boot <- cbind(1, boot_dat$father, boot_dat$parent_report)
  design_boot <- cbind(1, y_boot, x_boot)

  Phi_boot <- function(param) {
    eta <- design_boot %*% param
    pi_vec <- inv_link(eta)
    rw <- r_boot / as.vector(pi_vec)
    g.matrix <- (rw - 1) * h_x_boot
    return(g.matrix)
  }

  obj_boot <- function(param) {
    g.mat <- Phi_boot(param)
    G <- colMeans(g.mat)
    W_mat <- tryCatch(solve(crossprod(g.mat) / n_sub), error = function(e) diag(3))
    return(as.numeric(t(G) %*% W_mat %*% G))
  }

  tryCatch({
    opt <- optim(alpha_hat, obj_boot, method = "L-BFGS-B",
                 lower = rep(-10, 3), upper = rep(10, 3),
                 control = list(maxit = 1000))
    boot_estimates[b, ] <- opt$par
  }, error = function(e) {
    boot_estimates[b, ] <- NA
  })
}

# Remove failed bootstrap samples
boot_estimates <- boot_estimates[complete.cases(boot_estimates), ]
cat("Successful bootstrap samples:", nrow(boot_estimates), "\n\n")

se_boot <- apply(boot_estimates, 2, sd)

cat("SE comparison:\n")
cat("  Reported SE:  ", round(se_reported, 6), "\n")
cat("  Bootstrap SE: ", round(se_boot, 6), "\n")
cat("  Ratio (Rep/Boot):", round(se_reported / se_boot, 3), "\n\n")

#------------------------------------------------------------------------------#
# Check the Q matrix and M_n term
#------------------------------------------------------------------------------#
cat("--------------------------------------------------------------------------------\n")
cat("CHECK Q MATRIX AND M_n TERM\n")
cat("--------------------------------------------------------------------------------\n\n")

# Q = (I - H_0n M_n)^{-1} H_0n
# H_0n = -(Gamma' W Gamma)^{-1}

H_0n <- -bread
cat("H_0n (should be -bread):\n")
print(round(H_0n, 6))
cat("\n")

cat("Q matrix:\n")
print(round(Q, 6))
cat("\n")

# Check if Q ≈ H_0n (which happens when M_n ≈ 0)
cat("Q - H_0n (should be ~0 if M_n is small):\n")
print(round(Q - H_0n, 8))
cat("\n")

# The M_n term involves eta_s, which should be 0 at the solution
# So M_n should be 0
cat("eta_s (should be 0):", round(eta_s, 10), "\n")
cat("If eta_s = 0, then M_n = 0 and Q = H_0n\n\n")

#------------------------------------------------------------------------------#
# Summary
#------------------------------------------------------------------------------#
cat("================================================================================\n")
cat("SUMMARY\n")
cat("================================================================================\n\n")

cat("SE from different methods:\n")
cat(sprintf("  %-25s: %.6f, %.6f, %.6f\n", "Reported (gmm_fit)", se_reported[1], se_reported[2], se_reported[3]))
cat(sprintf("  %-25s: %.6f, %.6f, %.6f\n", "Manual (var(psi)/n)", se_manual[1], se_manual[2], se_manual[3]))
cat(sprintf("  %-25s: %.6f, %.6f, %.6f\n", "Sandwich formula", se_sandwich[1], se_sandwich[2], se_sandwich[3]))
cat(sprintf("  %-25s: %.6f, %.6f, %.6f\n", "Simple (bread/n)", se_simple[1], se_simple[2], se_simple[3]))
cat(sprintf("  %-25s: %.6f, %.6f, %.6f\n", "Bootstrap", se_boot[1], se_boot[2], se_boot[3]))

cat("\n")
cat("================================================================================\n")
