#------------------------------------------------------------------------------#
# Compare Centering Approach vs Current Approach for MNAR Model
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
cat("     COMPARE CENTERING APPROACH VS CURRENT APPROACH                            \n")
cat("================================================================================\n\n")

# Focus on Health = 1, Model 1 (o(y) + father)
health_val <- 1
subdat <- dat[dat$health == health_val, ]
n_sub <- nrow(subdat)

cat("Health =", health_val, ", n =", n_sub, "\n")
cat("Missing rate:", round(100 * mean(subdat$r == 0), 1), "%\n\n")

# Setup common elements
r_vec <- subdat$r
y_vec <- subdat$teacher_report
x_vec <- subdat$father
h_x <- cbind(1, subdat$father, subdat$parent_report)
inv_link <- function(eta) 1 / (1 + exp(eta))

# Compute y_bar (mean of y among observed)
y_bar <- mean(y_vec[r_vec == 1])
cat("Mean of y among observed (y_bar):", round(y_bar, 4), "\n\n")

#------------------------------------------------------------------------------#
# CURRENT APPROACH: design_matrix = (1, y, x)
#------------------------------------------------------------------------------#
cat("================================================================================\n")
cat("CURRENT APPROACH: design_matrix = (1, y, x)\n")
cat("================================================================================\n\n")

design_current <- cbind(1, y_vec, x_vec)

g_current <- function(alpha) {
  eta <- design_current %*% alpha
  pi_vec <- inv_link(eta)
  rw <- r_vec / as.vector(pi_vec)
  g.matrix <- (rw - 1) * h_x
  return(g.matrix)
}

obj_current <- function(alpha) {
  g.mat <- g_current(alpha)
  G <- colMeans(g.mat)
  W_mat <- tryCatch(solve(crossprod(g.mat) / n_sub), error = function(e) diag(3))
  val <- as.numeric(t(G) %*% W_mat %*% G)
  if (!is.finite(val)) return(1e10)
  return(val)
}

# Find solution
set.seed(12345)
best_obj_curr <- Inf
best_alpha_curr <- NULL

starts <- rbind(
  c(0, 0, 0), c(0, 1, 0), c(0, -1, 0), c(0, 2, -0.5), c(-0.5, 1.5, -0.3),
  c(0.1, -1.5, -0.4), c(-0.5, -2, 0.5)
)

for (i in 1:nrow(starts)) {
  opt <- optim(starts[i,], obj_current, method = "L-BFGS-B",
               lower = rep(-10, 3), upper = rep(10, 3),
               control = list(maxit = 2000))
  if (opt$value < best_obj_curr && opt$convergence == 0) {
    best_obj_curr <- opt$value
    best_alpha_curr <- opt$par
  }
}

cat("Solution: alpha =", round(best_alpha_curr, 4), "\n")
cat("Objective:", format(best_obj_curr, scientific = TRUE), "\n\n")

# Check derivative factor for each observation
pi_curr <- inv_link(design_current %*% best_alpha_curr)
deriv_factor_curr <- r_vec * (1 - pi_curr) / pi_curr * y_vec

cat("Derivative factor r*(1-pi)/pi*y:\n")
cat("  Among y=0 (observed): all =", unique(deriv_factor_curr[r_vec == 1 & y_vec == 0]), "\n")
cat("  Among y=1 (observed): mean =", round(mean(deriv_factor_curr[r_vec == 1 & y_vec == 1]), 4), "\n")
cat("  -> Only y=1 observations contribute to alpha_y!\n\n")

# Bootstrap SE
n_boot <- 500
boot_curr <- matrix(NA, n_boot, 3)

for (b in 1:n_boot) {
  idx <- sample(1:n_sub, replace = TRUE)
  r_b <- r_vec[idx]; design_b <- design_current[idx,]; h_x_b <- h_x[idx,]

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
    opt <- optim(best_alpha_curr, obj_b, method = "L-BFGS-B",
                 lower = rep(-10, 3), upper = rep(10, 3), control = list(maxit = 1000))
    if (opt$convergence == 0 && all(abs(opt$par) < 9.9)) boot_curr[b,] <- opt$par
  }, error = function(e) {})
}

valid_curr <- complete.cases(boot_curr)
se_boot_curr <- apply(boot_curr[valid_curr,], 2, sd)

cat("Bootstrap SE (n_valid =", sum(valid_curr), "):\n")
cat("  SE:", round(se_boot_curr, 4), "\n\n")

z_curr <- best_alpha_curr / se_boot_curr
p_curr <- 2 * pnorm(-abs(z_curr))

cat("Results:\n")
cat(sprintf("  Intercept: %.4f (SE: %.4f, z: %.3f, p: %.4f)\n",
            best_alpha_curr[1], se_boot_curr[1], z_curr[1], p_curr[1]))
cat(sprintf("  o(y):      %.4f (SE: %.4f, z: %.3f, p: %.4f)\n",
            best_alpha_curr[2], se_boot_curr[2], z_curr[2], p_curr[2]))
cat(sprintf("  father:    %.4f (SE: %.4f, z: %.3f, p: %.4f)\n",
            best_alpha_curr[3], se_boot_curr[3], z_curr[3], p_curr[3]))

#------------------------------------------------------------------------------#
# CENTERING APPROACH: design_matrix = (1, y - y_bar, x)
#------------------------------------------------------------------------------#
cat("\n================================================================================\n")
cat("CENTERING APPROACH: design_matrix = (1, y - y_bar, x)\n")
cat("================================================================================\n\n")

y_centered <- y_vec - y_bar
design_centered <- cbind(1, y_centered, x_vec)

cat("y_bar =", round(y_bar, 4), "\n")
cat("Centered y values among observed: y=0 ->", round(0 - y_bar, 4),
    ", y=1 ->", round(1 - y_bar, 4), "\n\n")

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
  val <- as.numeric(t(G) %*% W_mat %*% G)
  if (!is.finite(val)) return(1e10)
  return(val)
}

# Find solution
best_obj_cent <- Inf
best_alpha_cent <- NULL

for (i in 1:nrow(starts)) {
  opt <- optim(starts[i,], obj_centered, method = "L-BFGS-B",
               lower = rep(-10, 3), upper = rep(10, 3),
               control = list(maxit = 2000))
  if (opt$value < best_obj_cent && opt$convergence == 0) {
    best_obj_cent <- opt$value
    best_alpha_cent <- opt$par
  }
}

cat("Solution: alpha =", round(best_alpha_cent, 4), "\n")
cat("Objective:", format(best_obj_cent, scientific = TRUE), "\n\n")

# Check derivative factor for each observation
pi_cent <- inv_link(design_centered %*% best_alpha_cent)
deriv_factor_cent <- r_vec * (1 - pi_cent) / pi_cent * y_centered

cat("Derivative factor r*(1-pi)/pi*(y-y_bar):\n")
cat("  Among y=0 (observed): mean =", round(mean(deriv_factor_cent[r_vec == 1 & y_vec == 0]), 4), "\n")
cat("  Among y=1 (observed): mean =", round(mean(deriv_factor_cent[r_vec == 1 & y_vec == 1]), 4), "\n")
cat("  -> BOTH y=0 and y=1 observations contribute to alpha_y!\n\n")

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
    opt <- optim(best_alpha_cent, obj_b, method = "L-BFGS-B",
                 lower = rep(-10, 3), upper = rep(10, 3), control = list(maxit = 1000))
    if (opt$convergence == 0 && all(abs(opt$par) < 9.9)) boot_cent[b,] <- opt$par
  }, error = function(e) {})
}

valid_cent <- complete.cases(boot_cent)
se_boot_cent <- apply(boot_cent[valid_cent,], 2, sd)

cat("Bootstrap SE (n_valid =", sum(valid_cent), "):\n")
cat("  SE:", round(se_boot_cent, 4), "\n\n")

z_cent <- best_alpha_cent / se_boot_cent
p_cent <- 2 * pnorm(-abs(z_cent))

cat("Results:\n")
cat(sprintf("  Intercept:   %.4f (SE: %.4f, z: %.3f, p: %.4f)\n",
            best_alpha_cent[1], se_boot_cent[1], z_cent[1], p_cent[1]))
cat(sprintf("  o(y-y_bar):  %.4f (SE: %.4f, z: %.3f, p: %.4f)\n",
            best_alpha_cent[2], se_boot_cent[2], z_cent[2], p_cent[2]))
cat(sprintf("  father:      %.4f (SE: %.4f, z: %.3f, p: %.4f)\n",
            best_alpha_cent[3], se_boot_cent[3], z_cent[3], p_cent[3]))

#------------------------------------------------------------------------------#
# COMPARISON
#------------------------------------------------------------------------------#
cat("\n================================================================================\n")
cat("COMPARISON\n")
cat("================================================================================\n\n")

cat("                        CURRENT          CENTERED\n")
cat("                        (1, y, x)        (1, y-y_bar, x)\n")
cat("----------------------------------------------------------------\n")
cat(sprintf("Intercept:              %7.4f          %7.4f\n", best_alpha_curr[1], best_alpha_cent[1]))
cat(sprintf("Coef on y:              %7.4f          %7.4f\n", best_alpha_curr[2], best_alpha_cent[2]))
cat(sprintf("Coef on father:         %7.4f          %7.4f\n", best_alpha_curr[3], best_alpha_cent[3]))
cat("----------------------------------------------------------------\n")
cat(sprintf("SE(Intercept):          %7.4f          %7.4f\n", se_boot_curr[1], se_boot_cent[1]))
cat(sprintf("SE(y coef):             %7.4f          %7.4f\n", se_boot_curr[2], se_boot_cent[2]))
cat(sprintf("SE(father):             %7.4f          %7.4f\n", se_boot_curr[3], se_boot_cent[3]))
cat("----------------------------------------------------------------\n")
cat(sprintf("p-value(y coef):        %7.4f          %7.4f\n", p_curr[2], p_cent[2]))
cat("----------------------------------------------------------------\n")

# Check if propensity scores are the same
cat("\nPropensity score comparison:\n")
cat("  Correlation:", round(cor(as.vector(pi_curr), as.vector(pi_cent)), 6), "\n")
cat("  Max abs diff:", round(max(abs(pi_curr - pi_cent)), 6), "\n")

# The coefficients should give equivalent models (just reparameterized)
# alpha_0_new = alpha_0_old - alpha_1 * y_bar
# alpha_1_new = alpha_1_old
cat("\nReparameterization check:\n")
cat("  alpha_1 (current):", round(best_alpha_curr[2], 4), "\n")
cat("  alpha_1 (centered):", round(best_alpha_cent[2], 4), "\n")
cat("  Expected alpha_0 (centered) = alpha_0 - alpha_1 * y_bar:\n")
cat("    ", round(best_alpha_curr[1] - best_alpha_curr[2] * y_bar, 4), "\n")
cat("  Actual alpha_0 (centered):", round(best_alpha_cent[1], 4), "\n")

cat("\n================================================================================\n")
