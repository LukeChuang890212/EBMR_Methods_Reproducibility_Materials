setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
source("Simulation.r")
library(EBMRalgorithmFast4)
library(numDeriv)

ps_spec <- get_ps_spec("9-alt1")
ps_spec_13 <- list(
  formula.list = ps_spec[["formula.list"]][c(1, 3)],
  h_alpha.list = ps_spec[["h_alpha.list"]][c(1, 3)],
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]]
)

n_val <- 2000

for (rep_i in c(3, 2, 10)) {
  set.seed(12345 + rep_i)
  dat <- setting4.B1(n_val)

  r_vec <- as.vector(dat$r)
  n <- n_val

  ebmr <- EBMRAlgorithmFast4[["new"]]("y", ps_spec_13, dat,
    function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix)))
  ps_fit2 <- ebmr[["ps_fit.list"]][[2]]
  design_matrix <- ps_fit2[["design_matrix"]]
  h_x <- ps_fit2[["h_x"]]
  h_dim <- ncol(h_x)
  alpha_dim <- ncol(design_matrix)
  alpha_hat <- ps_fit2[["coefficients"]]

  eps_mix <- 0.05
  compute_pi <- function(eta) eps_mix + (1 - 2*eps_mix) / (1 + exp(eta))
  compute_f <- function(eta) 1 / (1 + exp(eta))

  g_func <- function(alpha) {
    eta <- as.vector(design_matrix %*% alpha)
    pi_hat <- compute_pi(eta)
    (r_vec / pi_hat - 1) * h_x
  }
  G_func <- function(alpha) colMeans(g_func(alpha))
  Q_func <- function(alpha) {
    g_mat <- g_func(alpha)
    G <- colMeans(g_mat)
    W <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(h_dim))
    as.numeric(t(G) %*% W %*% G)
  }

  cat(sprintf("\n============================================================\n"))
  cat(sprintf("=== Rep %d (||alpha||=%.2f) ===\n", rep_i, sqrt(sum(alpha_hat^2))))
  cat(sprintf("============================================================\n"))

  # --- At the GMM solution ---
  eta_hat <- as.vector(design_matrix %*% alpha_hat)
  pi_hat <- compute_pi(eta_hat)
  f_hat <- compute_f(eta_hat)

  # Gamma = E[dg/dalpha] = E[-r/pi^2 * dpi/dalpha * h_x]
  # dpi/dalpha = (1-2*eps) * f'(eta) * X = -(1-2*eps) * f*(1-f) * X  [logistic complement]
  fprime <- -f_hat * (1 - f_hat)  # f'(eta) for logistic complement
  dpi_dalpha <- (1 - 2*eps_mix) * fprime  # scalar factor per obs (without X)

  # Gamma[l, k] = mean(-r/pi^2 * dpi/dalpha * X[,k] * h_x[,l])
  common_factor <- -r_vec / (pi_hat^2) * dpi_dalpha  # n vector
  Gamma_mat <- matrix(0, h_dim, alpha_dim)
  for (k in 1:alpha_dim) {
    for (l in 1:h_dim) {
      Gamma_mat[l, k] <- mean(common_factor * design_matrix[, k] * h_x[, l])
    }
  }

  cat(sprintf("\n--- Gamma (dG/dalpha) ---\n"))
  cat(sprintf("  Frobenius norm: %.6f\n", sqrt(sum(Gamma_mat^2))))
  cat(sprintf("  Singular values: %s\n", paste(round(svd(Gamma_mat)$d, 6), collapse = ", ")))

  # Check f'(eta) magnitude
  cat(sprintf("\n--- f'(eta) at solution ---\n"))
  cat(sprintf("  |f'(eta)| quantiles: %.6f %.6f %.6f %.6f %.6f\n",
      quantile(abs(fprime), 0.01), quantile(abs(fprime), 0.10),
      quantile(abs(fprime), 0.50), quantile(abs(fprime), 0.90),
      quantile(abs(fprime), 0.99)))
  cat(sprintf("  |f'(eta)| < 1e-4: %d / %d\n", sum(abs(fprime) < 1e-4), n))
  cat(sprintf("  |f'(eta)| < 1e-2: %d / %d\n", sum(abs(fprime) < 1e-2), n))

  # G and W at solution
  g_mat <- g_func(alpha_hat)
  G <- colMeans(g_mat)
  W <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(h_dim))
  Q_val <- as.numeric(t(G) %*% W %*% G)

  cat(sprintf("\n--- FOC: 2*Gamma'*W*G ---\n"))
  foc <- 2 * as.vector(t(Gamma_mat) %*% W %*% G)
  cat(sprintf("  FOC: %s\n", paste(round(foc, 8), collapse = ", ")))
  cat(sprintf("  ||FOC||: %.8f\n", sqrt(sum(foc^2))))

  # Hessian of Q at solution (numerical)
  cat(sprintf("\n--- Hessian of Q(alpha) ---\n"))
  H <- hessian(Q_func, alpha_hat)
  eig_H <- eigen(H, only.values = TRUE)$values
  cat(sprintf("  Eigenvalues: %s\n", paste(round(eig_H, 6), collapse = ", ")))
  cat(sprintf("  Min eigenvalue: %.6f\n", min(eig_H)))
  cat(sprintf("  Condition number: %.2f\n", max(abs(eig_H)) / max(min(abs(eig_H)), 1e-15)))
  cat(sprintf("  Positive definite: %s\n", all(eig_H > 1e-8)))

  # Also check Hessian of Q with W=I (simpler: Hessian of G'G)
  Q_I_func <- function(alpha) {
    G <- G_func(alpha)
    sum(G^2)
  }
  H_I <- hessian(Q_I_func, alpha_hat)
  eig_H_I <- eigen(H_I, only.values = TRUE)$values
  cat(sprintf("\n--- Hessian of G'G (W=I) ---\n"))
  cat(sprintf("  Eigenvalues: %s\n", paste(round(eig_H_I, 6), collapse = ", ")))
  cat(sprintf("  Min eigenvalue: %.6f\n", min(eig_H_I)))
  cat(sprintf("  Positive definite: %s\n", all(eig_H_I > 1e-8)))

  cat(sprintf("\n  Q(alpha) = %.6f, ||G|| = %.6f\n", Q_val, sqrt(sum(G^2))))
}
