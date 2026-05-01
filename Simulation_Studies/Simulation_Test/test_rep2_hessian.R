setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
library(EBMRalgorithmFast4)
library(numDeriv)

ps_spec <- get_ps_spec("9-alt1")
ps_spec_13 <- list(
  formula.list = ps_spec[["formula.list"]][c(1, 3)],
  h_alpha.list = ps_spec[["h_alpha.list"]][c(1, 3)],
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]]
)
W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

set.seed(12345 + 2)
dat <- setting4.B1(2000)
ebmr <- EBMRAlgorithmFast4[["new"]]("y", ps_spec_13, dat, W_func)

for (j in 1:2) {
  ps_fit <- ebmr[["ps_fit.list"]][[j]]
  alpha <- ps_fit[["coefficients"]]
  gmm_fit <- ps_fit[["gmm_fit"]]
  opt <- gmm_fit[["opt"]]

  cat(sprintf("\n=== Model %d ===\n", j))
  cat("Alpha:", round(alpha, 4), "\n")
  cat("Objective:", formatC(opt[["objective"]], format = "e", digits = 6), "\n")
  cat("hessian_pd:", opt[["hessian_pd"]], "\n")

  # Reconstruct g function from the stored ps_fit components
  design_matrix <- ps_fit[["design_matrix"]]
  h_x <- ps_fit[["h_x"]]
  r_vec <- dat[["r"]]
  n <- length(r_vec)
  esteq_dim <- ncol(h_x)

  g_func <- function(a) {
    eta <- design_matrix %*% a
    eta_c <- pmin(pmax(eta, -20), 20)
    pi_val <- as.vector(1 / (1 + exp(eta_c)))
    (r_vec / pi_val - 1) * h_x
  }

  Q_full <- function(a) {
    g_m <- g_func(a)
    G_vec <- matrix(.colMeans(g_m, n, esteq_dim), esteq_dim, 1)
    W_hat <- tryCatch(solve(t(g_m) %*% g_m / n), error = function(e) diag(esteq_dim))
    v <- as.numeric(crossprod(G_vec, W_hat %*% G_vec))
    if (!is.finite(v)) 1e8 else v
  }

  H_Q <- hessian(Q_full, alpha)
  eig_vals <- eigen(H_Q, symmetric = TRUE, only.values = TRUE)$values

  cat("Hessian eigenvalues:\n")
  for (k in seq_along(eig_vals)) {
    cat(sprintf("  eig[%d] = %s\n", k, formatC(eig_vals[k], format = "e", digits = 6)))
  }
  cat("min(eig):", formatC(min(eig_vals), format = "e", digits = 6), "\n")
  cat("max(|eig|):", formatC(max(abs(eig_vals)), format = "e", digits = 6), "\n")
  cat("Ratio min/max(|eig|):", formatC(min(eig_vals) / max(abs(eig_vals)), format = "e", digits = 6), "\n")
  cat("Threshold (1e-3 * max):", formatC(1e-3 * max(abs(eig_vals)), format = "e", digits = 6), "\n")
  cat("Passes check:", min(eig_vals) > 1e-3 * max(abs(eig_vals)) && min(eig_vals) > 0, "\n")
}
