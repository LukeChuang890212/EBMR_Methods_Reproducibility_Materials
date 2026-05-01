setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
library(EBMRalgorithmFast4)

ps_spec <- get_ps_spec("9-alt1")
ps_spec_13 <- list(
  formula.list = ps_spec[["formula.list"]][c(1, 3)],
  h_alpha.list = ps_spec[["h_alpha.list"]][c(1, 3)],
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]]
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
n_val <- 2000

set.seed(12345 + 2)
dat <- setting4.B1(n_val)

ebmr <- EBMRAlgorithmFast4[["new"]]("y", ps_spec_13, dat, W_func)
ps_fit2 <- ebmr[["ps_fit.list"]][[2]]
alpha <- ps_fit2[["coefficients"]]
pi_hat <- ps_fit2[["fitted.values"]]
cat("Alpha:", round(alpha, 4), "\n")
cat("||alpha||:", round(sqrt(sum(alpha^2)), 3), "\n")
cat("Pi quantiles:", round(quantile(pi_hat, c(0.01, 0.05, 0.10, 0.50, 0.90, 0.95, 0.99)), 4), "\n")
cat("Pi in [0.06, 0.94]:", sum(pi_hat >= 0.06 & pi_hat <= 0.94), "/", n_val, "\n")

# Check Hessian eigenvalues
library(numDeriv)
design_matrix <- ps_fit2[["design_matrix"]]
h_x <- ps_fit2[["h_x"]]
r_vec <- as.vector(dat$r)
eps_mix <- 0.05
compute_pi <- function(eta) eps_mix + (1 - 2*eps_mix) / (1 + exp(eta))

g_func <- function(a) {
  eta <- as.vector(design_matrix %*% a)
  pi_hat <- compute_pi(eta)
  (r_vec / pi_hat - 1) * h_x
}
Q_func <- function(a) {
  g_mat <- g_func(a)
  G <- colMeans(g_mat)
  W <- tryCatch(solve(crossprod(g_mat) / n_val), error = function(e) diag(ncol(h_x)))
  as.numeric(t(G) %*% W %*% G)
}

H <- hessian(Q_func, alpha)
eig_H <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
cat("\nHessian eigenvalues:", round(eig_H, 6), "\n")
cat("Min/Max ratio:", round(min(eig_H) / max(abs(eig_H)), 6), "\n")
cat("Min eigenvalue:", round(min(eig_H), 6), "\n")
