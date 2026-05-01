setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
library(EBMRalgorithmFast4)
library(numDeriv)

ps_spec <- get_ps_spec("9-alt1")
ps_spec_2 <- list(
  formula.list = ps_spec[["formula.list"]][2],
  h_alpha.list = ps_spec[["h_alpha.list"]][2],
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]]
)
W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
n <- 2000

all_data <- readRDS("Simulation_Data/Setting3.B1_n2000_replicate1000.RDS")

# Get a single rep's fitted model to extract alpha, g, etc.
dat <- all_data[1:n, ]
ebmr <- EBMRAlgorithmFast4$new("y", ps_spec_2, dat, W_func)

ps_fit <- ebmr$ps_fit.list[[1]]
alpha <- ps_fit$coefficients
design_matrix <- ps_fit$design_matrix
h_x <- ps_fit$h_x
r_vec <- dat[["r"]]
esteq_dim <- ncol(h_x)
param_dim <- length(alpha)

# Reconstruct g and dg
g_func <- function(a) {
  eta <- design_matrix %*% a
  eta_c <- pmin(pmax(eta, -20), 20)
  pi_val <- as.vector(1 / (1 + exp(eta_c)))
  (r_vec / pi_val - 1) * h_x
}

dg_func <- function(a) {
  eta <- as.vector(design_matrix %*% a)
  pi_vec <- as.vector(1 / (1 + exp(pmin(pmax(eta, -20), 20))))
  cf <- r_vec * (1 - pi_vec) / pi_vec
  arr <- array(0, dim = c(n, param_dim, esteq_dim))
  for (l in 1:esteq_dim) arr[, , l] <- (cf * h_x[, l]) * design_matrix
  arr
}

cat("=== Timing individual components ===\n\n")

# 1. Time the full numerical CUE hessian (numDeriv::hessian)
Q_full <- function(a) {
  g_m <- g_func(a)
  G_v <- matrix(.colMeans(g_m, n, ncol(g_m)), ncol(g_m), 1)
  W_q <- tryCatch(solve(t(g_m) %*% g_m / n), error = function(e) diag(ncol(g_m)))
  v <- as.numeric(crossprod(G_v, W_q %*% G_v))
  if (!is.finite(v)) 1e8 else v
}

t0 <- Sys.time()
for (i in 1:5) H_num <- hessian(Q_full, alpha)
t_numerical <- as.numeric(difftime(Sys.time(), t0, units = "secs")) / 5
cat(sprintf("numDeriv::hessian(Q_full):  %.3f sec\n", t_numerical))

# 2. Time the analytical hessian (all 3 terms)
t0 <- Sys.time()
for (i in 1:5) {
  g_mat <- g_func(alpha)
  W_h <- tryCatch(solve(t(g_mat) %*% g_mat / n), error = function(e) diag(esteq_dim))
  dg_arr <- dg_func(alpha)
  Gamma_h <- matrix(0, esteq_dim, param_dim)
  for (j in 1:param_dim) Gamma_h[, j] <- .colMeans(dg_arr[, j, , drop = FALSE], n, esteq_dim)
  G_h <- .colMeans(g_mat, n, esteq_dim)
  wg <- as.vector(W_h %*% G_h)

  H_Q <- 2 * crossprod(Gamma_h, W_h %*% Gamma_h)

  GtW <- crossprod(Gamma_h, W_h)
  B_mat <- matrix(0, param_dim, param_dim)
  for (k in 1:param_dim) {
    dg_k <- dg_arr[, k, , drop = FALSE]
    dim(dg_k) <- c(n, esteq_dim)
    dOmega_k <- (crossprod(dg_k, g_mat) + crossprod(g_mat, dg_k)) / n
    B_mat[, k] <- -2 * as.vector(GtW %*% dOmega_k %*% wg)
  }
  H_Q <- H_Q + B_mat + t(B_mat)
}
t_analytical <- as.numeric(difftime(Sys.time(), t0, units = "secs")) / 5
cat(sprintf("Analytical hessian (3 terms): %.3f sec\n", t_analytical))

# 3. Time the GMM optimization alone (no hessian)
t0 <- Sys.time()
for (i in 1:5) {
  dat_i <- all_data[((i-1)*n + 1):(i*n), ]
  ebmr_i <- EBMRAlgorithmFast4$new("y", ps_spec_2, dat_i, W_func)
}
t_total <- as.numeric(difftime(Sys.time(), t0, units = "secs")) / 5
cat(sprintf("Full EBMRAlgorithmFast4$new: %.3f sec\n", t_total))

# 4. Time conv_grad evaluation
conv_grad <- function(param) {
  g_mat <- g_func(param)
  G_vec <- matrix(.colMeans(g_mat, n, esteq_dim), esteq_dim, 1)
  W_hat <- tryCatch(solve(t(g_mat) %*% g_mat / n), error = function(e) diag(esteq_dim))
  dg_arr <- dg_func(param)
  Gamma_mat <- matrix(0, esteq_dim, param_dim)
  for (j in 1:param_dim) Gamma_mat[, j] <- .colMeans(dg_arr[, j, , drop = FALSE], n, esteq_dim)
  2 * as.vector(crossprod(Gamma_mat, W_hat %*% G_vec))
}

t0 <- Sys.time()
for (i in 1:5) H_jac <- jacobian(conv_grad, alpha)
t_jacobian <- as.numeric(difftime(Sys.time(), t0, units = "secs")) / 5
cat(sprintf("numDeriv::jacobian(conv_grad): %.3f sec\n", t_jacobian))

cat(sprintf("\nSpeedup: numerical/analytical = %.1fx\n", t_numerical / t_analytical))
cat(sprintf("Speedup: jacobian/analytical = %.1fx\n", t_jacobian / t_analytical))
