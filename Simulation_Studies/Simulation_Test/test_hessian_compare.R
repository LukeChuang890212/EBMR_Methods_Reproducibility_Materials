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

for (rep_i in c(1, 2, 3, 5, 10, 15, 20)) {
  set.seed(12345 + rep_i)
  dat <- setting4.B1(2000)
  ebmr <- EBMRAlgorithmFast4[["new"]]("y", ps_spec_13, dat, W_func)

  for (j in 1:2) {
    ps_fit <- ebmr[["ps_fit.list"]][[j]]
    alpha <- ps_fit[["coefficients"]]
    design_matrix <- ps_fit[["design_matrix"]]
    h_x <- ps_fit[["h_x"]]
    r_vec <- dat[["r"]]
    n <- length(r_vec)
    esteq_dim <- ncol(h_x)
    param_dim <- length(alpha)

    # Reconstruct g, dg, d2g
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

    g_mat <- g_func(alpha)
    W_hat <- tryCatch(solve(t(g_mat) %*% g_mat / n), error = function(e) diag(esteq_dim))
    dg_arr <- dg_func(alpha)
    Gamma_hat <- matrix(0, esteq_dim, param_dim)
    for (jj in 1:param_dim) Gamma_hat[, jj] <- .colMeans(dg_arr[, jj, , drop = FALSE], n, esteq_dim)
    G_vec <- .colMeans(g_mat, n, esteq_dim)
    wg <- as.vector(W_hat %*% G_vec)

    # Full analytical: 2*Gamma'*W*Gamma + B + B'
    H_analytic <- 2 * crossprod(Gamma_hat, W_hat %*% Gamma_hat)

    # B + B' term
    GtW <- crossprod(Gamma_hat, W_hat)
    B_mat <- matrix(0, param_dim, param_dim)
    for (k in 1:param_dim) {
      dg_k <- dg_arr[, k, , drop = FALSE]
      dim(dg_k) <- c(n, esteq_dim)
      dOmega_k <- (crossprod(dg_k, g_mat) + crossprod(g_mat, dg_k)) / n
      B_mat[, k] <- -2 * as.vector(GtW %*% dOmega_k %*% wg)
    }
    H_analytic <- H_analytic + B_mat + t(B_mat)

    eig_analytic <- eigen(H_analytic, symmetric = TRUE, only.values = TRUE)$values

    # Numerical CUE hessian
    Q_full <- function(a) {
      g_m <- g_func(a)
      G_v <- matrix(.colMeans(g_m, n, ncol(g_m)), ncol(g_m), 1)
      W_q <- tryCatch(solve(t(g_m) %*% g_m / n), error = function(e) diag(ncol(g_m)))
      v <- as.numeric(crossprod(G_v, W_q %*% G_v))
      if (!is.finite(v)) 1e8 else v
    }
    H_numerical <- hessian(Q_full, alpha)
    eig_numerical <- eigen(H_numerical, symmetric = TRUE, only.values = TRUE)$values

    ratio_a <- min(eig_analytic) / max(abs(eig_analytic))
    ratio_n <- min(eig_numerical) / max(abs(eig_numerical))
    pd_a <- min(eig_analytic) > 2e-3 * max(abs(eig_analytic)) && min(eig_analytic) > 0
    pd_n <- min(eig_numerical) > 2e-3 * max(abs(eig_numerical)) && min(eig_numerical) > 0

    cat(sprintf("Rep %2d Model %d: ||a||=%.2f  analytic=[%.2e,%.2e] ratio=%.4f pd=%s | numerical=[%.2e,%.2e] ratio=%.4f pd=%s | %s\n",
        rep_i, j, sqrt(sum(alpha^2)),
        min(eig_analytic), max(eig_analytic), ratio_a, pd_a,
        min(eig_numerical), max(eig_numerical), ratio_n, pd_n,
        ifelse(pd_a == pd_n, "AGREE", "DISAGREE")))
  }
}
