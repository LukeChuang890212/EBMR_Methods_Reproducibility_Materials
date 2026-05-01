setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
source("Simulation.r")
library(EBMRalgorithmFast4)

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

  # Get design matrix and h_x from a quick package fit
  ebmr <- EBMRAlgorithmFast4[["new"]]("y", ps_spec_13, dat,
    function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix)))
  ps_fit2 <- ebmr[["ps_fit.list"]][[2]]
  design_matrix <- ps_fit2[["design_matrix"]]
  h_x <- ps_fit2[["h_x"]]
  h_dim <- ncol(h_x)
  alpha_dim <- ncol(design_matrix)

  cat(sprintf("\n============================================================\n"))
  cat(sprintf("=== Rep %d (package ||alpha||=%.2f) ===\n", rep_i,
      sqrt(sum(ps_fit2[["coefficients"]]^2))))
  cat(sprintf("============================================================\n"))

  # Mixture model
  eps_mix <- 0.05
  compute_pi <- function(eta) {
    eps_mix + (1 - 2*eps_mix) / (1 + exp(eta))
  }

  g_func <- function(alpha) {
    eta <- as.vector(design_matrix %*% alpha)
    pi_hat <- compute_pi(eta)
    (r_vec / pi_hat - 1) * h_x
  }

  G_func <- function(alpha) colMeans(g_func(alpha))

  W_func_local <- function(g_mat) {
    tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(h_dim))
  }

  # Step 0: W = I, optimize
  alpha <- rep(0, alpha_dim)
  W <- diag(h_dim)

  opt0 <- optim(alpha, function(a) {
    G <- G_func(a); as.numeric(t(G) %*% W %*% G)
  }, method = "L-BFGS-B", control = list(maxit = 1000))
  alpha <- opt0$par

  cat(sprintf("\n%4s %10s %12s %12s %10s %8s %s\n",
      "iter", "||alpha||", "G'WG", "G'G", "max_eig_W", "pi_bnd", "alpha"))

  # Print step 0
  g_mat <- g_func(alpha)
  G <- G_func(alpha)
  GG <- sum(G^2)
  pi_hat <- compute_pi(as.vector(design_matrix %*% alpha))
  n_bnd <- sum(pi_hat < 0.06 | pi_hat > 0.94)
  cat(sprintf("%4d %10.4f %12.6f %12.6f %10.4f %8d   %s\n",
      0, sqrt(sum(alpha^2)), opt0$value, GG, 1.0, n_bnd,
      paste(round(alpha, 3), collapse=", ")))

  # Iterative GMM
  for (iter in 1:30) {
    # Update W
    g_mat <- g_func(alpha)
    W <- W_func_local(g_mat)
    eig_max <- max(eigen(W, only.values = TRUE)$values)

    # Optimize with new W
    opt <- optim(alpha, function(a) {
      G <- G_func(a); as.numeric(t(G) %*% W %*% G)
    }, method = "L-BFGS-B", control = list(maxit = 1000))
    alpha <- opt$par

    # Diagnostics
    G <- G_func(alpha)
    GWG <- opt$value
    GG <- sum(G^2)
    pi_hat <- compute_pi(as.vector(design_matrix %*% alpha))
    n_bnd <- sum(pi_hat < 0.06 | pi_hat > 0.94)

    cat(sprintf("%4d %10.4f %12.6f %12.6f %10.4f %8d   %s\n",
        iter, sqrt(sum(alpha^2)), GWG, GG, eig_max, n_bnd,
        paste(round(alpha, 3), collapse=", ")))

    if (GWG < 1e-10) break
  }
}
