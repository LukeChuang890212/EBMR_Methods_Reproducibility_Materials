setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
library(EBMRalgorithmFast4)

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

# After fitting, extract the fitted alpha and reconstruct the problem
# to trace what GMM convergence looks like

for (rep_i in c(2, 4, 7, 10)) {
  dat <- all_data[((rep_i-1)*n + 1):(rep_i*n), ]
  ebmr <- EBMRAlgorithmFast4$new("y", ps_spec_2, dat, W_func)
  ps_fit <- ebmr$ps_fit.list[[1]]
  alpha_final <- ps_fit$coefficients
  opt <- ps_fit$gmm_fit$opt

  cat(sprintf("\n=== Rep %d ===\n", rep_i))
  cat(sprintf("iter=%d conv=%s hpd=%s obj=%.4e grad=%.2e ||a||=%.3f\n",
      opt$iterations, opt$converged, opt$hessian_pd,
      opt$objective, opt$final_grad_norm, sqrt(sum(alpha_final^2))))
  cat(sprintf("alpha: %s\n", paste(round(alpha_final, 4), collapse=", ")))

  # Now manually trace the CUE loop using the same components
  design_matrix <- ps_fit$design_matrix
  h_x <- ps_fit$h_x
  r_vec <- dat[["r"]]
  esteq_dim <- ncol(h_x)
  param_dim <- ncol(design_matrix)

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

  cue_obj_and_grad <- function(param) {
    g_mat <- g_func(param)
    G_vec <- .colMeans(g_mat, n, esteq_dim)
    W_hat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(esteq_dim))
    wG <- as.vector(W_hat %*% G_vec)
    obj <- sum(G_vec * wG)
    dg_arr <- dg_func(param)
    Gamma_mat <- matrix(0, esteq_dim, param_dim)
    for (j in 1:param_dim) Gamma_mat[, j] <- .colMeans(dg_arr[, j, , drop = FALSE], n, esteq_dim)
    grad <- 2 * as.vector(crossprod(Gamma_mat, wG))
    list(obj = obj, grad = grad, grad_norm = max(abs(grad)))
  }

  # Trace from init=0
  estimates <- rep(0, param_dim)
  W_hat_fixed <- diag(esteq_dim)

  obj_fixed <- function(p) {
    g_m <- g_func(p)
    G_v <- matrix(.colMeans(g_m, n, esteq_dim), esteq_dim, 1)
    v <- as.numeric(crossprod(G_v, W_hat_fixed %*% G_v))
    if (!is.finite(v)) 1e8 else v
  }
  grad_fixed <- function(p) {
    g_m <- g_func(p)
    G_v <- matrix(.colMeans(g_m, n, esteq_dim), esteq_dim, 1)
    dg_arr <- dg_func(p)
    Gamma_mat <- matrix(0, esteq_dim, param_dim)
    for (j in 1:param_dim) Gamma_mat[, j] <- .colMeans(dg_arr[, j, , drop = FALSE], n, esteq_dim)
    2 * as.vector(crossprod(Gamma_mat, W_hat_fixed %*% G_v))
  }

  # Init step
  opt0 <- optim(estimates, obj_fixed, gr = grad_fixed, method = "L-BFGS-B", control = list(maxit = 1000))
  estimates <- opt0$par

  cat(sprintf("\nAfter init (W=I): ||a||=%.3f\n", sqrt(sum(estimates^2))))
  cat(sprintf("%4s %12s %12s %12s %12s\n", "iter", "cue_grad", "cue_obj", "||alpha||", "delta_alpha"))

  tol <- 1e-6
  for (t in 1:200) {
    g_mat <- g_func(estimates)
    W_hat_fixed <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(esteq_dim))

    info <- cue_obj_and_grad(estimates)

    if (t <= 20 || t %% 25 == 0 || info$grad_norm < tol) {
      cat(sprintf("%4d %12.3e %12.4e %12.4f\n", t, info$grad_norm, info$obj, sqrt(sum(estimates^2))))
    }

    if (info$grad_norm < tol) {
      cat("  ** CONVERGED **\n")
      break
    }

    prev <- estimates
    opt_t <- optim(estimates, obj_fixed, gr = grad_fixed, method = "L-BFGS-B", control = list(maxit = 1000))
    estimates <- opt_t$par
  }
  if (t == 200) cat("  ** NOT CONVERGED in 200 iters **\n")
}
