setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

n_val <- 2000
ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

# Known degenerate reps from earlier diagnosis
degen_reps <- c(34, 35, 38, 68, 109)

model_idx <- 3

for (rep_i in degen_reps) {
  cat(sprintf("\n========== Rep %d ==========\n", rep_i))
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
  r_vec <- dat[["r"]]
  n <- n_val

  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL),
    optimizer = "constrained_nr"
  )

  # Build components manually to trace
  ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
  gmm_cnr <- ebmr$ps_fit.list[[1]]$gmm_fit
  design_mat <- ebmr$ps_fit.list[[1]]$design_matrix
  h_x <- ebmr$ps_fit.list[[1]]$h_x

  cat(sprintf("  CNR final: grad=%.2e, cond=%.2e, obj=%.6f, conv=%s, iter=%d\n",
      gmm_cnr$opt$final_grad_norm, gmm_cnr$opt$solution_cond,
      gmm_cnr$opt$objective, gmm_cnr$opt$converged, gmm_cnr$opt$iterations))
  cat(sprintf("  alpha = %s\n", paste(round(gmm_cnr$estimates, 4), collapse=", ")))

  pi_cnr <- ebmr$ps_fit.list[[1]]$fitted.values
  cat(sprintf("  pi: min=%.6f, max=%.6f, mean=%.4f\n", min(pi_cnr), max(pi_cnr), mean(pi_cnr)))
  cat(sprintf("  # pi > 0.999: %d, # pi > 0.9999: %d, # pi > 0.99999: %d\n",
      sum(pi_cnr > 0.999), sum(pi_cnr > 0.9999), sum(pi_cnr > 0.99999)))

  # Now trace what happened: run constrained_nr step by step
  # We need access to g, Gamma, W functions
  # Reconstruct from the package internals
  compute_pi <- function(eta) plogis(-eta)  # logistic_complement

  g_fn <- function(alpha) {
    eta <- as.vector(design_mat %*% alpha)
    pi_v <- compute_pi(eta)
    (r_vec / pi_v - 1) * h_x
  }

  Gamma_fn <- function(alpha) {
    eta <- as.vector(design_mat %*% alpha)
    pi_v <- compute_pi(eta)
    cf <- r_vec * (1 - pi_v) / pi_v
    crossprod(h_x * cf, design_mat) / n
  }

  cond_fn <- function(alpha) {
    Gamma_hat <- Gamma_fn(alpha)
    g_mat <- g_fn(alpha)
    W_hat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(ncol(h_x)))
    H <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
    eig <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    max(eig) / max(min(eig), 1e-15)
  }

  grad_fn <- function(alpha) {
    g_mat <- g_fn(alpha)
    G_vec <- matrix(colMeans(g_mat), ncol(h_x), 1)
    W_hat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(ncol(h_x)))
    Gamma_hat <- Gamma_fn(alpha)
    max(abs(2 * as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))))
  }

  obj_fn <- function(alpha) {
    g_mat <- g_fn(alpha)
    G_vec <- matrix(colMeans(g_mat), ncol(h_x), 1)
    W_hat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(ncol(h_x)))
    as.numeric(crossprod(G_vec, W_hat %*% G_vec))
  }

  # Trace the objective landscape along the path from init=0 to the degenerate solution
  alpha_degen <- gmm_cnr$estimates
  cat("\n  === Landscape from 0 to degenerate solution ===\n")
  for (frac in c(0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0)) {
    alpha_t <- frac * alpha_degen
    pi_t <- compute_pi(as.vector(design_mat %*% alpha_t))
    tryCatch({
      cat(sprintf("    t=%.1f: obj=%.6f, grad=%.2e, cond=%.2e, max_pi=%.6f, alpha2=%.3f\n",
          frac, obj_fn(alpha_t), grad_fn(alpha_t), cond_fn(alpha_t),
          max(pi_t), alpha_t[2]))
    }, error = function(e) {
      cat(sprintf("    t=%.1f: ERROR: %s\n", frac, e$message))
    })
  }

  # Also compare: L-BFGS-B solution
  single_ps$optimizer <- "L-BFGS-B"
  ebmr_lb <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
  gmm_lb <- ebmr_lb$ps_fit.list[[1]]$gmm_fit
  cat(sprintf("\n  L-BFGS-B: grad=%.2e, cond=%.2e, obj=%.6f\n",
      gmm_lb$opt$final_grad_norm, gmm_lb$opt$solution_cond, gmm_lb$opt$objective))
  cat(sprintf("  alpha = %s\n", paste(round(gmm_lb$estimates, 4), collapse=", ")))

  # Check: is the constrained_nr solution the same as L-BFGS-B or different?
  cat(sprintf("  ||alpha_cnr - alpha_lb|| = %.4f\n",
      sqrt(sum((gmm_cnr$estimates - gmm_lb$estimates)^2))))

  # Where does cond cross 1e8?
  cat("\n  === Finding cond=1e8 crossing ===\n")
  # Binary search along alpha[2] dimension, keeping others at cnr values
  alpha_base <- gmm_cnr$estimates
  for (a2 in seq(-3, alpha_base[2], length.out = 10)) {
    alpha_test <- alpha_base
    alpha_test[2] <- a2
    pi_test <- compute_pi(as.vector(design_mat %*% alpha_test))
    tryCatch({
      cc <- cond_fn(alpha_test)
      gg <- grad_fn(alpha_test)
      cat(sprintf("    alpha2=%.2f: cond=%.2e, grad=%.2e, max_pi=%.6f\n",
          a2, cc, gg, max(pi_test)))
    }, error = function(e) {
      cat(sprintf("    alpha2=%.2f: ERROR\n", a2))
    })
  }
}
