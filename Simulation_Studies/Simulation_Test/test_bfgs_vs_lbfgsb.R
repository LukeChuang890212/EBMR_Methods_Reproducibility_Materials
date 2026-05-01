setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
library(EBMRalgorithmFast4)

data_file <- misspecified_model_all_data_file.list$setting3$miss50[[1]]
all_data <- readRDS(data_file)
n_val <- 2000

ps_spec <- get_ps_spec("9-alt1")
subset_ps_spec <- list(
  formula.list = ps_spec$formula.list[2],
  h_alpha.list = ps_spec$h_alpha.list[2],
  inv_link = ps_spec$inv_link,
  outcome = ps_spec$outcome
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

test_reps <- c(929, 438, 430, 221, 180, 1, 2, 3, 4, 5)

cat(sprintf("%-6s %-8s %-8s %-10s %-12s %-40s\n",
            "Rep", "Type", "Method", "mu_ipw", "Objective", "Alpha"))

for (rep_i in test_reps) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]
  label <- if (rep_i %in% c(929, 438, 430, 221, 180)) "OUT" else "NORM"

  # Get L-BFGS-B result from package (current)
  ebmr <- EBMRAlgorithmFast4$new("y", subset_ps_spec, dat, W_func)
  fit <- ebmr$ps_fit.list[[1]]
  mu_lbfgsb <- mean(as.numeric(dat$r) / fit$fitted.values * dat$y)

  cat(sprintf("%-6d %-8s %-8s %-10.4f %-12.6e [%s]\n",
              rep_i, label, "L-BFGS-B", mu_lbfgsb, fit$gmm_fit$opt$objective,
              paste(round(fit$coefficients, 4), collapse=", ")))

  # Now run BFGS manually using same setup
  h_x <- fit$h_x
  design_matrix <- fit$design_matrix
  r_vec <- as.numeric(dat$r)
  n <- nrow(dat)
  esteq_dim <- ncol(h_x)
  param_dim <- ncol(design_matrix)

  eta_max <- 20
  compute_pi <- function(eta) {
    eta <- pmin(pmax(eta, -eta_max), eta_max)
    1 / (1 + exp(eta))
  }

  g_func <- function(param) {
    eta <- design_matrix %*% param
    pi_vec <- compute_pi(eta)
    rw <- r_vec / as.vector(pi_vec)
    (rw - 1) * h_x
  }

  G_func <- function(param) {
    g_mat <- g_func(param)
    matrix(colMeans(g_mat), esteq_dim, 1)
  }

  # Two-step GMM with BFGS
  # Step 1: W = I
  obj1 <- function(param) {
    G_bar <- G_func(param)
    as.numeric(crossprod(G_bar))
  }
  opt1 <- optim(rep(0, param_dim), obj1, method = "BFGS", control = list(maxit = 1000))

  # Step 2: iterative W updates
  estimates <- opt1$par
  for (t in 1:200) {
    g_mat <- g_func(estimates)
    W_hat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(esteq_dim))

    obj2 <- function(param) {
      G_bar <- G_func(param)
      as.numeric(crossprod(G_bar, W_hat %*% G_bar))
    }
    opt_t <- optim(estimates, obj2, method = "BFGS", control = list(maxit = 1000))
    if (max(abs(opt_t$par - estimates)) < 1e-6) break
    estimates <- opt_t$par
  }

  # Evaluate final objective
  g_final <- g_func(estimates)
  G_final <- matrix(colMeans(g_final), esteq_dim, 1)
  W_final <- tryCatch(solve(crossprod(g_final) / n), error = function(e) diag(esteq_dim))
  obj_final <- as.numeric(crossprod(G_final, W_final %*% G_final))

  pi_bfgs <- as.vector(compute_pi(design_matrix %*% estimates))
  mu_bfgs <- mean(r_vec / pi_bfgs * dat$y)

  cat(sprintf("%-6d %-8s %-8s %-10.4f %-12.6e [%s]\n\n",
              rep_i, label, "BFGS", mu_bfgs, obj_final,
              paste(round(estimates, 4), collapse=", ")))
}

cat("Done!\n")
