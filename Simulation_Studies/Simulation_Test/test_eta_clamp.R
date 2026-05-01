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

test_reps <- c(929, 438, 430, 1, 2, 3)

# For each eta_max, manually run GMM with that clamp
for (eta_max_val in c(20, 10, 5, 3)) {
  cat(sprintf("=== eta_max = %d  (min pi = %.2e) ===\n", eta_max_val,
              1 / (1 + exp(eta_max_val))))

  compute_pi_test <- function(eta) {
    eta <- pmin(pmax(eta, -eta_max_val), eta_max_val)
    1 / (1 + exp(eta))
  }

  for (rep_i in test_reps) {
    dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]
    label <- if (rep_i %in% c(929, 438, 430)) "OUT" else "NORM"

    # Get h_x and design_matrix from the package fit (use current eta_max=20)
    ebmr <- EBMRAlgorithmFast4$new("y", subset_ps_spec, dat, W_func)
    fit <- ebmr$ps_fit.list[[1]]
    h_x <- fit$h_x
    design_matrix <- fit$design_matrix
    r_vec <- as.numeric(dat$r)
    n <- nrow(dat)
    esteq_dim <- ncol(h_x)
    param_dim <- ncol(design_matrix)

    g_func <- function(param) {
      eta <- design_matrix %*% param
      pi_vec <- compute_pi_test(eta)
      rw <- r_vec / as.vector(pi_vec)
      (rw - 1) * h_x
    }

    G_func <- function(param) {
      g_mat <- g_func(param)
      matrix(colMeans(g_mat), esteq_dim, 1)
    }

    # Two-step GMM with L-BFGS-B, init = 0
    obj1 <- function(param) {
      G_bar <- G_func(param)
      as.numeric(crossprod(G_bar))
    }
    gr1 <- NULL  # use numerical gradient for simplicity
    opt1 <- optim(rep(0, param_dim), obj1, method = "L-BFGS-B",
                  lower = rep(-Inf, param_dim), upper = rep(Inf, param_dim),
                  control = list(maxit = 1000))

    estimates <- opt1$par
    for (t in 1:200) {
      g_mat <- g_func(estimates)
      W_hat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(esteq_dim))

      obj2 <- function(param) {
        G_bar <- G_func(param)
        as.numeric(crossprod(G_bar, W_hat %*% G_bar))
      }
      opt_t <- optim(estimates, obj2, method = "L-BFGS-B",
                     lower = rep(-Inf, param_dim), upper = rep(Inf, param_dim),
                     control = list(maxit = 1000))
      if (max(abs(opt_t$par - estimates)) < 1e-6) break
      estimates <- opt_t$par
    }

    pi_hat <- as.vector(compute_pi_test(design_matrix %*% estimates))
    mu_hat <- mean(r_vec / pi_hat * dat$y)

    cat(sprintf("  %-6d %-5s mu=%-10.4f min_pi=%-10.6f max_r/pi=%-10.1f alpha=[%s]\n",
                rep_i, label, mu_hat, min(pi_hat), max(r_vec / pi_hat),
                paste(round(estimates, 4), collapse=", ")))
  }
  cat("\n")
}

cat("Done!\n")
