setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})

n_val <- 2000
ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)

formula_j <- ps_spec[["formula.list"]][[3]]
h_alpha_vars <- ps_spec[["h_alpha.list"]][[3]]

# Trace a few reps through the iterative GMM to see when/if drift happens
test_reps <- c(1, 3, 7, 15, 34, 50)

for (rep_i in test_reps) {
  cat(sprintf("\n========== Rep %d ==========\n", rep_i))
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]

  r_vec <- dat[["r"]]
  design_mat <- model.matrix(formula_j, data=dat)
  h_alpha_mat <- as.matrix(dat[, h_alpha_vars, drop=FALSE])
  h_full <- cbind(1, h_alpha_mat)
  n <- n_val
  p <- ncol(design_mat)
  h_dim <- ncol(h_full)

  model_fn <- function(alpha) {
    eta <- as.vector(design_mat %*% alpha)
    1 / (1 + exp(-eta))
  }

  Phi_alpha <- function(alpha) {
    pi_hat <- model_fn(alpha)
    (r_vec / pi_hat - 1) * h_full
  }

  W_func <- function(g_mat) solve(crossprod(g_mat) / nrow(g_mat))

  Gamma_fn <- function(alpha) {
    pi_hat <- model_fn(alpha)
    common <- -r_vec * (1 - pi_hat) / pi_hat
    crossprod(h_full * common, design_mat) / n
  }

  # Step 1: W=I
  obj_I <- function(alpha) {
    G <- colMeans(Phi_alpha(alpha))
    sum(G^2)
  }
  opt0 <- optim(rep(0, p), obj_I, method = "L-BFGS-B", control = list(maxit = 5000))
  estimates <- opt0$par

  cat(sprintf("  Step1: alpha[2]=%.4f, mean(pi)=%.4f\n", estimates[2], mean(model_fn(estimates))))

  # Iterative steps
  for (t in 1:30) {
    g_mat <- Phi_alpha(estimates)
    G_vec <- colMeans(g_mat)
    W_hat <- tryCatch(W_func(g_mat), error = function(e) diag(h_dim))

    Gamma_hat <- Gamma_fn(estimates)
    conv_grad <- 2 * as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))
    grad_norm <- max(abs(conv_grad))

    obj_fixedW <- function(alpha) {
      G <- colMeans(Phi_alpha(alpha))
      as.numeric(t(G) %*% W_hat %*% G)
    }
    opt <- optim(estimates, obj_fixedW, method = "L-BFGS-B", control = list(maxit = 1000))
    estimates <- opt$par

    cat(sprintf("  t=%2d: alpha[2]=%8.4f, mean(pi)=%.4f, grad=%.2e, obj=%.6f\n",
        t, estimates[2], mean(model_fn(estimates)), grad_norm, opt$value))

    if (grad_norm < 1e-6) { cat("  CONVERGED\n"); break }
  }
}
