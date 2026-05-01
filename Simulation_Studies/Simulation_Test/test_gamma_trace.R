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

# Trace ||Gamma||, cond(Gamma'WGamma), mean(pi) through iterations
# for a "good" rep and a "bad" rep
for (rep_i in c(1, 15)) {
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

  Gamma_fn <- function(alpha) {
    pi_hat <- model_fn(alpha)
    common <- -r_vec * (1 - pi_hat) / pi_hat
    crossprod(h_full * common, design_mat) / n
  }

  W_func <- function(g_mat) solve(crossprod(g_mat) / nrow(g_mat))

  # Check at alpha = 0
  gamma_zero <- sqrt(sum(Gamma_fn(rep(0, p))^2))
  cat(sprintf("||Gamma|| at alpha=0: %.4f\n", gamma_zero))

  # Step 1: W=I
  obj_I <- function(alpha) { G <- colMeans(Phi_alpha(alpha)); sum(G^2) }
  opt0 <- optim(rep(0, p), obj_I, method = "L-BFGS-B", control = list(maxit = 5000))
  estimates <- opt0$par

  gamma_s1 <- sqrt(sum(Gamma_fn(estimates)^2))
  cat(sprintf("Step1: alpha[2]=%.4f, mean(pi)=%.4f, ||Gamma||=%.4f, ratio=%.3f\n",
      estimates[2], mean(model_fn(estimates)), gamma_s1, gamma_s1/gamma_zero))

  # Trace iterative GMM
  cat(sprintf("%4s  %8s  %8s  %10s  %10s  %10s  %12s\n",
      "t", "alpha[2]", "mean(pi)", "||Gamma||", "ratio", "cond(H)", "grad"))

  for (t in 1:20) {
    g_mat <- Phi_alpha(estimates)
    G_vec <- colMeans(g_mat)
    W_hat <- tryCatch(W_func(g_mat), error = function(e) diag(h_dim))

    Gamma_hat <- Gamma_fn(estimates)
    gamma_norm <- sqrt(sum(Gamma_hat^2))
    H_gn <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
    eig_H <- eigen(H_gn, symmetric = TRUE, only.values = TRUE)$values
    cond_H <- max(eig_H) / max(min(eig_H), 1e-15)

    conv_grad <- 2 * as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))
    grad_norm <- max(abs(conv_grad))

    cat(sprintf("%4d  %8.4f  %8.4f  %10.4f  %10.3f  %10.2e  %12.4e\n",
        t, estimates[2], mean(model_fn(estimates)), gamma_norm,
        gamma_norm/gamma_zero, cond_H, grad_norm))

    obj_fn <- function(alpha) {
      G <- colMeans(Phi_alpha(alpha))
      as.numeric(t(G) %*% W_hat %*% G)
    }
    opt <- optim(estimates, obj_fn, method = "L-BFGS-B", control = list(maxit = 1000))
    estimates <- opt$par
  }
}
