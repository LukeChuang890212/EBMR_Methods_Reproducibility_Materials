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

test_reps <- c(1, 3, 7, 15, 34)

for (rep_i in test_reps) {
  cat(sprintf("\n========== Rep %d ==========\n", rep_i))
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]

  r_vec <- dat[["r"]]
  design_mat <- model.matrix(formula_j, data=dat)
  h_alpha_mat <- as.matrix(dat[, h_alpha_vars, drop=FALSE])
  h_full <- cbind(1, h_alpha_mat)
  n <- n_val
  p <- ncol(design_mat)

  model_fn <- function(alpha) {
    eta <- as.vector(design_mat %*% alpha)
    1 / (1 + exp(-eta))
  }

  Phi_alpha <- function(alpha) {
    pi_hat <- model_fn(alpha)
    (r_vec / pi_hat - 1) * h_full
  }

  # Objective Q(alpha, W(alpha)) = G(alpha)' W(alpha) G(alpha)
  Q_full <- function(alpha) {
    g_mat <- Phi_alpha(alpha)
    G <- colMeans(g_mat)
    W <- tryCatch(solve(crossprod(g_mat)/n), error=function(e) diag(ncol(h_full)))
    as.numeric(t(G) %*% W %*% G)
  }

  # Objective Q(alpha, W_fixed) for fixed W
  Q_fixed <- function(alpha, W_fixed) {
    g_mat <- Phi_alpha(alpha)
    G <- colMeans(g_mat)
    as.numeric(t(G) %*% W_fixed %*% G)
  }

  # One-step NR iteration (same as the best algorithm)
  Gamma_fn <- function(alpha) {
    pi_hat <- model_fn(alpha)
    neg_r_pi2 <- -r_vec / (pi_hat^2)
    dpi <- design_mat * (pi_hat * (1 - pi_hat))
    crossprod(h_full * neg_r_pi2, dpi) / n
  }

  alpha_t <- rep(0, p)
  cat(sprintf("%4s  %12s  %12s  %12s  %12s  %8s\n",
      "t", "Q(a,W(a))", "Q(a,W_t)", "Q(a+,W_t)", "||step||", "alpha[2]"))

  for (t in 1:30) {
    g_mat <- Phi_alpha(alpha_t)
    G_t <- colMeans(g_mat)
    W_t <- tryCatch(solve(crossprod(g_mat)/n), error=function(e) diag(ncol(h_full)))
    Q_curr_full <- as.numeric(t(G_t) %*% W_t %*% G_t)
    Q_curr_fixed <- Q_curr_full  # same since W_t = W(alpha_t)

    # Gradient of Q(alpha, W_t) w.r.t. alpha
    Gamma_hat <- Gamma_fn(alpha_t)
    g_vec <- 2 * as.vector(crossprod(Gamma_hat, W_t %*% G_t))

    # Gauss-Newton direction
    H_gn <- crossprod(Gamma_hat, W_t %*% Gamma_hat)
    direction <- tryCatch(solve(H_gn, -g_vec), error=function(e) -g_vec)

    # Trust region
    step_norm <- sqrt(sum(direction^2))
    if (step_norm > 0.5) direction <- direction * (0.5 / step_norm)

    # Backtracking
    step <- 1.0
    for (ls in 1:20) {
      candidate <- alpha_t + step * direction
      obj_new <- Q_fixed(candidate, W_t)
      if (is.finite(obj_new) && obj_new < Q_curr_fixed - 1e-4 * step * sum(g_vec * direction)) break
      step <- step * 0.5
    }
    alpha_new <- alpha_t + step * direction
    Q_new_fixed <- Q_fixed(alpha_new, W_t)
    actual_step <- sqrt(sum((alpha_new - alpha_t)^2))

    cat(sprintf("%4d  %12.6f  %12.6f  %12.6f  %12.6f  %8.4f\n",
        t, Q_curr_full, Q_curr_fixed, Q_new_fixed, actual_step, alpha_t[2]))

    alpha_t <- alpha_new
  }
}
