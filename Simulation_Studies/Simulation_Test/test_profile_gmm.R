setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})

n_val <- 2000
n_reps <- 200
ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
mu_true <- get_mu_true("setting4")

formula_j <- ps_spec[["formula.list"]][[3]]
h_alpha_vars <- ps_spec[["h_alpha.list"]][[3]]

mu_vals <- rep(NA_real_, n_reps)
alpha_mat <- matrix(NA, n_reps, 4)
conv_vals <- rep(NA, n_reps)
grad_vals <- rep(NA_real_, n_reps)
iter_vals <- rep(NA_integer_, n_reps)

for (rep_i in 1:n_reps) {
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
  tryCatch({
    r_vec <- dat[["r"]]
    design_mat <- model.matrix(formula_j, data=dat)
    h_alpha_mat <- as.matrix(dat[, h_alpha_vars, drop=FALSE])
    h_full <- cbind(1, h_alpha_mat)
    n <- n_val
    p <- ncol(design_mat)  # 4: intercept, y, u2, z2
    h_dim <- ncol(h_full)  # 5: 1, u1, u2, z1, z2
    r_mean <- mean(r_vec)

    # X_rest: design matrix columns excluding intercept
    X_rest <- design_mat[, -1, drop=FALSE]

    # Solve intercept from mean(r/pi) = 1 given alpha_rest
    # For logistic pi = 1/(1+exp(-eta)): alpha0 = log(mean(r*exp(-X_rest*alpha_rest))) - log(1-mean(r))
    solve_intercept <- function(alpha_rest) {
      eta_rest <- as.vector(X_rest %*% alpha_rest)
      num <- mean(r_vec * exp(-eta_rest))
      denom <- 1 - r_mean
      if (num <= 0 || denom <= 0) return(0)
      log(num) - log(denom)
    }

    # Full alpha from alpha_rest
    full_alpha <- function(alpha_rest) {
      c(solve_intercept(alpha_rest), alpha_rest)
    }

    # Model
    model_fn <- function(alpha) {
      eta <- as.vector(design_mat %*% alpha)
      1 / (1 + exp(-eta))
    }

    # Moment function g_i(alpha) = (r/pi - 1) * h
    Phi_alpha <- function(alpha) {
      pi_hat <- model_fn(alpha)
      (r_vec / pi_hat - 1) * h_full
    }

    # Profile moment: g as function of alpha_rest
    Phi_profile <- function(alpha_rest) {
      alpha <- full_alpha(alpha_rest)
      Phi_alpha(alpha)
    }

    # Profile Gamma (Jacobian) via numerical differentiation
    G_profile <- function(alpha_rest) {
      g_mat <- Phi_profile(alpha_rest)
      colMeans(g_mat)
    }

    Gamma_profile <- function(alpha_rest) {
      G0 <- G_profile(alpha_rest)
      p_rest <- length(alpha_rest)
      Gamma_mat <- matrix(0, h_dim, p_rest)
      eps <- 1e-7
      for (j in 1:p_rest) {
        e_j <- rep(0, p_rest); e_j[j] <- eps
        Gamma_mat[, j] <- (G_profile(alpha_rest + e_j) - G0) / eps
      }
      Gamma_mat
    }

    # Iterative GMM on alpha_rest
    W_func <- function(g_mat) solve(crossprod(g_mat) / nrow(g_mat))

    alpha_rest <- rep(0, p - 1)
    max_outer <- 200
    tol <- 1e-6
    best_grad <- Inf
    best_est <- alpha_rest

    for (t in 1:max_outer) {
      g_mat <- Phi_profile(alpha_rest)
      G_vec <- colMeans(g_mat)
      W_hat <- tryCatch(W_func(g_mat), error=function(e) diag(h_dim))

      # Convergence gradient
      Gamma_hat <- Gamma_profile(alpha_rest)
      conv_grad <- 2 * as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))
      grad_norm <- max(abs(conv_grad))

      if (grad_norm < best_grad) {
        best_grad <- grad_norm
        best_est <- alpha_rest
      }

      if (grad_norm < tol) break

      # Inner minimization: minimize Q(alpha_rest) = G'WG for fixed W
      obj_fn <- function(a_rest) {
        g_m <- Phi_profile(a_rest)
        G_v <- colMeans(g_m)
        as.numeric(t(G_v) %*% W_hat %*% G_v)
      }
      grad_fn <- function(a_rest) {
        g_m <- Phi_profile(a_rest)
        G_v <- colMeans(g_m)
        Gam <- Gamma_profile(a_rest)
        2 * as.vector(crossprod(Gam, W_hat %*% G_v))
      }
      opt <- optim(alpha_rest, obj_fn, gr = grad_fn, method = "L-BFGS-B",
                   control = list(maxit = 1000))
      alpha_rest <- opt$par
    }

    alpha_final <- full_alpha(best_est)
    pi_hat <- model_fn(alpha_final)

    # HT estimator
    mu_vals[rep_i] <- mean(r_vec * dat[["y"]] / pi_hat)
    alpha_mat[rep_i, ] <- alpha_final
    conv_vals[rep_i] <- best_grad < tol
    grad_vals[rep_i] <- best_grad
    iter_vals[rep_i] <- t
  }, error = function(e) NULL)
}

valid <- !is.na(mu_vals)
conv <- conv_vals[valid]

cat(sprintf("mu.true = %.4f\n", mu_true))
cat(sprintf("Total valid: %d, converged: %d, not converged: %d\n\n", sum(valid), sum(conv), sum(!conv)))

esd <- sd(mu_vals[valid])
cat(sprintf("=== Overall ===\n"))
cat(sprintf("Bias=%.4f, ESD=%.4f\n\n", mean(mu_vals[valid]) - mu_true, esd))

cat(sprintf("=== alpha ===\n"))
colnames(alpha_mat) <- c("intercept", "y", "u2", "z2")
for (j in 1:4) {
  cat(sprintf("alpha[%d] (%s): mean=%.4f, sd=%.4f, range=[%.4f, %.4f]\n",
      j, colnames(alpha_mat)[j],
      mean(alpha_mat[valid, j]), sd(alpha_mat[valid, j]),
      min(alpha_mat[valid, j]), max(alpha_mat[valid, j])))
}

cat(sprintf("\n=== Convergence ===\n"))
cat(sprintf("Converged: %d, Not: %d\n", sum(conv), sum(!conv)))
cat(sprintf("Grad norm: mean=%.2e, max=%.2e\n", mean(grad_vals[valid]), max(grad_vals[valid])))

cat(sprintf("\n=== mean(pi) ===\n"))
# Check mean pi at final estimates
pi_means <- rep(NA, n_reps)
for (i in which(valid)) {
  dat <- all_data[((i-1)*n_val + 1):(i*n_val), ]
  dm <- model.matrix(formula_j, data=dat)
  eta <- as.vector(dm %*% alpha_mat[i,])
  pi_means[i] <- mean(1/(1+exp(-eta)))
}
cat(sprintf("mean(pi): mean=%.4f, sd=%.4f, range=[%.4f, %.4f]\n",
    mean(pi_means[valid]), sd(pi_means[valid]),
    min(pi_means[valid]), max(pi_means[valid])))
cat(sprintf("mean(r) = %.4f\n", mean(all_data[["r"]][1:n_val])))
