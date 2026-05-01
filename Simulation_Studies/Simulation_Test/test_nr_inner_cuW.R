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

for (rep_i in 1:n_reps) {
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
  tryCatch({
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

    # NR solver with simultaneous W(alpha) update at each step
    # + trust region + local preference
    nr_cu_inner <- function(start, max_inner = 200L) {
      alpha <- start
      init_inner <- start
      best_inner_grad <- Inf
      best_inner <- alpha

      for (k in 1:max_inner) {
        g_mat <- Phi_alpha(alpha)
        G_vec <- colMeans(g_mat)
        W_hat <- tryCatch(W_func(g_mat), error = function(e) diag(h_dim))

        Gamma_hat <- Gamma_fn(alpha)
        g_vec <- 2 * as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))
        inner_grad <- max(abs(g_vec))

        # Track best within distance 5 of inner start
        d <- sqrt(sum((alpha - init_inner)^2))
        if (d < 5 && inner_grad < best_inner_grad) {
          best_inner_grad <- inner_grad
          best_inner <- alpha
        }

        if (inner_grad < 1e-8) break

        H_gn <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
        direction <- tryCatch(solve(H_gn, -g_vec), error = function(e) -g_vec)

        # Trust region: cap step norm to 0.5
        step_norm <- sqrt(sum(direction^2))
        if (step_norm > 0.5) direction <- direction * (0.5 / step_norm)

        # Backtracking line search on Q(alpha, W(alpha_current))
        obj_cur <- as.numeric(t(G_vec) %*% W_hat %*% G_vec)
        step <- 1.0
        for (ls in 1:20) {
          candidate <- alpha + step * direction
          g_new <- Phi_alpha(candidate)
          G_new <- colMeans(g_new)
          obj_new <- as.numeric(t(G_new) %*% W_hat %*% G_new)
          if (is.finite(obj_new) && obj_new < obj_cur - 1e-4 * step * sum(g_vec * direction)) break
          step <- step * 0.5
        }
        alpha <- alpha + step * direction
      }

      # If drifted far, use local best
      d_final <- sqrt(sum((alpha - init_inner)^2))
      if (d_final > 5 && best_inner_grad < 0.1) alpha <- best_inner
      alpha
    }

    # Single call: full NR minimization with CU-style W updates
    estimates <- nr_cu_inner(rep(0, p))

    # Check convergence
    g_mat <- Phi_alpha(estimates)
    G_vec <- colMeans(g_mat)
    W_hat <- tryCatch(W_func(g_mat), error = function(e) diag(h_dim))
    Gamma_hat <- Gamma_fn(estimates)
    conv_grad <- 2 * as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))
    grad_norm <- max(abs(conv_grad))
    converged <- grad_norm < 1e-6

    pi_hat <- model_fn(estimates)
    mu_vals[rep_i] <- mean(r_vec * dat[["y"]] / pi_hat)
    alpha_mat[rep_i, ] <- estimates
    conv_vals[rep_i] <- converged
    grad_vals[rep_i] <- grad_norm
  }, error = function(e) NULL)
}

valid <- !is.na(mu_vals)
conv <- conv_vals[valid]

cat(sprintf("mu.true = %.4f\n", mu_true))
cat(sprintf("Total valid: %d, converged: %d, not converged: %d\n\n", sum(valid), sum(conv), sum(!conv)))

esd <- sd(mu_vals[valid])
cat(sprintf("=== Overall ===\n"))
cat(sprintf("Bias=%.4f, ESD=%.4f\n\n", mean(mu_vals[valid]) - mu_true, esd))

colnames(alpha_mat) <- c("intercept", "y", "u2", "z2")
cat(sprintf("=== alpha ===\n"))
for (j in 1:4) {
  a <- alpha_mat[valid, j]
  cat(sprintf("alpha[%d] (%s): mean=%.4f, sd=%.4f, range=[%.4f, %.4f]\n",
      j, colnames(alpha_mat)[j], mean(a), sd(a), min(a), max(a)))
}

cat(sprintf("\n=== Grad norm ===\n"))
cat(sprintf("Converged:     mean=%.2e, max=%.2e\n",
    mean(grad_vals[valid][conv]), max(grad_vals[valid][conv])))
if (any(!conv)) {
  cat(sprintf("Not converged: mean=%.2e, max=%.2e\n",
      mean(grad_vals[valid][!conv]), max(grad_vals[valid][!conv])))
}

cat(sprintf("\n=== mean(pi) ===\n"))
pi_means <- rep(NA, n_reps)
for (i in which(valid)) {
  dat <- all_data[((i-1)*n_val + 1):(i*n_val), ]
  dm <- model.matrix(formula_j, data=dat)
  eta <- as.vector(dm %*% alpha_mat[i,])
  pi_means[i] <- mean(1/(1+exp(-eta)))
}
cat(sprintf("mean=%.4f, sd=%.4f, range=[%.4f, %.4f]\n",
    mean(pi_means[valid]), sd(pi_means[valid]),
    min(pi_means[valid]), max(pi_means[valid])))
