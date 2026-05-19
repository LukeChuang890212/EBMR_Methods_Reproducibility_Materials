## Test different inner-loop optimizers for iterative CUE
## Keep the outer loop (update W, then minimize G'WG for fixed W),
## but try different optimizers for the inner minimization step
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)
library(numDeriv)

ps_spec_base <- get_ps_spec("9-alt1")
h_alpha_fn <- function(dat) cbind(
  u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2,
  u1u2 = dat$u1*dat$u2, z1z2 = dat$z1*dat$z2,
  u1z1 = dat$u1*dat$z1, z1u2 = dat$z1*dat$u2,
  u1z2 = dat$u1*dat$z2, u2z2 = dat$u2*dat$z2)

data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]]
all_data <- NULL
for (fi in seq_along(data_file)) {
  if (file.exists(data_file[[fi]])) {
    test <- readRDS(data_file[[fi]])
    if (nrow(test) / 1000 == 2000) { all_data <- test; break }
  }
}

dat <- all_data[1:2000, ]
design_mat <- model.matrix(ps_spec_base$formula.list[[3]], dat)
r_vec <- as.vector(dat$r)
n <- 2000
h_x <- cbind(1, h_alpha_fn(dat))
h_dim <- ncol(h_x)
alpha_dim <- ncol(design_mat)

compute_pi <- function(eta) 1/(1+exp(eta))

g_fn <- function(alpha) {
  eta <- as.vector(design_mat %*% alpha)
  eta <- pmin(pmax(eta, -20), 20)
  pi_vec <- compute_pi(eta)
  as.vector(r_vec/pi_vec - 1) * h_x
}
G_fn <- function(alpha) colMeans(g_fn(alpha))

# Inner objective for fixed W: Q_W(alpha) = G(alpha)' W G(alpha)
make_inner_obj <- function(W_fixed) {
  function(alpha) {
    G <- matrix(G_fn(alpha), ncol = 1)
    as.numeric(t(G) %*% W_fixed %*% G)
  }
}

# Run iterative CUE with a given inner optimizer
run_iterative_cue <- function(inner_optimizer, max_outer = 100, stall_limit = 20) {
  alpha <- rep(0, alpha_dim)
  best_grad <- Inf
  best_alpha <- alpha
  no_improve <- 0L

  for (t in 1:max_outer) {
    g_mat <- g_fn(alpha)
    W_hat <- tryCatch(solve(var(g_mat)), error = function(e) diag(h_dim))

    # Inner optimization with the specified optimizer
    alpha <- inner_optimizer(alpha, W_hat)

    # Check convergence: gradient of CUE obj at current alpha with current W
    G_new <- G_fn(alpha)
    g_new <- g_fn(alpha)
    W_new <- tryCatch(solve(var(g_new)), error = function(e) diag(h_dim))
    # The "gradient" we care about is the iterative stationarity:
    # 2 * Gamma' W G where W is evaluated at current alpha
    grad_norm <- max(abs(G_new))  # simplified convergence check

    if (grad_norm < best_grad) {
      best_grad <- grad_norm
      best_alpha <- alpha
      no_improve <- 0L
    } else {
      no_improve <- no_improve + 1L
      if (no_improve >= stall_limit) {
        alpha <- best_alpha
        break
      }
    }
    if (grad_norm < 1e-6) break
  }

  g_final <- g_fn(best_alpha)
  G_final <- G_fn(best_alpha)
  W_final <- tryCatch(solve(var(g_final)), error = function(e) diag(h_dim))
  obj_final <- as.numeric(t(G_final) %*% W_final %*% G_final)

  list(alpha = best_alpha, grad_norm = best_grad, obj = obj_final, iter = t)
}

cat("=== Iterative CUE with different inner-loop optimizers ===\n\n")

# 1. L-BFGS-B (current default)
cat("--- 1. L-BFGS-B (current) ---\n")
inner_lbfgsb <- function(alpha, W) {
  obj <- make_inner_obj(W)
  opt <- optim(alpha, obj, method = "L-BFGS-B", control = list(maxit = 1000))
  opt$par
}
res1 <- run_iterative_cue(inner_lbfgsb)
cat(sprintf("   iter=%d, obj=%.6e, max|G|=%.4e, alpha=(%s)\n\n",
    res1$iter, res1$obj, res1$grad_norm, paste(round(res1$alpha, 4), collapse=", ")))

# 2. BFGS
cat("--- 2. BFGS ---\n")
inner_bfgs <- function(alpha, W) {
  obj <- make_inner_obj(W)
  opt <- optim(alpha, obj, method = "BFGS", control = list(maxit = 1000))
  opt$par
}
res2 <- run_iterative_cue(inner_bfgs)
cat(sprintf("   iter=%d, obj=%.6e, max|G|=%.4e, alpha=(%s)\n\n",
    res2$iter, res2$obj, res2$grad_norm, paste(round(res2$alpha, 4), collapse=", ")))

# 3. Nelder-Mead
cat("--- 3. Nelder-Mead ---\n")
inner_nm <- function(alpha, W) {
  obj <- make_inner_obj(W)
  opt <- optim(alpha, obj, method = "Nelder-Mead", control = list(maxit = 5000))
  opt$par
}
res3 <- run_iterative_cue(inner_nm)
cat(sprintf("   iter=%d, obj=%.6e, max|G|=%.4e, alpha=(%s)\n\n",
    res3$iter, res3$obj, res3$grad_norm, paste(round(res3$alpha, 4), collapse=", ")))

# 4. CG (conjugate gradient)
cat("--- 4. CG ---\n")
inner_cg <- function(alpha, W) {
  obj <- make_inner_obj(W)
  opt <- optim(alpha, obj, method = "CG", control = list(maxit = 1000))
  opt$par
}
res4 <- run_iterative_cue(inner_cg)
cat(sprintf("   iter=%d, obj=%.6e, max|G|=%.4e, alpha=(%s)\n\n",
    res4$iter, res4$obj, res4$grad_norm, paste(round(res4$alpha, 4), collapse=", ")))

# 5. nlminb
cat("--- 5. nlminb ---\n")
inner_nlminb <- function(alpha, W) {
  obj <- make_inner_obj(W)
  res <- nlminb(alpha, obj, control = list(iter.max = 1000))
  res$par
}
res5 <- run_iterative_cue(inner_nlminb)
cat(sprintf("   iter=%d, obj=%.6e, max|G|=%.4e, alpha=(%s)\n\n",
    res5$iter, res5$obj, res5$grad_norm, paste(round(res5$alpha, 4), collapse=", ")))

# 6. L-BFGS-B with analytical gradient
cat("--- 6. L-BFGS-B with analytical gradient ---\n")
inner_lbfgsb_grad <- function(alpha, W) {
  obj <- make_inner_obj(W)
  gr <- function(a) {
    g_mat <- g_fn(a)
    G <- colMeans(g_mat)
    # Gamma = E[dg/dalpha]
    pi_vec <- compute_pi(pmin(pmax(as.vector(design_mat %*% a), -20), 20))
    common <- -r_vec / (pi_vec^2) * (-pi_vec * (1 - pi_vec))
    Gamma_mat <- matrix(0, h_dim, alpha_dim)
    for (j in 1:alpha_dim) {
      Gamma_mat[, j] <- colMeans(common * design_mat[, j] * h_x)
    }
    2 * as.vector(t(Gamma_mat) %*% W %*% G)
  }
  opt <- optim(alpha, obj, gr, method = "L-BFGS-B", control = list(maxit = 1000))
  opt$par
}
res6 <- run_iterative_cue(inner_lbfgsb_grad)
cat(sprintf("   iter=%d, obj=%.6e, max|G|=%.4e, alpha=(%s)\n\n",
    res6$iter, res6$obj, res6$grad_norm, paste(round(res6$alpha, 4), collapse=", ")))

# 7. SANN (simulated annealing)
cat("--- 7. SANN (simulated annealing inner) ---\n")
inner_sann <- function(alpha, W) {
  obj <- make_inner_obj(W)
  opt <- optim(alpha, obj, method = "SANN", control = list(maxit = 5000, temp = 1))
  opt$par
}
res7 <- run_iterative_cue(inner_sann, max_outer = 30)
cat(sprintf("   iter=%d, obj=%.6e, max|G|=%.4e, alpha=(%s)\n\n",
    res7$iter, res7$obj, res7$grad_norm, paste(round(res7$alpha, 4), collapse=", ")))
