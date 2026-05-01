setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Data_Generation.r")
source("config/scenarios.R")
library(EBMRalgorithmFast4)
library(numDeriv)
source("Basic_setup.r")

data_file <- misspecified_model_all_data_file.list$setting4$miss50[[1]]
all_data <- readRDS(data_file)
n_val <- 2000

ps_spec <- get_ps_spec("9-alt1")
subset_ps_spec <- list(
  formula.list = ps_spec$formula.list[c(1, 3)],
  h_alpha.list = ps_spec$h_alpha.list[c(1, 3)],
  inv_link = ps_spec$inv_link,
  outcome = ps_spec$outcome
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

rep_i <- 27
dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]
ebmr <- EBMRAlgorithmFast4$new("y", subset_ps_spec, dat, W_func)

J <- 2
ps.matrix <- do.call(cbind, lapply(ebmr$ps_fit.list, function(f) f$fitted.values))
h_x2 <- cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2, u1_u2 = dat$u1*dat$u2)
d <- ps.matrix[, -1, drop = FALSE] - ps.matrix[, 1]
h_x <- cbind(d, h_x2)
h_dim <- ncol(h_x)
r_vec <- as.numeric(dat$r)
n <- nrow(dat)

Phi_nu <- function(nu) {
  ps_nu <- as.vector(ps.matrix %*% nu)
  (r_vec / ps_nu - 1) * h_x
}

gmm_obj_with_W <- function(nu, W_mat) {
  g_mat <- Phi_nu(nu)
  G_val <- matrix(colMeans(g_mat), h_dim, 1)
  as.numeric(crossprod(G_val, W_mat %*% G_val))
}

cat("=== PACKAGE-STYLE GMM (W=I init, then iterative) ===\n\n")

for (nm in c("uniform", "model1", "model2")) {
  init <- switch(nm,
    "uniform" = c(0.5, 0.5),
    "model1"  = c(1, 0),
    "model2"  = c(0, 1)
  )

  # Step 1: W = I
  W_hat <- diag(h_dim)
  obj_fn <- function(nu) {
    g_m <- Phi_nu(nu)
    G_v <- matrix(colMeans(g_m), h_dim, 1)
    v <- as.numeric(crossprod(G_v, W_hat %*% G_v))
    if (!is.finite(v)) 1e8 else v
  }
  opt_init <- optim(init, obj_fn, method = "L-BFGS-B",
                    lower = rep(-Inf, J), upper = rep(Inf, J),
                    control = list(maxit = 1000))
  estimates <- opt_init$par
  cat(sprintf("init=%-8s Step1(W=I): nu=[%7.4f,%7.4f], obj=%.6e\n",
              nm, estimates[1], estimates[2], opt_init$value))

  # Step 2: Iterative
  for (t in 1:200) {
    g_mat <- Phi_nu(estimates)
    W_hat <- tryCatch(solve(t(g_mat) %*% g_mat / n), error = function(e) diag(h_dim))

    grad_val <- numDeriv::grad(obj_fn, estimates)
    grad_norm <- max(abs(grad_val))
    if (grad_norm < 1e-6) {
      cat(sprintf("  iter %3d: CONVERGED (grad_norm=%.2e)\n", t, grad_norm))
      break
    }

    opt_t <- optim(estimates, obj_fn, method = "L-BFGS-B",
                   lower = rep(-Inf, J), upper = rep(Inf, J),
                   control = list(maxit = 1000))
    conv_err <- max(abs(opt_t$par - estimates))
    estimates <- opt_t$par

    if (t <= 5 || t %% 50 == 0) {
      cat(sprintf("  iter %3d: nu=[%7.4f,%7.4f], obj=%.6e, grad=%.2e, conv=%.2e\n",
                  t, estimates[1], estimates[2], opt_t$value, grad_norm, conv_err))
    }

    if (conv_err < 1e-6) {
      cat(sprintf("  iter %3d: CONVERGED (conv_err=%.2e)\n", t, conv_err))
      break
    }
  }
  w <- estimates^2 / sum(estimates^2)
  cat(sprintf("  FINAL: nu=[%7.4f,%7.4f], w=[%.4f,%.4f]\n\n",
              estimates[1], estimates[2], w[1], w[2]))
}

cat("\n=== MANUAL GMM (direct iterative, no W=I init) ===\n\n")

for (nm in c("uniform", "model1", "model2")) {
  init <- switch(nm,
    "uniform" = c(0.5, 0.5),
    "model1"  = c(1, 0.001),
    "model2"  = c(0.001, 1)
  )
  estimates <- init

  for (t in 1:200) {
    g_mat <- Phi_nu(estimates)
    W_hat <- tryCatch(solve(t(g_mat) %*% g_mat / n), error = function(e) diag(h_dim))

    obj_fn <- function(nu) {
      g_m <- Phi_nu(nu)
      G_v <- matrix(colMeans(g_m), h_dim, 1)
      v <- as.numeric(crossprod(G_v, W_hat %*% G_v))
      if (!is.finite(v)) 1e8 else v
    }

    opt <- optim(estimates, obj_fn, method = "L-BFGS-B",
                 lower = rep(-Inf, J), upper = rep(Inf, J),
                 control = list(maxit = 1000))
    conv_err <- max(abs(opt$par - estimates))
    estimates <- opt$par

    if (t <= 5 || t %% 50 == 0) {
      cat(sprintf("  iter %3d: nu=[%7.4f,%7.4f], obj=%.6e, conv=%.2e\n",
                  t, estimates[1], estimates[2], opt$value, conv_err))
    }

    if (conv_err < 1e-6) {
      cat(sprintf("  iter %3d: CONVERGED (conv_err=%.2e)\n", t, conv_err))
      break
    }
  }
  w <- estimates^2 / sum(estimates^2)
  cat(sprintf("  init=%-8s FINAL: nu=[%7.4f,%7.4f], w=[%.4f,%.4f]\n\n",
              nm, estimates[1], estimates[2], w[1], w[2]))
}

cat("Done!\n")
