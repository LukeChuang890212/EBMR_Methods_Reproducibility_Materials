setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Data_Generation.r")
source("config/scenarios.R")
library(EBMRalgorithmFast4)
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
h_nu_func <- function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2, u1_u2 = dat$u1*dat$u2)

rep_i <- 27
dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]
ebmr <- EBMRAlgorithmFast4$new("y", subset_ps_spec, dat, W_func)

J <- 2
ps.matrix <- do.call(cbind, lapply(ebmr$ps_fit.list, function(f) f$fitted.values))
r_vec <- as.numeric(dat$r)
n <- nrow(dat)

h_nu_vals <- h_nu_func(dat)

# PACKAGE's h_x: intercept + continuous variables
# (since all h_nu vars are continuous, separate_variable_types puts them in x2)
h_x_pkg <- cbind(1, h_nu_vals)  # intercept + u1, u2, z1, z2, u1*u2
cat("Package h_x: intercept + 5 covariates =", ncol(h_x_pkg), "columns\n")

# MANUAL h_x: ps differences + continuous variables
d_manual <- ps.matrix[, -1, drop = FALSE] - ps.matrix[, 1]
h_x_manual <- cbind(d_manual, h_nu_vals)
cat("Manual h_x: ps_diff + 5 covariates =", ncol(h_x_manual), "columns\n\n")

run_gmm <- function(h_x, label) {
  h_dim <- ncol(h_x)

  Phi_nu <- function(nu) {
    ps_nu <- as.vector(ps.matrix %*% nu)
    (r_vec / ps_nu - 1) * h_x
  }

  cat(sprintf("=== %s (h_dim=%d) ===\n", label, h_dim))

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
      if (conv_err < 1e-6) break
    }
    w <- estimates^2 / sum(estimates^2)
    cat(sprintf("  init=%-8s: nu=[%7.4f,%7.4f], w=[%.4f,%.4f], obj=%.6e, iter=%d\n",
                nm, estimates[1], estimates[2], w[1], w[2], opt$value, t))
  }
  cat("\n")
}

run_gmm(h_x_pkg, "PACKAGE h_x (intercept)")
run_gmm(h_x_manual, "MANUAL h_x (ps_diff)")

cat("Done!\n")
