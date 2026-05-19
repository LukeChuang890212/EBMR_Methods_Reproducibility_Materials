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

# Rep 27 (outlier)
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

# Phi_nu (same as package)
Phi_nu <- function(nu) {
  ps_nu <- as.vector(ps.matrix %*% nu)
  (r_vec / ps_nu - 1) * h_x
}

# Analytical dg_nu (same as package)
dg_nu_analytical <- function(param) {
  ps_nu <- as.vector(ps.matrix %*% param)
  neg_r_ps2 <- -r_vec / (ps_nu * ps_nu)
  base_mat <- neg_r_ps2 * ps.matrix
  Gamma_arr <- array(0, dim = c(n, J, h_dim))
  for (l in 1:h_dim) {
    Gamma_arr[, , l] <- base_mat * h_x[, l]
  }
  return(Gamma_arr)
}

# Numerical dg_nu
dg_nu_numerical <- function(param) {
  Gamma_arr <- array(NA, dim = c(n, J, h_dim))
  for (l in 1:h_dim) {
    Gamma_arr[, , l] <- jacobian(function(p) Phi_nu(p)[, l], param)
  }
  return(Gamma_arr)
}

# Test at several nu values
test_points <- list(
  "uniform" = c(0.5, 0.5),
  "model1"  = c(1, 0.01),
  "model2"  = c(0.01, 1),
  "pkg_sol" = c(0.0245, 0.9628),
  "correct" = c(1.0626, -0.0494)
)

cat("=== Comparing analytical vs numerical dg_nu ===\n\n")

for (nm in names(test_points)) {
  nu <- test_points[[nm]]

  dg_a <- dg_nu_analytical(nu)
  dg_n <- dg_nu_numerical(nu)

  max_diff <- max(abs(dg_a - dg_n))
  rel_diff <- max(abs(dg_a - dg_n) / (abs(dg_n) + 1e-10))

  cat(sprintf("nu=%-12s: max_abs_diff=%.2e, max_rel_diff=%.2e\n",
              nm, max_diff, rel_diff))

  if (max_diff > 1e-4) {
    cat("  WARNING: Large gradient discrepancy!\n")
    # Show where the differences are largest
    for (j in 1:J) {
      for (l in 1:h_dim) {
        diff_jl <- max(abs(dg_a[, j, l] - dg_n[, j, l]))
        if (diff_jl > 1e-4) {
          cat(sprintf("  j=%d, l=%d: max_diff=%.6f\n", j, l, diff_jl))
        }
      }
    }
  }
}

# Now check the GMM gradient (the full objective gradient, not just dg)
cat("\n=== Comparing GMM objective gradient ===\n\n")

# The GMM gradient is: 2 * t(Gamma) %*% W %*% G
# where Gamma = mean(dg), G = mean(g), W = (g'g/n)^{-1}

for (nm in names(test_points)) {
  nu <- test_points[[nm]]

  g_mat <- Phi_nu(nu)
  G_val <- matrix(colMeans(g_mat), h_dim, 1)
  W_hat <- tryCatch(solve(t(g_mat) %*% g_mat / n), error = function(e) diag(h_dim))

  # Analytical Gamma
  dg_a <- dg_nu_analytical(nu)
  Gamma_a <- matrix(0, h_dim, J)
  for (j in 1:J) Gamma_a[, j] <- colMeans(dg_a[, j, , drop = FALSE])
  grad_a <- 2 * as.vector(crossprod(Gamma_a, W_hat %*% G_val))

  # Numerical gradient of GMM objective
  obj_fn <- function(p) {
    gm <- Phi_nu(p)
    Gv <- matrix(colMeans(gm), h_dim, 1)
    as.numeric(crossprod(Gv, W_hat %*% Gv))
  }
  grad_n <- as.vector(numDeriv::grad(obj_fn, nu))

  cat(sprintf("nu=%-12s:\n", nm))
  cat(sprintf("  grad_analytical: [%s]\n", paste(sprintf("%10.6f", grad_a), collapse=", ")))
  cat(sprintf("  grad_numerical:  [%s]\n", paste(sprintf("%10.6f", grad_n), collapse=", ")))
  cat(sprintf("  diff:            [%s]\n", paste(sprintf("%10.6f", grad_a - grad_n), collapse=", ")))
}

cat("\nDone!\n")
