setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
source("Data_Generation.r")
source("config/scenarios.R")

old_wd <- getwd()
setwd("../EBMRalgorithmFast4")
source("R/EBMRAlgorithm.r")
setwd(old_wd)

library(numDeriv)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
inv_link_fn <- function(eta) 1 / (1 + exp(eta))

ps_spec_m2 <- list(
  formula.list = list(FORMULAS$u1_z2),
  h_alpha.list = list(H_ALPHA$full),
  inv_link = inv_link_fn,
  outcome = "y",
  alpha_init.list = list(NULL)
)

set.seed(999)
big <- setting3.A1(1e6)
mu_true <- mean(big$y)
rm(big); gc()
cat(sprintf("mu_true = %.6f\n\n", mu_true))

# Test on a few representative reps
set.seed(42)
n_val <- 2000
test_reps <- c(1, 5, 10, 50, 92, 119, 139)  # include the 3 baseline outlier reps

# Generate datasets up to max(test_reps)
max_rep <- max(test_reps)
datasets <- vector("list", max_rep)
for (i in 1:max_rep) {
  datasets[[i]] <- setting3.B1(n_val)
}

# CU-GMM objective Q(alpha) = G(alpha)'W(alpha)G(alpha)
Q_cugmm <- function(alpha, g_func, n, esteq_dim) {
  g_mat <- g_func(alpha)
  G_hat <- matrix(.colMeans(g_mat, n, esteq_dim), esteq_dim, 1)
  W_hat <- tryCatch(solve(t(g_mat) %*% g_mat / n), error = function(e) diag(esteq_dim))
  as.numeric(crossprod(G_hat, W_hat %*% G_hat))
}

cat("======================================================================\n")
cat("Checking FOC of Q(alpha) = G(alpha)'W(alpha)G(alpha) at GMM solutions\n")
cat("FOC: dQ/dalpha = 0  (numerical gradient via numDeriv)\n")
cat("======================================================================\n\n")

for (rep_i in test_reps) {
  dat <- datasets[[rep_i]]
  cat(sprintf("--- Rep %d ---\n", rep_i))

  # Fit with baseline (no W' gradient term)
  # Temporarily remove dg to force numerical gradient in optim
  # Actually, we just need to compare. Let's fit with the current code (which has W' term)
  ebmr <- EBMRAlgorithmFast4$new("y", ps_spec_m2, dat, W_func, bounds = 5)
  fit <- ebmr$ps_fit.list[[1]]
  alpha_with_Wprime <- fit$coefficients
  gmm_fit <- fit$gmm_fit

  # Get the g function for this dataset
  # We need to reconstruct the g function. Let's use the internal gmm_fit info.
  # Instead, let's build it from the EBMRAlgorithm internals.
  # The g function is Phi_alpha from WangShaoKim2014.

  # Extract what we need from the fit object
  n <- nrow(dat)
  esteq_dim <- ncol(gmm_fit$g.matrix)

  # Build g function manually
  formula <- FORMULAS$u1_z2
  h_alpha_names <- H_ALPHA$full
  outcome <- "y"
  r <- dat$r
  y <- dat[[outcome]]

  # Parse formula to get design matrix
  terms_obj <- terms(formula, data = dat)
  X <- model.matrix(terms_obj, data = dat)
  alpha_dim <- ncol(X)

  # h_alpha matrix
  h_x <- as.matrix(dat[, h_alpha_names, drop = FALSE])
  h_dim <- ncol(h_x)

  g_func <- function(alpha) {
    eta <- as.vector(X %*% alpha)
    pi_hat <- inv_link_fn(eta)
    g_mat <- matrix(0, n, h_dim)
    for (l in 1:h_dim) {
      g_mat[, l] <- (r - pi_hat) * h_x[, l]
    }
    g_mat
  }

  # Compute numerical gradient of Q at the solution
  grad_Q <- grad(function(a) Q_cugmm(a, g_func, n, h_dim), alpha_with_Wprime)

  # Also compute the two components separately
  g_mat <- g_func(alpha_with_Wprime)
  G_hat <- .colMeans(g_mat, n, h_dim)
  W_hat <- solve(t(g_mat) %*% g_mat / n)
  Q_val <- as.numeric(crossprod(G_hat, W_hat %*% G_hat))

  # Compute Gamma(alpha) numerically
  Gamma_hat <- jacobian(function(a) {
    gm <- g_func(a)
    .colMeans(gm, n, h_dim)
  }, alpha_with_Wprime)

  # Term 1: 2 Gamma' W G
  term1 <- 2 * as.vector(t(Gamma_hat) %*% W_hat %*% G_hat)

  # Term 2: G' W'(alpha) G (from numerical gradient minus term1)
  term2 <- grad_Q - term1

  mu_hat <- mean(r / inv_link_fn(X %*% alpha_with_Wprime) * y)

  cat(sprintf("  alpha = [%s]\n", paste(round(alpha_with_Wprime, 4), collapse=", ")))
  cat(sprintf("  Q(alpha) = %.8f\n", Q_val))
  cat(sprintf("  mu_hat = %.4f\n", mu_hat))
  cat(sprintf("  ||G|| = %.6f\n", sqrt(sum(G_hat^2))))
  cat(sprintf("  Numerical dQ/dalpha = [%s]\n", paste(round(grad_Q, 6), collapse=", ")))
  cat(sprintf("  max|dQ/dalpha| = %.6f\n", max(abs(grad_Q))))
  cat(sprintf("  Term1 (2*Gamma'WG) = [%s]\n", paste(round(term1, 6), collapse=", ")))
  cat(sprintf("  Term2 (G'W'G)      = [%s]\n", paste(round(term2, 6), collapse=", ")))
  cat(sprintf("  max|Term1| = %.6f,  max|Term2| = %.6f\n", max(abs(term1)), max(abs(term2))))
  cat(sprintf("  |Term2|/|Term1| ratio = %.4f\n\n",
              max(abs(term2)) / max(abs(term1) + 1e-16)))
}

# Now also check what the baseline (no W' term) solutions look like
cat("\n======================================================================\n")
cat("Now checking baseline solutions (numerical grad, no W' term)\n")
cat("We re-fit without the analytical gradient\n")
cat("======================================================================\n\n")

# To get baseline solutions, we need to temporarily disable the analytical gradient.
# The simplest way: modify gmm to not use grad_obj by setting dg = NULL effectively.
# Instead, let's just source a version without the W' term.
# Actually, simpler: let's just check the iterated GMM FOC (Gamma'W_final G = 0)
# AND the CU-GMM FOC at the baseline solution.

# Reload without W' term: use the old GMM by setting dg=NULL in the algorithm
# Actually we can't easily do that. Let's just numerically check the FOC
# at the iterated GMM solution from the previous diagnostic run.

# Re-run with a modified version that skips the W' term
# Simplest: create a temporary gmm without analytical gradient
gmm_no_Wprime <- function(g_func_local, n_local, h_dim_local, alpha_dim_local, bounds_local) {
  W_local <- function(gm) solve(t(gm) %*% gm / n_local)

  G_local <- function(param) {
    gm <- g_func_local(param)
    matrix(.colMeans(gm, n_local, h_dim_local), h_dim_local, 1)
  }

  obj_local <- function(param) {
    gm <- g_func_local(param)
    G_hat <- matrix(.colMeans(gm, n_local, h_dim_local), h_dim_local, 1)
    W_hat <- state_local$W.hat
    if (is.null(W_hat)) {
      as.numeric(crossprod(G_hat))
    } else {
      as.numeric(crossprod(G_hat, W_hat %*% G_hat))
    }
  }

  state_local <- new.env(parent = emptyenv())
  state_local$W.hat <- NULL

  init <- rep(0, alpha_dim_local)
  sol <- init
  conv_err <- 1e8
  prev_conv_err <- 1e8

  for (t in 1:1000) {
    opt <- optim(sol, obj_local, method = "L-BFGS-B",
                 lower = rep(-bounds_local, alpha_dim_local),
                 upper = rep(bounds_local, alpha_dim_local),
                 control = list(maxit = 1000))
    new_sol <- opt$par
    gm <- g_func_local(new_sol)
    state_local$W.hat <- W_local(gm)
    conv_err <- max(abs(new_sol - sol))
    sol <- new_sol

    if (conv_err < 1e-6) break
    if (t > 3) {
      rel_change <- abs(conv_err - prev_conv_err) / (prev_conv_err + 1e-16)
      if (rel_change < 0.01) break
    }
    prev_conv_err <- conv_err
  }
  list(estimates = sol, iterations = t, converged = conv_err < 1e-6)
}

for (rep_i in test_reps) {
  dat <- datasets[[rep_i]]
  cat(sprintf("--- Rep %d (baseline, no W' term) ---\n", rep_i))

  n <- nrow(dat)
  r <- dat$r
  y <- dat[["y"]]

  terms_obj <- terms(FORMULAS$u1_z2, data = dat)
  X <- model.matrix(terms_obj, data = dat)
  alpha_dim <- ncol(X)

  h_x <- as.matrix(dat[, H_ALPHA$full, drop = FALSE])
  h_dim <- ncol(h_x)

  g_func <- function(alpha) {
    eta <- as.vector(X %*% alpha)
    pi_hat <- inv_link_fn(eta)
    g_mat <- matrix(0, n, h_dim)
    for (l in 1:h_dim) {
      g_mat[, l] <- (r - pi_hat) * h_x[, l]
    }
    g_mat
  }

  fit_base <- gmm_no_Wprime(g_func, n, h_dim, alpha_dim, 5)
  alpha_base <- fit_base$estimates

  # Check CU-GMM FOC at baseline solution
  grad_Q <- grad(function(a) Q_cugmm(a, g_func, n, h_dim), alpha_base)

  g_mat <- g_func(alpha_base)
  G_hat <- .colMeans(g_mat, n, h_dim)
  W_hat <- solve(t(g_mat) %*% g_mat / n)
  Q_val <- as.numeric(crossprod(G_hat, W_hat %*% G_hat))

  Gamma_hat <- jacobian(function(a) {
    gm <- g_func(a)
    .colMeans(gm, n, h_dim)
  }, alpha_base)

  term1 <- 2 * as.vector(t(Gamma_hat) %*% W_hat %*% G_hat)
  term2 <- grad_Q - term1

  mu_hat <- mean(r / inv_link_fn(X %*% alpha_base) * y)

  cat(sprintf("  alpha = [%s]\n", paste(round(alpha_base, 4), collapse=", ")))
  cat(sprintf("  Q(alpha) = %.8f\n", Q_val))
  cat(sprintf("  mu_hat = %.4f\n", mu_hat))
  cat(sprintf("  conv=%s, iter=%d\n", fit_base$converged, fit_base$iterations))
  cat(sprintf("  ||G|| = %.6f\n", sqrt(sum(G_hat^2))))
  cat(sprintf("  Numerical dQ/dalpha = [%s]\n", paste(round(grad_Q, 6), collapse=", ")))
  cat(sprintf("  max|dQ/dalpha| = %.6f\n", max(abs(grad_Q))))
  cat(sprintf("  Term1 (2*Gamma'WG) = [%s]\n", paste(round(term1, 6), collapse=", ")))
  cat(sprintf("  Term2 (G'W'G)      = [%s]\n", paste(round(term2, 6), collapse=", ")))
  cat(sprintf("  max|Term1| = %.6f,  max|Term2| = %.6f\n", max(abs(term1)), max(abs(term2))))
  cat(sprintf("  |Term2|/|Term1| ratio = %.4f\n\n",
              max(abs(term2)) / (max(abs(term1)) + 1e-16)))
}

cat("Done!\n")
