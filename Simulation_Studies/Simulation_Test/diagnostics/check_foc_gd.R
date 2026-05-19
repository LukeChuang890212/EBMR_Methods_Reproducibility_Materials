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

# Generate datasets
set.seed(42)
n_val <- 2000
test_reps <- c(1, 5, 10, 50, 92, 100, 119, 139, 150)
max_rep <- max(test_reps)
datasets <- vector("list", max_rep)
for (i in 1:max_rep) datasets[[i]] <- setting3.B1(n_val)

# CU-GMM objective Q(alpha) = G(alpha)'W(alpha)G(alpha)
Q_cugmm <- function(alpha, g_func, n, esteq_dim) {
  g_mat <- g_func(alpha)
  G_hat <- matrix(.colMeans(g_mat, n, esteq_dim), esteq_dim, 1)
  W_hat <- tryCatch(solve(t(g_mat) %*% g_mat / n), error = function(e) diag(esteq_dim))
  as.numeric(crossprod(G_hat, W_hat %*% G_hat))
}

methods <- c("iterated", "gd")

cat("======================================================================\n")
cat("Checking FOC of Q(alpha) = G(alpha)'W(alpha)G(alpha) at solutions\n")
cat("FOC: dQ/dalpha = 0  (numerical gradient via numDeriv)\n")
cat("======================================================================\n\n")

for (rep_i in test_reps) {
  dat <- datasets[[rep_i]]
  n <- nrow(dat)
  r <- dat$r
  y <- dat[["y"]]

  # Build g function
  terms_obj <- terms(FORMULAS$u1_z2, data = dat)
  X <- model.matrix(terms_obj, data = dat)
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

  cat(sprintf("=== Rep %d ===\n", rep_i))

  for (m in methods) {
    tryCatch({
      ebmr <- EBMRAlgorithmFast4$new("y", ps_spec_m2, dat, W_func, bounds = 5, gmm_method = m)
      fit <- ebmr$ps_fit.list[[1]]
      alpha_hat <- fit$coefficients
      gmm_fit <- fit$gmm_fit

      # Numerical gradient of CU-GMM objective at solution
      grad_Q <- grad(function(a) Q_cugmm(a, g_func, n, h_dim), alpha_hat)

      # Q value
      g_mat <- g_func(alpha_hat)
      G_hat <- .colMeans(g_mat, n, h_dim)
      W_hat <- solve(t(g_mat) %*% g_mat / n)
      Q_val <- as.numeric(crossprod(G_hat, W_hat %*% G_hat))

      # Decompose: term1 = 2*Gamma'WG, term2 = G'W'G
      Gamma_hat <- jacobian(function(a) {
        gm <- g_func(a)
        .colMeans(gm, n, h_dim)
      }, alpha_hat)
      term1 <- 2 * as.vector(t(Gamma_hat) %*% W_hat %*% G_hat)
      term2 <- grad_Q - term1

      ps_vals <- inv_link_fn(as.vector(X %*% alpha_hat))
      mu_hat <- mean(r / ps_vals * y)

      cat(sprintf("  [%s] alpha=[%s]\n", toupper(m), paste(round(alpha_hat, 4), collapse=", ")))
      cat(sprintf("         Q=%.8f, mu=%.4f, conv=%s, iter=%d\n",
                  Q_val, mu_hat, gmm_fit$opt$converged, gmm_fit$opt$iterations))
      cat(sprintf("         max|dQ/dalpha|=%.6f\n", max(abs(grad_Q))))
      cat(sprintf("         max|term1(Gamma'WG)|=%.6f, max|term2(G'W'G)|=%.6f\n",
                  max(abs(term1)), max(abs(term2))))
      cat(sprintf("         dQ/dalpha=[%s]\n", paste(round(grad_Q, 6), collapse=", ")))
    }, error = function(e) {
      cat(sprintf("  [%s] ERROR: %s\n", toupper(m), conditionMessage(e)))
    })
  }
  cat("\n")
}

cat("Done!\n")
