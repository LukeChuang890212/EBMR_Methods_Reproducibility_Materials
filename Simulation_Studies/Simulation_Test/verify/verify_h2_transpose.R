setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EPS", quiet = TRUE)

n_val <- 2000
ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)

W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

# Use model 3 (the problematic one), constrained_nr, rep 1
model_idx <- 3
dat <- all_data[1:n_val, ]

single_ps <- list(
  formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
  h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]],
  alpha_init.list = list(NULL),
  optimizer = "constrained_nr"
)

ebmr <- EPS$new("y", single_ps, dat, W_fn)
gmm_fit <- ebmr$ps_fit.list[[1]]$gmm_fit
design_mat <- ebmr$ps_fit.list[[1]]$design_matrix
link_type <- ebmr$ps_fit.list[[1]]$link_type
alpha_hat <- gmm_fit$estimates
n <- n_val
k <- length(alpha_hat)  # param_dim
r_vec <- dat[["r"]]; y_vec <- dat[["y"]]

cat(sprintf("Model 3, constrained_nr, n=%d, k=%d\n", n, k))
cat(sprintf("alpha_hat = %s\n", paste(round(alpha_hat, 4), collapse=", ")))
cat(sprintf("FOC norm = %.2e\n", max(abs(crossprod(gmm_fit$Gamma.hat, gmm_fit$W.hat) %*% gmm_fit$eta_s))))

# ============================================================
# Build the score function s(theta) = Gamma(theta)'W(theta)G(theta)
# where G = mean(g(theta)), W = (g'g/n)^{-1}, Gamma = mean(dg/dtheta)
# ============================================================

# Functions to evaluate at arbitrary alpha
h_alpha_vars <- ps_spec[["h_alpha.list"]][[model_idx]]
h_x <- cbind(1, as.matrix(dat[, h_alpha_vars, drop = FALSE]))
h_dim <- ncol(h_x)

compute_pi_fn <- function(eta) plogis(-eta)  # logistic_complement

g_fn <- function(alpha) {
  eta <- as.vector(design_mat %*% alpha)
  pi_v <- compute_pi_fn(eta)
  (r_vec / pi_v - 1) * h_x  # n x h
}

score_fn <- function(alpha) {
  g_mat <- g_fn(alpha)
  G_bar <- colMeans(g_mat)
  W_mat <- solve(crossprod(g_mat) / n)

  # Gamma = mean of individual Jacobians
  eta <- as.vector(design_mat %*% alpha)
  pi_v <- compute_pi_fn(eta)
  cf <- r_vec * (1 - pi_v) / pi_v  # common factor for logistic_complement
  # Gamma[l, j] = mean(cf * h_x[,l] * X[,j])
  Gamma_mat <- crossprod(h_x * cf, design_mat) / n  # h x k

  as.vector(crossprod(Gamma_mat, W_mat %*% G_bar))  # k x 1
}

# ============================================================
# Numerical H2: ds/dalpha' via central finite differences
# ============================================================
eps <- 1e-5
H2_num <- matrix(0, k, k)
for (j in 1:k) {
  ej <- rep(0, k); ej[j] <- eps
  H2_num[, j] <- (score_fn(alpha_hat + ej) - score_fn(alpha_hat - ej)) / (2 * eps)
}

# ============================================================
# Analytical H2 from the package
# ============================================================
Gamma_hat <- gmm_fit$Gamma.hat
g_mat <- gmm_fit$g.matrix
W_hat <- gmm_fit$W.hat
eta_s <- gmm_fit$eta_s
GtW <- crossprod(Gamma_hat, W_hat)
GtWG <- GtW %*% Gamma_hat

# Reconstruct H2 from package components
# Need R_mat and S_mat
# R_mat: use d2g_alpha
pi_hat <- ebmr$ps_fit.list[[1]]$fitted.values
eta_vals <- as.vector(design_mat %*% alpha_hat)
dc_deta <- r_vec * (1 - pi_hat) / pi_hat  # logistic_complement
dc_X <- dc_deta * design_mat  # n x k
R_mat <- matrix(0, h_dim * k, k)
for (j in 1:k) {
  Xj_h <- design_mat[, j] * h_x  # n x h
  block <- crossprod(dc_X, Xj_h) / n  # k x h
  for (l in 1:h_dim) {
    R_mat[(l-1)*k + j, ] <- block[, l]
  }
}

# S_mat: use dg_alpha
# dg_j[i, l] = dg_il/dalpha_j = cf_i * X[i,j] * h_x[i,l]
cf <- r_vec * (1 - pi_hat) / pi_hat
S_mat <- matrix(0, h_dim^2, k)
for (j in 1:k) {
  dg_j <- cf * design_mat[, j] * h_x  # n x h (dg_il/dalpha_j)
  cross <- (crossprod(dg_j, g_mat) + crossprod(g_mat, dg_j)) / n  # h x h
  S_mat[, j] <- as.vector(cross)
}

GW_row <- as.vector(t(eta_s) %*% W_hat)  # G'W (1 x h)
GW_kron_Ik <- kronecker(matrix(GW_row, 1, h_dim), diag(k))  # k x (h*k)
GW_kron_GtW <- kronecker(matrix(GW_row, 1, h_dim), GtW)     # k x (h*h)
H2_analytical <- GtWG + GW_kron_Ik %*% R_mat - GW_kron_GtW %*% S_mat

cat("\n=== H2 comparison ===\n")
cat(sprintf("||H2_num - H2_analytical|| = %.6e\n", norm(H2_num - H2_analytical, "F")))
cat(sprintf("||H2_num|| = %.6f\n", norm(H2_num, "F")))
cat(sprintf("||H2_analytical|| = %.6f\n", norm(H2_analytical, "F")))
cat(sprintf("Relative diff = %.6e\n", norm(H2_num - H2_analytical, "F") / norm(H2_num, "F")))

cat("\nH2_num:\n"); print(round(H2_num, 4))
cat("\nH2_analytical:\n"); print(round(H2_analytical, 4))
cat("\nDifference:\n"); print(round(H2_num - H2_analytical, 6))

# ============================================================
# Now check: what if we swap the Kronecker order?
# Try (I_k ⊗ G'W) R instead of (G'W ⊗ I_k) R
# Try (Γ'W ⊗ G'W) S instead of (G'W ⊗ Γ'W) S
# ============================================================

# Alt 1: swap both Kronecker orders
GW_kron_Ik_alt <- kronecker(diag(k), matrix(GW_row, 1, h_dim))  # k x (k*h)
GW_kron_GtW_alt <- kronecker(GtW, matrix(GW_row, 1, h_dim))     # k x (h*h)

# For alt Kronecker, we also need to reorder R_mat and S_mat
# If using (I_k ⊗ G'W), the vec convention changes to vec(Γ) instead of vec(Γ')
# vec(Γ) has Γ[l,j] at position l + (j-1)*h, i.e., column-major of Γ (h x k)
# This gives R_alt[l + (j-1)*h, m] = dΓ[l,j]/dtheta_m
R_alt <- matrix(0, h_dim * k, k)
for (j in 1:k) {
  Xj_h <- design_mat[, j] * h_x  # n x h
  block <- crossprod(dc_X, Xj_h) / n  # k x h
  for (l in 1:h_dim) {
    R_alt[l + (j-1)*h_dim, ] <- block[, l]
  }
}

H2_alt1 <- GtWG + GW_kron_Ik_alt %*% R_alt - GW_kron_GtW_alt %*% S_mat
cat("\n=== Alt H2 (swapped Kronecker, vec(Γ) convention) ===\n")
cat(sprintf("||H2_num - H2_alt1|| = %.6e\n", norm(H2_num - H2_alt1, "F")))
cat(sprintf("Relative diff = %.6e\n", norm(H2_num - H2_alt1, "F") / norm(H2_num, "F")))

# Alt 2: just transpose the S Kronecker
H2_alt2 <- GtWG + GW_kron_Ik %*% R_mat - GW_kron_GtW_alt %*% S_mat
cat("\n=== Alt H2 (only S Kronecker swapped) ===\n")
cat(sprintf("||H2_num - H2_alt2|| = %.6e\n", norm(H2_num - H2_alt2, "F")))

# Alt 3: just transpose the R Kronecker
H2_alt3 <- GtWG + GW_kron_Ik_alt %*% R_alt - GW_kron_GtW %*% S_mat
cat("\n=== Alt H2 (only R Kronecker swapped) ===\n")
cat(sprintf("||H2_num - H2_alt3|| = %.6e\n", norm(H2_num - H2_alt3, "F")))

# ============================================================
# Numerical Q_i verification
# Q_i should be the "score" contribution of observation i
# s(θ) = (1/n) Σ_i q_i(θ) where q_i captures obs i's contribution
# through g_i, Γ_i, and W^{-1}_{contribution_from_i}
# ============================================================
cat("\n=== Q_i comparison (first 3 obs) ===\n")

WG <- as.vector(W_hat %*% eta_s)
GtW_g <- GtW %*% t(g_mat)           # k x n
g_WG <- as.vector(g_mat %*% WG)      # n-vector

# Individual Gamma_i
cf_all <- r_vec * (1 - pi_hat) / pi_hat
# Gamma_i[l, j] = cf_i * h_x[i,l] * X[i,j]

# GammaT_WG: Gamma_i' W G for each i
# Gamma_i' is k x h, W is h x h, G is h x 1
# (Gamma_i' W G)[j] = sum_l Gamma_i[l,j] * (WG)_l = cf_i * X[i,j] * sum_l h_x[i,l] * WG[l]
hx_WG <- as.vector(h_x %*% WG)  # n-vector
cf_hxWG <- cf_all * hx_WG
GammaT_WG <- t(design_mat * cf_hxWG)  # k x n

Q_analytical <- GtW_g + GammaT_WG - sweep(GtW_g, 2, g_WG, `*`)

# Numerical Q_i: perturb observation weight
# If we up-weight obs i by epsilon, the score changes by epsilon * Q_i / n
# More precisely: s(θ; w) where w_i = 1 + ε, w_j = 1 for j≠i
# ds/dw_i ≈ Q_i / n (at w = 1)
# But this is expensive. Instead, verify via leave-one-out:

for (i_test in 1:3) {
  # Numerical: s with all obs vs s without obs i
  # s_full = (1/n) Σ Gamma_full' W_full G_full
  # s_{-i} ≈ ... complicated, just do finite diff on weight

  # Leave-one-out approach: compute everything without obs i
  idx <- setdiff(1:n, i_test)
  n_m1 <- n - 1
  g_loo <- g_fn(alpha_hat)[idx, ]
  G_loo <- colMeans(g_loo)
  W_loo <- solve(crossprod(g_loo) / n_m1)

  eta_loo <- as.vector(design_mat[idx, ] %*% alpha_hat)
  pi_loo <- compute_pi_fn(eta_loo)
  cf_loo <- r_vec[idx] * (1 - pi_loo) / pi_loo
  Gamma_loo <- crossprod(h_x[idx, ] * cf_loo, design_mat[idx, ]) / n_m1

  s_loo <- as.vector(crossprod(Gamma_loo, W_loo %*% G_loo))
  s_full <- score_fn(alpha_hat)

  # Approximate Q_i: n * (s_full - (n-1)/n * s_loo)
  # Note: s_full = (1/n)Σ q_i, s_loo = (1/(n-1))Σ_{j≠i} q_j
  # So n*s_full = Σ q_i, (n-1)*s_loo = Σ_{j≠i} q_j
  # q_i = n*s_full - (n-1)*s_loo
  Q_i_num <- n * s_full - (n-1) * s_loo
  Q_i_ana <- Q_analytical[, i_test]

  cat(sprintf("  Obs %d: Q_num = %s\n", i_test, paste(round(Q_i_num, 6), collapse=", ")))
  cat(sprintf("         Q_ana = %s\n", paste(round(Q_i_ana, 6), collapse=", ")))
  cat(sprintf("         diff  = %s\n", paste(round(Q_i_num - Q_i_ana, 6), collapse=", ")))
  cat(sprintf("         rel   = %.6e\n\n", sqrt(sum((Q_i_num - Q_i_ana)^2)) / max(sqrt(sum(Q_i_ana^2)), 1e-10)))
}
