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
alpha_hat <- gmm_fit$estimates
pi_hat <- ebmr$ps_fit.list[[1]]$fitted.values
n <- n_val; k <- length(alpha_hat)
r_vec <- dat[["r"]]; y_vec <- dat[["y"]]

# USE THE PACKAGE'S h_x from ps_fit result!
h_x <- ebmr$ps_fit.list[[1]]$h_x
h <- ncol(h_x)
cat(sprintf("Using package h_x: %d x %d\n", nrow(h_x), h))
cat(sprintf("h_x colnames: %s\n", paste(colnames(h_x), collapse=", ")))

compute_pi_fn <- function(eta) plogis(-eta)

# Build g_fn using the PACKAGE's h_x
g_fn <- function(alpha) {
  eta <- as.vector(design_mat %*% alpha)
  pi_v <- compute_pi_fn(eta)
  (r_vec / pi_v - 1) * h_x
}

# Verify g matches
g_mine <- g_fn(alpha_hat)
g_pkg <- gmm_fit$g.matrix
cat(sprintf("g match: max diff = %.2e\n", max(abs(g_mine - g_pkg))))

# Score function
score_fn <- function(alpha) {
  g_mat <- g_fn(alpha)
  G_bar <- colMeans(g_mat)
  W_mat <- solve(crossprod(g_mat) / n)
  eta <- as.vector(design_mat %*% alpha)
  pi_v <- compute_pi_fn(eta)
  cf <- r_vec * (1 - pi_v) / pi_v
  Gamma_mat <- crossprod(h_x * cf, design_mat) / n  # h x k
  as.vector(crossprod(Gamma_mat, W_mat %*% G_bar))  # k x 1
}

# Numerical H2
eps <- 1e-5
H2_num <- matrix(0, k, k)
for (j in 1:k) {
  ej <- rep(0, k); ej[j] <- eps
  H2_num[, j] <- (score_fn(alpha_hat + ej) - score_fn(alpha_hat - ej)) / (2 * eps)
}

# Analytical H2
Gamma_hat <- gmm_fit$Gamma.hat
g_mat <- gmm_fit$g.matrix
W_hat <- gmm_fit$W.hat
eta_s <- gmm_fit$eta_s
GtW <- crossprod(Gamma_hat, W_hat)
GtWG <- GtW %*% Gamma_hat

cf <- r_vec * (1 - pi_hat) / pi_hat
dc_deta <- cf  # For logistic_complement, dc/deta = c
dc_X <- dc_deta * design_mat

# R_mat using correct h_x
R_mat <- matrix(0, h * k, k)
for (j in 1:k) {
  Xj_h <- design_mat[, j] * h_x
  block <- crossprod(dc_X, Xj_h) / n
  for (l in 1:h) R_mat[(l-1)*k + j, ] <- block[, l]
}

# S_mat using correct h_x and g
S_mat <- matrix(0, h^2, k)
for (j in 1:k) {
  dg_j <- cf * design_mat[, j] * h_x
  cross <- (crossprod(dg_j, g_mat) + crossprod(g_mat, dg_j)) / n
  S_mat[, j] <- as.vector(cross)
}

# Verify R_mat and S_mat numerically
cat("\n=== R_mat verification ===\n")
Gamma_at <- function(alpha) {
  eta <- as.vector(design_mat %*% alpha)
  pi_v <- compute_pi_fn(eta)
  cf_l <- r_vec * (1 - pi_v) / pi_v
  crossprod(h_x * cf_l, design_mat) / n
}
R_num <- matrix(0, h * k, k)
for (j in 1:k) {
  ej <- rep(0, k); ej[j] <- eps
  R_num[, j] <- as.vector(t((Gamma_at(alpha_hat + ej) - Gamma_at(alpha_hat - ej)) / (2*eps)))
}
cat(sprintf("R: ||ana - num|| / ||num|| = %.2e\n",
    norm(R_mat - R_num, "F") / norm(R_num, "F")))

cat("\n=== S_mat verification ===\n")
S_num <- matrix(0, h^2, k)
for (j in 1:k) {
  ej <- rep(0, k); ej[j] <- eps
  g_p <- g_fn(alpha_hat + ej); g_m <- g_fn(alpha_hat - ej)
  S_num[, j] <- as.vector((crossprod(g_p) - crossprod(g_m)) / (2*n*eps))
}
cat(sprintf("S: ||ana - num|| / ||num|| = %.2e\n",
    norm(S_mat - S_num, "F") / norm(S_num, "F")))

# Build H2
GW_row <- as.vector(t(eta_s) %*% W_hat)
GW_kron_Ik <- kronecker(matrix(GW_row, 1, h), diag(k))
GW_kron_GtW <- kronecker(matrix(GW_row, 1, h), GtW)
H2_analytical <- GtWG + GW_kron_Ik %*% R_mat - GW_kron_GtW %*% S_mat

cat("\n=== H2 comparison ===\n")
cat(sprintf("||H2_num - H2_analytical|| / ||H2_num|| = %.6e\n",
    norm(H2_num - H2_analytical, "F") / norm(H2_num, "F")))

cat("\nH2_num:\n"); print(round(H2_num, 6))
cat("\nH2_analytical:\n"); print(round(H2_analytical, 6))
cat("\nDifference:\n"); print(round(H2_num - H2_analytical, 8))

# Also compare with numerical R and S
H2_with_num_RS <- GtWG + GW_kron_Ik %*% R_num - GW_kron_GtW %*% S_num
cat(sprintf("\n||H2_num - H2_with_num_RS|| / ||H2_num|| = %.6e\n",
    norm(H2_num - H2_with_num_RS, "F") / norm(H2_num, "F")))

# ============================================================
# Now: numerical Q_i verification with correct h_x
# ============================================================
cat("\n=== Q_i verification (correct h_x) ===\n")
WG <- as.vector(W_hat %*% eta_s)
GtW_g <- GtW %*% t(g_mat)
g_WG <- as.vector(g_mat %*% WG)

# Gamma_i'WG for each i
cf_hxWG <- cf * as.vector(h_x %*% WG)
GammaT_WG <- t(design_mat * cf_hxWG)  # k x n

Q_analytical <- GtW_g + GammaT_WG - sweep(GtW_g, 2, g_WG, `*`)

# Numerical Q_i via leave-one-out
s_full <- score_fn(alpha_hat)
for (i_test in c(1, 10, 100, 500)) {
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

  Q_i_num <- n * s_full - (n-1) * s_loo
  Q_i_ana <- Q_analytical[, i_test]
  rel <- sqrt(sum((Q_i_num - Q_i_ana)^2)) / sqrt(sum(Q_i_ana^2))
  cat(sprintf("  Obs %4d: rel_diff=%.6e, Q_num=%s\n", i_test, rel,
      paste(round(Q_i_num, 4), collapse=", ")))
  cat(sprintf("            %s Q_ana=%s\n", paste(rep(" ", nchar(sprintf("Obs %4d", i_test))), collapse=""),
      paste(round(Q_i_ana, 4), collapse=", ")))
}

# ============================================================
# Compare SE2 from analytical H2 vs numerical H2
# ============================================================
cat("\n=== SE2 from different H2 ===\n")

# H_alpha for mu
if (ebmr$ps_fit.list[[1]]$link_type == "logistic_complement") {
  dot_pi <- -design_mat * (pi_hat * (1 - pi_hat))
} else {
  dot_pi <- design_mat * (pi_hat * (1 - pi_hat))
}
ry_ps_inv2 <- as.vector(r_vec * y_vec * (pi_hat^(-2)))
H_alpha <- colMeans(dot_pi * ry_ps_inv2)

# With analytical H2
H2_inv_ana <- solve(H2_analytical)
psi2_ana <- -(H2_inv_ana %*% Q_analytical)
mu_iid_ana <- as.vector(t(r_vec/pi_hat*y_vec) - t(H_alpha) %*% psi2_ana)
se_ana <- sqrt(var(mu_iid_ana)/n)

# With numerical H2
H2_inv_num <- solve(H2_num)
psi2_num <- -(H2_inv_num %*% Q_analytical)
mu_iid_num <- as.vector(t(r_vec/pi_hat*y_vec) - t(H_alpha) %*% psi2_num)
se_num <- sqrt(var(mu_iid_num)/n)

# With package psi2 (should match analytical)
psi2_pkg <- gmm_fit$psi2
mu_iid_pkg <- as.vector(t(r_vec/pi_hat*y_vec) - t(H_alpha) %*% psi2_pkg)
se_pkg <- sqrt(var(mu_iid_pkg)/n)

cat(sprintf("  SE (analytical H2): %.6f\n", se_ana))
cat(sprintf("  SE (numerical H2):  %.6f\n", se_num))
cat(sprintf("  SE (package psi2):  %.6f\n", se_pkg))

# Also try: numerical H2 + numerical Q
cat("\n=== Fully numerical influence function ===\n")
# Numerical Q via leave-one-out for ALL obs (expensive for n=2000, sample 200)
set.seed(42)
sample_idx <- sort(sample(1:n, 200))
Q_num_sample <- matrix(0, k, length(sample_idx))
for (ii in seq_along(sample_idx)) {
  i_test <- sample_idx[ii]
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
  Q_num_sample[, ii] <- n * s_full - (n-1) * s_loo
}

# Compare Q variances
Q_ana_sample <- Q_analytical[, sample_idx]
cat(sprintf("  Var of Q (analytical): %s\n", paste(round(apply(Q_ana_sample, 1, var), 4), collapse=", ")))
cat(sprintf("  Var of Q (numerical):  %s\n", paste(round(apply(Q_num_sample, 1, var), 4), collapse=", ")))

# SE from fully numerical
psi2_full_num <- -(H2_inv_num %*% Q_num_sample)
psi2_full_ana <- -(H2_inv_ana %*% Q_ana_sample)
mu_iid_full_num <- as.vector(t(r_vec[sample_idx]/pi_hat[sample_idx]*y_vec[sample_idx]) - t(H_alpha) %*% psi2_full_num)
mu_iid_full_ana <- as.vector(t(r_vec[sample_idx]/pi_hat[sample_idx]*y_vec[sample_idx]) - t(H_alpha) %*% psi2_full_ana)
cat(sprintf("  Var(mu_iid) numerical: %.6f\n", var(mu_iid_full_num)))
cat(sprintf("  Var(mu_iid) analytical: %.6f\n", var(mu_iid_full_ana)))
