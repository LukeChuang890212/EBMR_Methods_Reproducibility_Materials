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
model_idx <- 3

single_ps <- list(
  formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
  h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]],
  alpha_init.list = list(NULL),
  optimizer = "constrained_nr"
)
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

# Fit one rep
dat <- all_data[1:n_val, ]
ebmr <- EPS$new("y", single_ps, dat, W_fn)
gmm_fit <- ebmr$ps_fit.list[[1]]$gmm_fit

alpha_hat <- gmm_fit$estimates
Gamma_hat <- gmm_fit$Gamma.hat
g_mat <- gmm_fit$g.matrix
W_hat <- gmm_fit$W.hat
eta_s <- gmm_fit$eta_s
n <- n_val
k <- length(alpha_hat)  # param_dim
h <- nrow(Gamma_hat)    # esteq_dim

cat(sprintf("alpha_hat: %s\n", paste(sprintf("%.6f", alpha_hat), collapse=", ")))
cat(sprintf("||eta_s||: %.6e\n", sqrt(sum(eta_s^2))))
cat(sprintf("k=%d, h=%d\n\n", k, h))

# ==========================================
# 1. Verify H2_mat numerically
# ==========================================
# H2_mat should equal d/dtheta' [Gamma(theta)' W(theta) G(theta)]
# We can check by finite differences of the FOC: foc(theta) = Gamma'WG

# Need the g function and Gamma_fast function
# Extract from the package internals
design_mat <- ebmr$ps_fit.list[[1]]$design_matrix
link_type <- ebmr$ps_fit.list[[1]]$link_type
r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
h_alpha_vars <- ps_spec[["h_alpha.list"]][[model_idx]]
h_x <- cbind(1, as.matrix(dat[, h_alpha_vars, drop = FALSE]))

inv_link <- ps_spec[["inv_link"]]

# Reconstruct g and Gamma_fast functions
get_pi <- function(alpha) {
  eta <- as.vector(design_mat %*% alpha)
  inv_link(eta)
}

g_fn <- function(alpha) {
  pi_vec <- get_pi(alpha)
  rw <- r_vec / pi_vec
  (rw - 1) * h_x
}

Gamma_fast_fn <- function(alpha) {
  pi_vec <- get_pi(alpha)
  if (link_type == "logistic_complement") {
    cf <- r_vec * (1 - pi_vec) / pi_vec
  } else {
    cf <- -r_vec * (1 - pi_vec) / pi_vec
  }
  crossprod(h_x * cf, design_mat) / n
}

# FOC function: Gamma'WG
foc <- function(alpha) {
  gm <- g_fn(alpha)
  G <- colMeans(gm)
  Wh <- tryCatch(W_fn(gm), error = function(e) diag(h))
  Gam <- Gamma_fast_fn(alpha)
  as.vector(crossprod(Gam, Wh %*% G))
}

# Numerical Jacobian of FOC
eps <- 1e-5
H2_numerical <- matrix(0, k, k)
for (j in 1:k) {
  ej <- rep(0, k); ej[j] <- eps
  H2_numerical[, j] <- (foc(alpha_hat + ej) - foc(alpha_hat - ej)) / (2 * eps)
}

# Package H2_mat (reconstructed from the se2 internals)
# Need to recompute R_mat and S_mat
library(numDeriv)

# Get t_Gamma_arr
dg_fn <- function(alpha) {
  pi_vec <- get_pi(alpha)
  if (link_type == "logistic_complement") {
    cf <- r_vec * (1 - pi_vec) / pi_vec
  } else {
    cf <- -r_vec * (1 - pi_vec) / pi_vec
  }
  Gamma_arr <- array(0, dim = c(n, k, h))
  for (l in 1:h) {
    Gamma_arr[, , l] <- (cf * h_x[, l]) * design_mat
  }
  Gamma_arr
}

t_Gamma_arr <- dg_fn(alpha_hat)

# R_mat from analytical d2g
if (link_type == "logistic_complement" || link_type == "logistic") {
  pi_vec <- get_pi(alpha_hat)
  dc_deta <- r_vec * (1 - pi_vec) / pi_vec
}
dc_X <- dc_deta * design_mat
R_mat_analytical <- matrix(0, h * k, k)
for (j in 1:k) {
  Xj_h <- design_mat[, j] * h_x
  block <- crossprod(dc_X, Xj_h) / n
  for (l in 1:h) {
    R_mat_analytical[(l-1)*k + j, ] <- block[, l]
  }
}

# R_mat numerical
R_mat_numerical <- matrix(0, k * h, k)
for (j in 1:k) {
  ej <- rep(0, k); ej[j] <- eps
  Gp <- Gamma_fast_fn(alpha_hat + ej)
  Gm <- Gamma_fast_fn(alpha_hat - ej)
  R_mat_numerical[, j] <- (as.vector(t(Gp)) - as.vector(t(Gm))) / (2 * eps)
}

cat("=== R_mat comparison ===\n")
cat(sprintf("Max diff (analytical vs numerical): %.4e\n", max(abs(R_mat_analytical - R_mat_numerical))))

# S_mat
S_mat <- matrix(0, h^2, k)
for (j in 1:k) {
  dg_j <- t_Gamma_arr[, j, ]
  cross <- (crossprod(dg_j, g_mat) + crossprod(g_mat, dg_j)) / n
  S_mat[, j] <- as.vector(cross)
}

# S_mat numerical
S_mat_numerical <- matrix(0, h^2, k)
for (j in 1:k) {
  ej <- rep(0, k); ej[j] <- eps
  gp <- g_fn(alpha_hat + ej)
  gm <- g_fn(alpha_hat - ej)
  Omega_p <- crossprod(gp) / n
  Omega_m <- crossprod(gm) / n
  S_mat_numerical[, j] <- as.vector((Omega_p - Omega_m) / (2 * eps))
}

cat(sprintf("Max diff S_mat (analytical vs numerical): %.4e\n", max(abs(S_mat - S_mat_numerical))))

# Now compute H2_mat
GtW <- crossprod(Gamma_hat, W_hat)
GtWG <- GtW %*% Gamma_hat
GW_row <- as.vector(t(eta_s) %*% W_hat)
GW_kron_Ik <- kronecker(matrix(GW_row, 1, h), diag(k))
GW_kron_GtW <- kronecker(matrix(GW_row, 1, h), GtW)
H2_analytical <- GtWG + GW_kron_Ik %*% R_mat_analytical - GW_kron_GtW %*% S_mat

cat("\n=== H2_mat comparison ===\n")
cat("H2 analytical:\n")
print(round(H2_analytical, 6))
cat("\nH2 numerical (finite diff of FOC):\n")
print(round(H2_numerical, 6))
cat(sprintf("\nMax diff H2 (analytical vs numerical): %.4e\n", max(abs(H2_analytical - H2_numerical))))
cat(sprintf("Relative max diff: %.4e\n", max(abs(H2_analytical - H2_numerical)) / max(abs(H2_numerical))))

# ==========================================
# 2. Verify Q_i numerically
# ==========================================
# Q_i = Gamma'W g_i + Gamma_i'WG - (Gamma'W g_i)(g_i'WG)

WG <- as.vector(W_hat %*% eta_s)
GtW_g <- GtW %*% t(g_mat)  # k x n
g_WG <- as.vector(g_mat %*% WG)  # n-vector

# Gamma_i'WG for each i
dim_arr <- dim(t_Gamma_arr)
t_Gamma_mat <- matrix(t_Gamma_arr, nrow = dim_arr[1] * dim_arr[2], ncol = dim_arr[3])
GammaT_WG <- t(matrix(t_Gamma_mat %*% WG, n, k))  # k x n

Q_mat <- GtW_g + GammaT_WG - sweep(GtW_g, 2, g_WG, `*`)

cat("\n=== Q_i verification ===\n")
cat(sprintf("mean(Q_i): %s\n", paste(sprintf("%.6e", rowMeans(Q_mat)), collapse=", ")))
cat(sprintf("Expected (Gamma'WG): %s\n", paste(sprintf("%.6e", as.vector(GtW %*% eta_s)), collapse=", ")))

# ==========================================
# 3. Compare SE2 from package with manual reconstruction
# ==========================================
H2_inv <- solve(H2_analytical)
psi2_manual <- -(H2_inv %*% Q_mat)  # k x n
se2_manual <- sqrt(diag(var(t(psi2_manual)) / n))

cat("\n=== SE2 comparison ===\n")
cat(sprintf("Package se2: %s\n", paste(sprintf("%.6f", gmm_fit$se2), collapse=", ")))
cat(sprintf("Manual se2:  %s\n", paste(sprintf("%.6f", se2_manual), collapse=", ")))

# Also with numerical H2
H2_inv_num <- solve(H2_numerical)
psi2_numH <- -(H2_inv_num %*% Q_mat)
se2_numH <- sqrt(diag(var(t(psi2_numH)) / n))
cat(sprintf("Num-H se2:  %s\n", paste(sprintf("%.6f", se2_numH), collapse=", ")))

# ==========================================
# 4. Eigenanalysis of H2_mat
# ==========================================
cat("\n=== H2_mat eigenvalues ===\n")
eig_H2 <- eigen(H2_analytical, symmetric = FALSE)
cat(sprintf("Eigenvalues: %s\n", paste(sprintf("%.6e", eig_H2$values), collapse=", ")))
cat(sprintf("Condition number: %.4e\n", max(Mod(eig_H2$values)) / min(Mod(eig_H2$values))))

eig_GtWG <- eigen(GtWG, symmetric = TRUE)
cat(sprintf("\nGtWG eigenvalues: %s\n", paste(sprintf("%.6e", eig_GtWG$values), collapse=", ")))
cat(sprintf("GtWG condition: %.4e\n", max(eig_GtWG$values) / min(eig_GtWG$values)))

# ==========================================
# 5. Contribution of each H2 term
# ==========================================
R_term <- GW_kron_Ik %*% R_mat_analytical
S_term <- GW_kron_GtW %*% S_mat

cat("\n=== H2 term norms ===\n")
cat(sprintf("||GtWG||_F: %.6e\n", norm(GtWG, "F")))
cat(sprintf("||(G'W⊗I)R||_F: %.6e\n", norm(R_term, "F")))
cat(sprintf("||(G'W⊗Γ'W)S||_F: %.6e\n", norm(S_term, "F")))
cat(sprintf("||R_term||/||GtWG||: %.4f\n", norm(R_term, "F") / norm(GtWG, "F")))
cat(sprintf("||S_term||/||GtWG||: %.4f\n", norm(S_term, "F") / norm(GtWG, "F")))
