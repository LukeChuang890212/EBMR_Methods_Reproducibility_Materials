setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

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

ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
gmm_fit <- ebmr$ps_fit.list[[1]]$gmm_fit
design_mat <- ebmr$ps_fit.list[[1]]$design_matrix
link_type <- ebmr$ps_fit.list[[1]]$link_type
alpha_hat <- gmm_fit$estimates
n <- n_val; k <- length(alpha_hat)
r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
pi_hat <- ebmr$ps_fit.list[[1]]$fitted.values

h_alpha_vars <- ps_spec[["h_alpha.list"]][[model_idx]]
h_x <- cbind(1, as.matrix(dat[, h_alpha_vars, drop = FALSE]))
h_dim <- ncol(h_x)

compute_pi_fn <- function(eta) plogis(-eta)

g_fn <- function(alpha) {
  eta <- as.vector(design_mat %*% alpha)
  pi_v <- compute_pi_fn(eta)
  (r_vec / pi_v - 1) * h_x
}

score_fn <- function(alpha) {
  g_mat <- g_fn(alpha)
  G_bar <- colMeans(g_mat)
  W_mat <- solve(crossprod(g_mat) / n)
  eta <- as.vector(design_mat %*% alpha)
  pi_v <- compute_pi_fn(eta)
  cf <- r_vec * (1 - pi_v) / pi_v
  Gamma_mat <- crossprod(h_x * cf, design_mat) / n
  as.vector(crossprod(Gamma_mat, W_mat %*% G_bar))
}

# ============================================================
# Convergence of numerical H2 with different eps
# ============================================================
cat("=== H2 numerical convergence ===\n")
for (eps in c(1e-3, 1e-4, 1e-5, 1e-6, 1e-7)) {
  H2_num <- matrix(0, k, k)
  for (j in 1:k) {
    ej <- rep(0, k); ej[j] <- eps
    H2_num[, j] <- (score_fn(alpha_hat + ej) - score_fn(alpha_hat - ej)) / (2 * eps)
  }
  cat(sprintf("eps=%.0e: H2_num[4,4]=%.6f, norm=%.6f\n", eps, H2_num[4,4], norm(H2_num, "F")))
}

# ============================================================
# Analytical H2 components
# ============================================================
Gamma_hat <- gmm_fit$Gamma.hat
g_mat <- gmm_fit$g.matrix
W_hat <- gmm_fit$W.hat
eta_s <- gmm_fit$eta_s
GtW <- crossprod(Gamma_hat, W_hat)
GtWG <- GtW %*% Gamma_hat

dc_deta <- r_vec * (1 - pi_hat) / pi_hat
dc_X <- dc_deta * design_mat
R_mat <- matrix(0, h_dim * k, k)
for (j in 1:k) {
  Xj_h <- design_mat[, j] * h_x
  block <- crossprod(dc_X, Xj_h) / n
  for (l in 1:h_dim) R_mat[(l-1)*k + j, ] <- block[, l]
}

cf_all <- r_vec * (1 - pi_hat) / pi_hat
S_mat <- matrix(0, h_dim^2, k)
for (j in 1:k) {
  dg_j <- cf_all * design_mat[, j] * h_x
  cross <- (crossprod(dg_j, g_mat) + crossprod(g_mat, dg_j)) / n
  S_mat[, j] <- as.vector(cross)
}

GW_row <- as.vector(t(eta_s) %*% W_hat)
GW_kron_Ik <- kronecker(matrix(GW_row, 1, h_dim), diag(k))
GW_kron_GtW <- kronecker(matrix(GW_row, 1, h_dim), GtW)

term1 <- GtWG
term2 <- GW_kron_Ik %*% R_mat
term3 <- GW_kron_GtW %*% S_mat
H2_analytical <- term1 + term2 - term3

cat(sprintf("\nH2_analytical[4,4]=%.6f\n", H2_analytical[4,4]))
cat("\nContributions to H2:\n")
cat(sprintf("  GtWG (term1):\n")); print(round(term1, 4))
cat(sprintf("  (G'W ⊗ I) R (term2):\n")); print(round(term2, 4))
cat(sprintf("  (G'W ⊗ Γ'W) S (term3):\n")); print(round(term3, 4))

# ============================================================
# Check: does the term3 (S contribution) have a sign/scale issue?
# Verify S_mat numerically
# ============================================================
cat("\n=== S_mat numerical verification ===\n")
eps_s <- 1e-5
S_num <- matrix(0, h_dim^2, k)
for (j in 1:k) {
  ej <- rep(0, k); ej[j] <- eps_s
  g_plus <- g_fn(alpha_hat + ej)
  g_minus <- g_fn(alpha_hat - ej)
  Winv_plus <- crossprod(g_plus) / n
  Winv_minus <- crossprod(g_minus) / n
  S_num[, j] <- as.vector((Winv_plus - Winv_minus) / (2 * eps_s))
}
cat(sprintf("||S_num - S_analytical|| = %.6e\n", norm(S_num - S_mat, "F")))
cat(sprintf("||S_num|| = %.6f\n", norm(S_num, "F")))
cat(sprintf("Relative diff = %.6e\n", norm(S_num - S_mat, "F") / norm(S_num, "F")))

# ============================================================
# Check R_mat numerically
# ============================================================
cat("\n=== R_mat numerical verification ===\n")
eps_r <- 1e-5
Gamma_at <- function(alpha) {
  eta <- as.vector(design_mat %*% alpha)
  pi_v <- compute_pi_fn(eta)
  cf <- r_vec * (1 - pi_v) / pi_v
  crossprod(h_x * cf, design_mat) / n  # h x k = Gamma
}
R_num <- matrix(0, h_dim * k, k)
for (j in 1:k) {
  ej <- rep(0, k); ej[j] <- eps_r
  Gp <- Gamma_at(alpha_hat + ej)
  Gm <- Gamma_at(alpha_hat - ej)
  # d Gamma / d alpha_j -> d vec(Gamma') / d alpha_j
  R_num[, j] <- as.vector(t((Gp - Gm) / (2 * eps_r)))
}
cat(sprintf("||R_num - R_analytical|| = %.6e\n", norm(R_num - R_mat, "F")))
cat(sprintf("||R_num|| = %.6f\n", norm(R_num, "F")))
cat(sprintf("Relative diff = %.6e\n", norm(R_num - R_mat, "F") / norm(R_num, "F")))

# Now recompute H2 with numerical R and S
H2_with_num_RS <- term1 + GW_kron_Ik %*% R_num - GW_kron_GtW %*% S_num

eps <- 1e-5
H2_num <- matrix(0, k, k)
for (j in 1:k) {
  ej <- rep(0, k); ej[j] <- eps
  H2_num[, j] <- (score_fn(alpha_hat + ej) - score_fn(alpha_hat - ej)) / (2 * eps)
}

cat(sprintf("\n||H2_num - H2_with_num_RS|| = %.6e  (rel=%.6e)\n",
    norm(H2_num - H2_with_num_RS, "F"),
    norm(H2_num - H2_with_num_RS, "F") / norm(H2_num, "F")))
cat(sprintf("||H2_num - H2_analytical||  = %.6e  (rel=%.6e)\n",
    norm(H2_num - H2_analytical, "F"),
    norm(H2_num - H2_analytical, "F") / norm(H2_num, "F")))
