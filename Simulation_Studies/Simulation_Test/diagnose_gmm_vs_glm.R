setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
library(EBMRalgorithmFast4)

data_file <- misspecified_model_all_data_file.list$setting3$miss50[[1]]
all_data <- readRDS(data_file)
n_val <- 2000

ps_spec <- get_ps_spec("9-alt1")
subset_ps_spec <- list(
  formula.list = ps_spec$formula.list[2],
  h_alpha.list = ps_spec$h_alpha.list[2],
  inv_link = ps_spec$inv_link,
  outcome = ps_spec$outcome
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

eta_max <- 20
compute_pi <- function(eta) {
  eta <- pmin(pmax(eta, -eta_max), eta_max)
  1 / (1 + exp(eta))
}

rep_i <- 929
dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]

# GMM solution from package
ebmr <- EBMRAlgorithmFast4$new("y", subset_ps_spec, dat, W_func)
fit <- ebmr$ps_fit.list[[1]]
alpha_gmm <- fit$coefficients
h_x <- fit$h_x
design_matrix <- fit$design_matrix
r_vec <- as.numeric(dat$r)
n <- nrow(dat)
esteq_dim <- ncol(h_x)
param_dim <- ncol(design_matrix)

# GLM solution (uses y for r=0 — not available in practice)
glm_fit <- glm(r ~ y + u1 + z2, family = binomial(link = "logit"), data = dat)
alpha_glm <- coef(glm_fit)

cat(sprintf("Alpha (GMM): [%s]\n", paste(round(alpha_gmm, 6), collapse=", ")))
cat(sprintf("Alpha (GLM): [%s]\n\n", paste(round(alpha_glm, 6), collapse=", ")))

# Helper functions
g_func <- function(param) {
  eta <- design_matrix %*% param
  pi_vec <- compute_pi(eta)
  rw <- r_vec / as.vector(pi_vec)
  (rw - 1) * h_x
}

G_func <- function(param) {
  g_mat <- g_func(param)
  matrix(colMeans(g_mat), esteq_dim, 1)
}

# Compute Q(alpha, W) = G'WG with W = (g'g/n)^{-1} evaluated at alpha itself
Q_self <- function(param) {
  g_mat <- g_func(param)
  G_bar <- matrix(colMeans(g_mat), esteq_dim, 1)
  W <- solve(crossprod(g_mat) / n)
  as.numeric(crossprod(G_bar, W %*% G_bar))
}

# Compute Q with a FIXED W
Q_fixed <- function(param, W_fixed) {
  G_bar <- G_func(param)
  as.numeric(crossprod(G_bar, W_fixed %*% G_bar))
}

# Gradient of Q with fixed W
grad_Q <- function(param, W_fixed) {
  numDeriv::grad(function(p) Q_fixed(p, W_fixed), param)
}

# Hessian of Q with fixed W
hess_Q <- function(param, W_fixed) {
  numDeriv::hessian(function(p) Q_fixed(p, W_fixed), param)
}

# === Compute W at each solution ===
g_gmm <- g_func(alpha_gmm)
W_gmm <- solve(crossprod(g_gmm) / n)

g_glm <- g_func(alpha_glm)
W_glm <- solve(crossprod(g_glm) / n)

# === 1. Objective values ===
cat("=== OBJECTIVE Q(alpha, W) ===\n")
cat(sprintf("Q(alpha_gmm, W_gmm) = %.6e\n", Q_fixed(alpha_gmm, W_gmm)))
cat(sprintf("Q(alpha_glm, W_glm) = %.6e\n", Q_fixed(alpha_glm, W_glm)))
cat(sprintf("Q(alpha_glm, W_gmm) = %.6e  (GLM alpha, GMM's W)\n", Q_fixed(alpha_glm, W_gmm)))
cat(sprintf("Q(alpha_gmm, W_glm) = %.6e  (GMM alpha, GLM's W)\n\n", Q_fixed(alpha_gmm, W_glm)))

# === 2. Gradient ===
cat("=== GRADIENT of Q ===\n")
grad_gmm_Wgmm <- grad_Q(alpha_gmm, W_gmm)
grad_glm_Wglm <- grad_Q(alpha_glm, W_glm)
grad_glm_Wgmm <- grad_Q(alpha_glm, W_gmm)
grad_gmm_Wglm <- grad_Q(alpha_gmm, W_glm)

cat(sprintf("grad Q(alpha_gmm, W_gmm) = [%s]  ||.||=%.2e\n",
            paste(round(grad_gmm_Wgmm, 6), collapse=", "), max(abs(grad_gmm_Wgmm))))
cat(sprintf("grad Q(alpha_glm, W_glm) = [%s]  ||.||=%.2e\n",
            paste(round(grad_glm_Wglm, 6), collapse=", "), max(abs(grad_glm_Wglm))))
cat(sprintf("grad Q(alpha_glm, W_gmm) = [%s]  ||.||=%.2e\n",
            paste(round(grad_glm_Wgmm, 6), collapse=", "), max(abs(grad_glm_Wgmm))))
cat(sprintf("grad Q(alpha_gmm, W_glm) = [%s]  ||.||=%.2e\n\n",
            paste(round(grad_gmm_Wglm, 6), collapse=", "), max(abs(grad_gmm_Wglm))))

# === 3. Hessian eigenvalues ===
cat("=== HESSIAN eigenvalues ===\n")
H_gmm_Wgmm <- hess_Q(alpha_gmm, W_gmm)
H_glm_Wglm <- hess_Q(alpha_glm, W_glm)
H_glm_Wgmm <- hess_Q(alpha_glm, W_gmm)
H_gmm_Wglm <- hess_Q(alpha_gmm, W_glm)

eig_gmm <- eigen(H_gmm_Wgmm)$values
eig_glm <- eigen(H_glm_Wglm)$values
eig_glm_Wgmm <- eigen(H_glm_Wgmm)$values
eig_gmm_Wglm <- eigen(H_gmm_Wglm)$values

cat(sprintf("H(alpha_gmm, W_gmm) eig: [%s]  all>0: %s\n",
            paste(round(eig_gmm, 6), collapse=", "), all(eig_gmm > 0)))
cat(sprintf("H(alpha_glm, W_glm) eig: [%s]  all>0: %s\n",
            paste(round(eig_glm, 6), collapse=", "), all(eig_glm > 0)))
cat(sprintf("H(alpha_glm, W_gmm) eig: [%s]  all>0: %s\n",
            paste(round(eig_glm_Wgmm, 6), collapse=", "), all(eig_glm_Wgmm > 0)))
cat(sprintf("H(alpha_gmm, W_glm) eig: [%s]  all>0: %s\n\n",
            paste(round(eig_gmm_Wglm, 6), collapse=", "), all(eig_gmm_Wglm > 0)))

# === 4. Moment conditions ===
cat("=== MOMENT CONDITIONS G_bar ===\n")
G_gmm <- as.vector(G_func(alpha_gmm))
G_glm <- as.vector(G_func(alpha_glm))
cat(sprintf("G(alpha_gmm): [%s]\n", paste(round(G_gmm, 6), collapse=", ")))
cat(sprintf("G(alpha_glm): [%s]\n\n", paste(round(G_glm, 6), collapse=", ")))

# === 5. W matrices ===
cat("=== W MATRICES (diagonal) ===\n")
cat(sprintf("diag(W_gmm): [%s]\n", paste(round(diag(W_gmm), 4), collapse=", ")))
cat(sprintf("diag(W_glm): [%s]\n\n", paste(round(diag(W_glm), 4), collapse=", ")))

# === 6. Variance of g ===
cat("=== VARIANCE of g (diag of g'g/n) ===\n")
var_gmm <- diag(crossprod(g_gmm) / n)
var_glm <- diag(crossprod(g_glm) / n)
cat(sprintf("Var(g) at GMM: [%s]\n", paste(round(var_gmm, 2), collapse=", ")))
cat(sprintf("Var(g) at GLM: [%s]\n", paste(round(var_glm, 2), collapse=", ")))
cat(sprintf("Ratio GMM/GLM: [%s]\n\n", paste(round(var_gmm / var_glm, 2), collapse=", ")))

# === 7. IPW estimates ===
pi_gmm <- as.vector(compute_pi(design_matrix %*% alpha_gmm))
pi_glm <- as.vector(compute_pi(design_matrix %*% alpha_glm))
cat("=== IPW estimates ===\n")
cat(sprintf("mu_ipw at GMM: %.4f\n", mean(r_vec / pi_gmm * dat$y)))
cat(sprintf("mu_ipw at GLM: %.4f\n", mean(r_vec / pi_glm * dat$y)))
cat(sprintf("min(pi) GMM: %.6f,  GLM: %.6f\n", min(pi_gmm), min(pi_glm)))
cat(sprintf("max(r/pi) GMM: %.1f,  GLM: %.1f\n", max(r_vec/pi_gmm), max(r_vec/pi_glm)))

cat("\nDone!\n")
