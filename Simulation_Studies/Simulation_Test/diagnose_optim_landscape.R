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

# Pick an outlier rep
rep_i <- 929
dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]

# === Solution 1: zero init (current method) ===
ebmr <- EBMRAlgorithmFast4$new("y", subset_ps_spec, dat, W_func)
fit <- ebmr$ps_fit.list[[1]]
alpha_zero <- fit$coefficients
h_x <- fit$h_x
design_matrix <- fit$design_matrix
r_vec <- as.numeric(dat$r)
n <- nrow(dat)
esteq_dim <- ncol(h_x)
param_dim <- ncol(design_matrix)

# === Solution 2: GLM init (not available in practice, but gives "good" solution) ===
# GLM uses y for all obs including r=0 — cheating, but gives reference solution
glm_fit <- glm(r ~ y + u1 + z2, family = binomial(link = "logit"), data = dat)
alpha_glm <- coef(glm_fit)

# Now run GMM starting from GLM solution
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

# Two-step GMM from GLM init
obj_identity <- function(param) {
  G_bar <- G_func(param)
  as.numeric(crossprod(G_bar))
}
opt1 <- optim(alpha_glm, obj_identity, method = "L-BFGS-B",
              lower = rep(-Inf, param_dim), upper = rep(Inf, param_dim),
              control = list(maxit = 1000))
alpha_glm_gmm <- opt1$par

for (t in 1:200) {
  g_mat <- g_func(alpha_glm_gmm)
  W_hat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(esteq_dim))
  obj2 <- function(param) {
    G_bar <- G_func(param)
    as.numeric(crossprod(G_bar, W_hat %*% G_bar))
  }
  opt_t <- optim(alpha_glm_gmm, obj2, method = "L-BFGS-B",
                 lower = rep(-Inf, param_dim), upper = rep(Inf, param_dim),
                 control = list(maxit = 1000))
  if (max(abs(opt_t$par - alpha_glm_gmm)) < 1e-6) break
  alpha_glm_gmm <- opt_t$par
}

# === Compare the two solutions ===
cat("=== Rep 929: Zero-init vs GLM-init solutions ===\n\n")
cat(sprintf("Alpha (zero-init): [%s]\n", paste(round(alpha_zero, 6), collapse=", ")))
cat(sprintf("Alpha (GLM-init):  [%s]\n", paste(round(alpha_glm_gmm, 6), collapse=", ")))
cat(sprintf("Alpha (GLM raw):   [%s]\n\n", paste(round(alpha_glm, 6), collapse=", ")))

# Propensity scores
pi_zero <- as.vector(compute_pi(design_matrix %*% alpha_zero))
pi_glm <- as.vector(compute_pi(design_matrix %*% alpha_glm_gmm))
cat(sprintf("min(pi) zero-init: %.6f,  GLM-init: %.6f\n", min(pi_zero), min(pi_glm)))
cat(sprintf("max(r/pi) zero-init: %.1f,  GLM-init: %.1f\n", max(r_vec/pi_zero), max(r_vec/pi_glm)))

mu_zero <- mean(r_vec / pi_zero * dat$y)
mu_glm <- mean(r_vec / pi_glm * dat$y)
cat(sprintf("mu_ipw zero-init: %.4f,  GLM-init: %.4f\n\n", mu_zero, mu_glm))

# === Objective values at both solutions ===
# Evaluate GMM objective Q(alpha) = G'WG at each solution
# Use the W computed at each solution's own g
g_zero <- g_func(alpha_zero)
G_zero <- matrix(colMeans(g_zero), esteq_dim, 1)
W_zero <- solve(crossprod(g_zero) / n)
obj_zero_self <- as.numeric(crossprod(G_zero, W_zero %*% G_zero))

g_glm2 <- g_func(alpha_glm_gmm)
G_glm2 <- matrix(colMeans(g_glm2), esteq_dim, 1)
W_glm2 <- solve(crossprod(g_glm2) / n)
obj_glm_self <- as.numeric(crossprod(G_glm2, W_glm2 %*% G_glm2))

# Cross-evaluate: use W from zero-init to evaluate GLM solution
obj_glm_at_Wzero <- as.numeric(crossprod(G_glm2, W_zero %*% G_glm2))
obj_zero_at_Wglm <- as.numeric(crossprod(G_zero, W_glm2 %*% G_zero))

cat("=== Objective values ===\n")
cat(sprintf("Q(alpha_zero, W_zero) = %.6e  (self-consistent)\n", obj_zero_self))
cat(sprintf("Q(alpha_glm,  W_glm)  = %.6e  (self-consistent)\n", obj_glm_self))
cat(sprintf("Q(alpha_glm,  W_zero) = %.6e  (GLM sol evaluated at zero-init's W)\n", obj_glm_at_Wzero))
cat(sprintf("Q(alpha_zero, W_glm)  = %.6e  (zero sol evaluated at GLM-init's W)\n\n", obj_zero_at_Wglm))

# === Gradient at both solutions ===
# Numerical gradient of Q(alpha) = G(alpha)' W G(alpha) where W is FIXED
# Use W from zero-init solution
grad_Q <- function(param, W_fixed) {
  numDeriv::grad(function(p) {
    G_bar <- G_func(p)
    as.numeric(crossprod(G_bar, W_fixed %*% G_bar))
  }, param)
}

grad_zero_Wzero <- grad_Q(alpha_zero, W_zero)
grad_glm_Wzero <- grad_Q(alpha_glm_gmm, W_zero)
grad_zero_Wglm <- grad_Q(alpha_zero, W_glm2)
grad_glm_Wglm <- grad_Q(alpha_glm_gmm, W_glm2)

cat("=== Gradient norm at each solution ===\n")
cat(sprintf("||grad Q(alpha_zero, W_zero)|| = %.6e\n", max(abs(grad_zero_Wzero))))
cat(sprintf("||grad Q(alpha_glm,  W_zero)|| = %.6e\n", max(abs(grad_glm_Wzero))))
cat(sprintf("||grad Q(alpha_zero, W_glm) || = %.6e\n", max(abs(grad_zero_Wglm))))
cat(sprintf("||grad Q(alpha_glm,  W_glm) || = %.6e\n\n", max(abs(grad_glm_Wglm))))

# === Hessian at both solutions ===
hess_Q <- function(param, W_fixed) {
  numDeriv::hessian(function(p) {
    G_bar <- G_func(p)
    as.numeric(crossprod(G_bar, W_fixed %*% G_bar))
  }, param)
}

cat("=== Hessian eigenvalues at zero-init solution (W_zero) ===\n")
H_zero <- hess_Q(alpha_zero, W_zero)
eig_zero <- eigen(H_zero)$values
cat(sprintf("  Eigenvalues: [%s]\n", paste(round(eig_zero, 6), collapse=", ")))
cat(sprintf("  All positive (local min)? %s\n\n", all(eig_zero > 0)))

cat("=== Hessian eigenvalues at GLM-init solution (W_glm) ===\n")
H_glm <- hess_Q(alpha_glm_gmm, W_glm2)
eig_glm <- eigen(H_glm)$values
cat(sprintf("  Eigenvalues: [%s]\n", paste(round(eig_glm, 6), collapse=", ")))
cat(sprintf("  All positive (local min)? %s\n\n", all(eig_glm > 0)))

# === Cross-check: Hessian at zero-init solution using W_glm ===
cat("=== Hessian eigenvalues at zero-init solution (W_glm) ===\n")
H_zero_Wglm <- hess_Q(alpha_zero, W_glm2)
eig_zero_Wglm <- eigen(H_zero_Wglm)$values
cat(sprintf("  Eigenvalues: [%s]\n", paste(round(eig_zero_Wglm, 6), collapse=", ")))
cat(sprintf("  All positive? %s\n\n", all(eig_zero_Wglm > 0)))

# === Moment conditions comparison ===
cat("=== Moment conditions G_bar at each solution ===\n")
cat(sprintf("G(alpha_zero): [%s]\n", paste(round(as.vector(G_zero), 6), collapse=", ")))
cat(sprintf("G(alpha_glm):  [%s]\n\n", paste(round(as.vector(G_glm2), 6), collapse=", ")))

# === Variance of g (diagonal of g'g/n) ===
var_g_zero <- diag(crossprod(g_zero) / n)
var_g_glm <- diag(crossprod(g_glm2) / n)
cat("=== Variance of moment conditions (diag of g'g/n) ===\n")
cat(sprintf("Var(g) zero-init: [%s]\n", paste(round(var_g_zero, 4), collapse=", ")))
cat(sprintf("Var(g) GLM-init:  [%s]\n\n", paste(round(var_g_glm, 4), collapse=", ")))

# === Path analysis: what happens during optimization from zero? ===
cat("=== Step-by-step optimization from zero ===\n")
alpha_path <- rep(0, param_dim)

# Step 1: W=I
obj_I <- function(param) {
  G_bar <- G_func(param)
  as.numeric(crossprod(G_bar))
}
opt_step1 <- optim(alpha_path, obj_I, method = "L-BFGS-B",
                   lower = rep(-Inf, param_dim), upper = rep(Inf, param_dim),
                   control = list(maxit = 1000))
alpha_step1 <- opt_step1$par
cat(sprintf("After W=I step:  alpha=[%s], obj=%.6e\n",
            paste(round(alpha_step1, 4), collapse=", "), opt_step1$value))

# Step 2-5: iterative W updates
alpha_curr <- alpha_step1
for (step in 2:5) {
  g_curr <- g_func(alpha_curr)
  W_curr <- tryCatch(solve(crossprod(g_curr) / n), error = function(e) diag(esteq_dim))
  obj_curr <- function(param) {
    G_bar <- G_func(param)
    as.numeric(crossprod(G_bar, W_curr %*% G_bar))
  }
  opt_curr <- optim(alpha_curr, obj_curr, method = "L-BFGS-B",
                    lower = rep(-Inf, param_dim), upper = rep(Inf, param_dim),
                    control = list(maxit = 1000))
  alpha_curr <- opt_curr$par

  pi_curr <- as.vector(compute_pi(design_matrix %*% alpha_curr))
  mu_curr <- mean(r_vec / pi_curr * dat$y)
  cat(sprintf("After step %d:    alpha=[%s], obj=%.6e, mu=%.4f, min_pi=%.6f\n",
              step, paste(round(alpha_curr, 4), collapse=", "),
              opt_curr$value, mu_curr, min(pi_curr)))
}

cat("\nDone!\n")
