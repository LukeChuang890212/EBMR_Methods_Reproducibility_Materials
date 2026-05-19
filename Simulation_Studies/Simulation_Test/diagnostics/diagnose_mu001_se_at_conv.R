## Compute SE at the adaptively-damped converged solution directly
## Without letting the package re-iterate
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EPS", quiet = TRUE)
library(numDeriv)

ps_spec_base <- get_ps_spec("9-alt1")
h_alpha_fn <- function(dat) cbind(
  u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2,
  u1u2 = dat$u1*dat$u2, z1z2 = dat$z1*dat$z2,
  u1z1 = dat$u1*dat$z1, z1u2 = dat$z1*dat$u2,
  u1z2 = dat$u1*dat$z2, u2z2 = dat$u2*dat$z2)

data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]]
all_data <- NULL
for (fi in seq_along(data_file)) {
  if (file.exists(data_file[[fi]])) {
    test <- readRDS(data_file[[fi]])
    if (nrow(test) / 1000 == 2000) { all_data <- test; break }
  }
}

dat <- all_data[1:2000, ]
design_mat <- model.matrix(ps_spec_base$formula.list[[3]], dat)
r_vec <- as.vector(dat$r)
y_vec <- as.vector(dat$y)
n <- 2000
h_x <- cbind(1, h_alpha_fn(dat))
h_dim <- ncol(h_x)
alpha_dim <- ncol(design_mat)

compute_pi <- function(eta) 1/(1+exp(eta))

g_fn <- function(alpha) {
  eta <- as.vector(design_mat %*% alpha)
  eta <- pmin(pmax(eta, -20), 20)
  pi_vec <- compute_pi(eta)
  as.vector(r_vec/pi_vec - 1) * h_x
}
G_fn <- function(alpha) colMeans(g_fn(alpha))

# Run adaptive damping to get converged alpha
inner_opt <- function(alpha, W_hat) {
  obj <- function(a) {
    G <- matrix(G_fn(a), ncol=1)
    as.numeric(t(G) %*% W_hat %*% G)
  }
  optim(alpha, obj, method = "L-BFGS-B", control = list(maxit = 1000))$par
}

Gamma_fn <- function(alpha) {
  eta <- as.vector(design_mat %*% alpha)
  eta <- pmin(pmax(eta, -20), 20)
  pi_vec <- compute_pi(eta)
  common <- r_vec * (1 - pi_vec) / pi_vec
  Gamma_mat <- matrix(0, h_dim, alpha_dim)
  for (j in 1:alpha_dim) {
    Gamma_mat[, j] <- colMeans(common * design_mat[, j] * h_x)
  }
  Gamma_mat
}

conv_check <- function(alpha) {
  g_mat <- g_fn(alpha)
  G <- colMeans(g_mat)
  W <- tryCatch(solve(var(g_mat)), error = function(e) diag(h_dim))
  Gamma_mat <- Gamma_fn(alpha)
  grad <- 2 * as.vector(t(Gamma_mat) %*% W %*% G)
  max(abs(grad))
}

cat("=== Running adaptive damping to convergence ===\n")
alpha <- rep(0, alpha_dim)
lambda <- 0.0
W_prev <- diag(h_dim)
prev_gn <- Inf
for (t in 1:1000) {
  g_mat <- g_fn(alpha)
  W_new <- tryCatch(solve(var(g_mat)), error = function(e) diag(h_dim))
  W_hat <- lambda * W_prev + (1 - lambda) * W_new
  W_prev <- W_hat
  alpha <- inner_opt(alpha, W_hat)
  gn <- conv_check(alpha)
  if (gn >= prev_gn) {
    lambda <- min(lambda + 0.05, 0.95)
  } else {
    lambda <- max(lambda - 0.01, 0.0)
  }
  prev_gn <- gn
  if (gn < 1e-6) { cat(sprintf("  Converged at iter %d, grad=%.2e\n", t, gn)); break }
}
if (gn >= 1e-6) cat(sprintf("  Did not converge: iter=%d, grad=%.2e\n", t, gn))

alpha_conv <- alpha
cat(sprintf("  alpha = (%s)\n", paste(round(alpha_conv, 4), collapse=", ")))
cat(sprintf("  max|G| = %.4e\n\n", max(abs(G_fn(alpha_conv)))))

# Now compute SE at this converged point using the formula manually
cat("=== Computing SE at converged point ===\n\n")

g_mat <- g_fn(alpha_conv)
G_vec <- G_fn(alpha_conv)
W_hat <- solve(var(g_mat))
Gamma_hat <- Gamma_fn(alpha_conv)
pi_hat <- compute_pi(pmin(pmax(as.vector(design_mat %*% alpha_conv), -20), 20))

# dg/dalpha for each observation: n x k x h array
dg_arr <- array(0, dim = c(n, alpha_dim, h_dim))
common_factor <- r_vec * (1 - pi_hat) / pi_hat
for (l in 1:h_dim) {
  dg_arr[, , l] <- common_factor * h_x[, l] * design_mat
}

# R_mat: d vec(Gamma') / d theta  [k*h x k]
# Use finite differences on Gamma
R_mat <- matrix(0, alpha_dim * h_dim, alpha_dim)
eps <- 1e-5
for (j in 1:alpha_dim) {
  ej <- rep(0, alpha_dim); ej[j] <- eps
  G_plus <- Gamma_fn(alpha_conv + ej)
  G_minus <- Gamma_fn(alpha_conv - ej)
  R_mat[, j] <- (as.vector(t(G_plus)) - as.vector(t(G_minus))) / (2 * eps)
}

# S_mat: d vec(W^{-1}) / d theta  [h^2 x k]
# W^{-1} = var(g) = (1/(n-1)) (g-gbar)'(g-gbar)
g_centered <- scale(g_mat, center = TRUE, scale = FALSE)
S_mat <- matrix(0, h_dim^2, alpha_dim)
for (j in 1:alpha_dim) {
  dg_j <- dg_arr[, j, ]  # n x h_dim
  dg_j_centered <- scale(dg_j, center = TRUE, scale = FALSE)
  cross <- (crossprod(dg_j_centered, g_centered) + crossprod(g_centered, dg_j_centered)) / (n - 1)
  S_mat[, j] <- as.vector(cross)
}

# H = Gamma'W Gamma + (G'W ⊗ I_k) R - (G'W ⊗ Gamma'W) S
GtW <- t(Gamma_hat) %*% W_hat
GtWG <- GtW %*% Gamma_hat
GW_row <- as.vector(t(G_vec) %*% W_hat)
GW_kron_Ik <- kronecker(matrix(GW_row, 1, h_dim), diag(alpha_dim))
GW_kron_GtW <- kronecker(matrix(GW_row, 1, h_dim), GtW)
H_mat <- GtWG + GW_kron_Ik %*% R_mat - GW_kron_GtW %*% S_mat

cat("H matrix:\n")
cat(sprintf("  cond(H) = %.2e\n", kappa(H_mat)))
cat(sprintf("  cond(GtWG) = %.2e\n", kappa(GtWG)))
cat(sprintf("  norm(GW_kron_Ik %%*%% R) / norm(GtWG) = %.4f\n",
    norm(GW_kron_Ik %*% R_mat, "F") / norm(GtWG, "F")))
cat(sprintf("  norm(GW_kron_GtW %%*%% S) / norm(GtWG) = %.4f\n\n",
    norm(GW_kron_GtW %*% S_mat, "F") / norm(GtWG, "F")))

# Q_i = Gamma'W g_i + Gamma_i'WG - (Gamma'W g_i)(g_i'WG)
WG <- as.vector(W_hat %*% G_vec)
GtW_g <- GtW %*% t(g_mat)  # k x n
g_WG <- as.vector(g_mat %*% WG)  # n-vector

# Gamma_i'WG
t_Gamma_mat <- matrix(dg_arr, nrow = n * alpha_dim, ncol = h_dim)
GammaT_WG <- t(matrix(t_Gamma_mat %*% WG, n, alpha_dim))  # k x n

Q_mat <- GtW_g + GammaT_WG - sweep(GtW_g, 2, g_WG, `*`)  # k x n

# psi = -H^{-1} Q
H_inv <- solve(H_mat)
psi <- -(H_inv %*% Q_mat)  # k x n

alpha_se <- sqrt(diag(var(t(psi)) / n))
cat(sprintf("Alpha SE at converged point: %s\n", paste(round(alpha_se, 6), collapse=", ")))
cat(sprintf("max|psi|: %.4f\n\n", max(abs(psi))))

# Now compute IPW SE
mu_ipw <- mean(r_vec * y_vec / pi_hat)
ry_pi2 <- r_vec * y_vec / (pi_hat^2)
dot_pi <- -design_mat * (pi_hat * (1 - pi_hat))
H_alpha_w <- colMeans(ry_pi2 * dot_pi)

mu_iid <- as.vector(r_vec * y_vec / pi_hat) - as.vector(t(H_alpha_w) %*% psi)
se_ipw <- sqrt(var(mu_iid) / n)

cat(sprintf("mu_ipw = %.6f\n", mu_ipw))
cat(sprintf("SE at converged point = %.6f\n", se_ipw))
cat(sprintf("SE naive (no alpha uncert) = %.6f\n\n", sqrt(var(r_vec * y_vec / pi_hat) / n)))

# Compare with package (which doesn't converge)
cat("--- Package result (non-converged) ---\n")
ps_spec <- list(
  formula.list = ps_spec_base$formula.list,
  h_alpha.list = list(h_alpha_fn, h_alpha_fn, h_alpha_fn),
  inv_link = ps_spec_base$inv_link,
  outcome = ps_spec_base$outcome
)
W_fn <- function(g.matrix) solve(var(g.matrix))
ebmr <- EPS$new("y", ps_spec, dat, W_fn)
ipw_pkg <- ebmr$EBMR_IPW(h_alpha_fn, model_indices = 3, se.fit = TRUE)
cat(sprintf("mu_ipw = %.6f, SE = %.6f\n", ipw_pkg$mu_ipw, ipw_pkg$se_ipw))
