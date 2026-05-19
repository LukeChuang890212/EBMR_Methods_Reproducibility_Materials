## Try mechanisms to escape the limit cycle for M3 CUE with W=solve(var(g))
## The issue: iterative CUE oscillates between two points
## Ideas: damping, averaging, relaxation, partial W update
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EPS", quiet = TRUE)

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

Gamma_fn <- function(alpha) {
  eta <- as.vector(design_mat %*% alpha)
  eta <- pmin(pmax(eta, -20), 20)
  pi_vec <- compute_pi(eta)
  common <- r_vec * (1 - pi_vec) / pi_vec  # = -r/pi^2 * dpi/deta for logistic_complement
  Gamma_mat <- matrix(0, h_dim, alpha_dim)
  for (j in 1:alpha_dim) {
    Gamma_mat[, j] <- colMeans(common * design_mat[, j] * h_x)
  }
  Gamma_mat
}

inner_opt <- function(alpha, W_hat) {
  obj <- function(a) {
    G <- matrix(G_fn(a), ncol=1)
    as.numeric(t(G) %*% W_hat %*% G)
  }
  optim(alpha, obj, method = "L-BFGS-B", control = list(maxit = 1000))$par
}

conv_check <- function(alpha) {
  g_mat <- g_fn(alpha)
  G <- colMeans(g_mat)
  W <- tryCatch(solve(var(g_mat)), error = function(e) diag(h_dim))
  Gamma_mat <- Gamma_fn(alpha)
  grad <- 2 * as.vector(t(Gamma_mat) %*% W %*% G)
  max(abs(grad))
}

cat("=== Mechanisms to escape limit cycle ===\n\n")

# Strategy 1: Damped W update (convex combination of old and new W)
cat("--- 1. Damped W (lambda=0.5) ---\n")
alpha <- rep(0, alpha_dim)
W_prev <- diag(h_dim)
for (t in 1:200) {
  g_mat <- g_fn(alpha)
  W_new <- tryCatch(solve(var(g_mat)), error = function(e) diag(h_dim))
  W_hat <- 0.5 * W_prev + 0.5 * W_new
  W_prev <- W_hat
  alpha <- inner_opt(alpha, W_hat)
  gn <- conv_check(alpha)
  if (gn < 1e-6) { cat(sprintf("  Converged at iter %d!\n", t)); break }
}
cat(sprintf("  iter=%d, grad_norm=%.4e, max|G|=%.4e\n\n", t, gn, max(abs(G_fn(alpha)))))

# Strategy 2: Damped W (lambda=0.9 — more aggressive damping)
cat("--- 2. Damped W (lambda=0.9) ---\n")
alpha <- rep(0, alpha_dim)
W_prev <- diag(h_dim)
for (t in 1:200) {
  g_mat <- g_fn(alpha)
  W_new <- tryCatch(solve(var(g_mat)), error = function(e) diag(h_dim))
  W_hat <- 0.9 * W_prev + 0.1 * W_new
  W_prev <- W_hat
  alpha <- inner_opt(alpha, W_hat)
  gn <- conv_check(alpha)
  if (gn < 1e-6) { cat(sprintf("  Converged at iter %d!\n", t)); break }
}
cat(sprintf("  iter=%d, grad_norm=%.4e, max|G|=%.4e\n\n", t, gn, max(abs(G_fn(alpha)))))

# Strategy 3: Damped alpha update (average of old and new alpha)
cat("--- 3. Damped alpha (average old and new) ---\n")
alpha <- rep(0, alpha_dim)
for (t in 1:200) {
  g_mat <- g_fn(alpha)
  W_hat <- tryCatch(solve(var(g_mat)), error = function(e) diag(h_dim))
  alpha_new <- inner_opt(alpha, W_hat)
  alpha <- 0.5 * alpha + 0.5 * alpha_new  # average
  gn <- conv_check(alpha)
  if (gn < 1e-6) { cat(sprintf("  Converged at iter %d!\n", t)); break }
}
cat(sprintf("  iter=%d, grad_norm=%.4e, max|G|=%.4e\n\n", t, gn, max(abs(G_fn(alpha)))))

# Strategy 4: Damped alpha (0.7 old + 0.3 new)
cat("--- 4. Damped alpha (0.7 old + 0.3 new) ---\n")
alpha <- rep(0, alpha_dim)
for (t in 1:200) {
  g_mat <- g_fn(alpha)
  W_hat <- tryCatch(solve(var(g_mat)), error = function(e) diag(h_dim))
  alpha_new <- inner_opt(alpha, W_hat)
  alpha <- 0.7 * alpha + 0.3 * alpha_new
  gn <- conv_check(alpha)
  if (gn < 1e-6) { cat(sprintf("  Converged at iter %d!\n", t)); break }
}
cat(sprintf("  iter=%d, grad_norm=%.4e, max|G|=%.4e\n\n", t, gn, max(abs(G_fn(alpha)))))

# Strategy 5: Line search on the CUE gradient along the update direction
cat("--- 5. Line search on update direction ---\n")
alpha <- rep(0, alpha_dim)
for (t in 1:200) {
  g_mat <- g_fn(alpha)
  W_hat <- tryCatch(solve(var(g_mat)), error = function(e) diag(h_dim))
  alpha_new <- inner_opt(alpha, W_hat)
  direction <- alpha_new - alpha
  # Line search: find step that minimizes conv_check
  best_step <- 1.0
  best_gn <- conv_check(alpha_new)
  for (s in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)) {
    cand <- alpha + s * direction
    gn_s <- conv_check(cand)
    if (gn_s < best_gn) { best_gn <- gn_s; best_step <- s }
  }
  alpha <- alpha + best_step * direction
  gn <- best_gn
  if (gn < 1e-6) { cat(sprintf("  Converged at iter %d!\n", t)); break }
}
cat(sprintf("  iter=%d, grad_norm=%.4e, max|G|=%.4e\n\n", t, gn, max(abs(G_fn(alpha)))))

# Strategy 6: Use running average of W from last K iterations
cat("--- 6. Running average W (last 5 iterations) ---\n")
alpha <- rep(0, alpha_dim)
W_history <- list()
for (t in 1:200) {
  g_mat <- g_fn(alpha)
  W_new <- tryCatch(solve(var(g_mat)), error = function(e) diag(h_dim))
  W_history <- c(W_history, list(W_new))
  if (length(W_history) > 5) W_history <- W_history[(length(W_history)-4):length(W_history)]
  W_hat <- Reduce("+", W_history) / length(W_history)
  alpha <- inner_opt(alpha, W_hat)
  gn <- conv_check(alpha)
  if (gn < 1e-6) { cat(sprintf("  Converged at iter %d!\n", t)); break }
}
cat(sprintf("  iter=%d, grad_norm=%.4e, max|G|=%.4e\n\n", t, gn, max(abs(G_fn(alpha)))))

# Strategy 7: Adaptive damping — increase damping when oscillation detected
cat("--- 7. Adaptive damping ---\n")
alpha <- rep(0, alpha_dim)
lambda <- 0.0  # start with no damping
W_prev <- diag(h_dim)
prev_gn <- Inf
for (t in 1:200) {
  g_mat <- g_fn(alpha)
  W_new <- tryCatch(solve(var(g_mat)), error = function(e) diag(h_dim))
  W_hat <- lambda * W_prev + (1 - lambda) * W_new
  W_prev <- W_hat
  alpha <- inner_opt(alpha, W_hat)
  gn <- conv_check(alpha)
  # If not improving, increase damping
  if (gn >= prev_gn) {
    lambda <- min(lambda + 0.05, 0.95)
  } else {
    lambda <- max(lambda - 0.01, 0.0)
  }
  prev_gn <- gn
  if (gn < 1e-6) { cat(sprintf("  Converged at iter %d!\n", t)); break }
}
cat(sprintf("  iter=%d, grad_norm=%.4e, max|G|=%.4e, final_lambda=%.2f\n\n",
    t, gn, max(abs(G_fn(alpha))), lambda))
