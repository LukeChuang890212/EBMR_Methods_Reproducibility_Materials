## Check if inner-loop L-BFGS-B converges at each outer iteration
## during the adaptive damping procedure
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

cat("=== Adaptive damping: checking inner-loop convergence ===\n\n")
cat(sprintf("%4s | %5s | %12s | %12s | %12s | %5s\n",
    "iter", "lambda", "inner_obj", "outer_grad", "max|G|", "inner_conv"))
cat(strrep("-", 75), "\n")

alpha <- rep(0, alpha_dim)
lambda <- 0.0
W_prev <- diag(h_dim)
prev_gn <- Inf

for (t in 1:50) {
  g_mat <- g_fn(alpha)
  W_new <- tryCatch(solve(var(g_mat)), error = function(e) diag(h_dim))
  W_hat <- lambda * W_prev + (1 - lambda) * W_new
  W_prev <- W_hat

  # Inner optimization with full details
  obj_fn <- function(a) {
    G <- matrix(G_fn(a), ncol=1)
    as.numeric(t(G) %*% W_hat %*% G)
  }
  opt <- optim(alpha, obj_fn, method = "L-BFGS-B", control = list(maxit = 1000))
  alpha <- opt$par
  inner_conv <- opt$convergence  # 0 = converged

  gn <- conv_check(alpha)

  if (gn >= prev_gn) {
    lambda <- min(lambda + 0.05, 0.95)
  } else {
    lambda <- max(lambda - 0.01, 0.0)
  }
  prev_gn <- gn

  cat(sprintf("%4d | %5.2f | %12.6e | %12.4e | %12.4e | %5d\n",
      t, lambda, opt$value, gn, max(abs(G_fn(alpha))), inner_conv))

  if (gn < 1e-6) { cat("  *** CONVERGED ***\n"); break }
}
