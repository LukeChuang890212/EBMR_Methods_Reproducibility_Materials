## Test M3 without stall limit — let it run many outer iterations
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

# Manually run iterative CUE GMM for M3 without stall limit
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

# Start from zeros
alpha <- rep(0, alpha_dim)

cat("=== Running 200 outer iterations without stall limit ===\n\n")
cat(sprintf("%4s | %12s | %12s | %12s\n", "iter", "obj", "max|G|", "grad_norm"))
cat(strrep("-", 55), "\n")

best_grad <- Inf
best_alpha <- alpha

for (t in 1:200) {
  g_mat <- g_fn(alpha)
  W_hat <- tryCatch(solve(var(g_mat)), error = function(e) diag(h_dim))
  G_vec <- colMeans(g_mat)
  obj <- as.numeric(t(G_vec) %*% W_hat %*% G_vec)

  opt <- optim(alpha, function(a) {
    G <- matrix(G_fn(a), ncol=1)
    as.numeric(t(G) %*% W_hat %*% G)
  }, method = "L-BFGS-B", control = list(maxit = 1000))
  alpha <- opt$par

  # Compute conv grad at new alpha with updated W
  g_new <- g_fn(alpha)
  W_new <- tryCatch(solve(var(g_new)), error = function(e) diag(h_dim))
  G_new <- colMeans(g_new)
  grad_norm <- max(abs(G_new))

  if (grad_norm < best_grad) {
    best_grad <- grad_norm
    best_alpha <- alpha
  }

  if (t <= 30 || t %% 10 == 0) {
    cat(sprintf("%4d | %12.6e | %12.6e | %12.6e\n", t, obj, max(abs(G_vec)), grad_norm))
  }
}

cat(sprintf("\nBest grad_norm achieved: %.6e\n", best_grad))
cat(sprintf("Best alpha: %s\n", paste(round(best_alpha, 4), collapse=", ")))

# Final check
g_best <- g_fn(best_alpha)
G_best <- colMeans(g_best)
cat(sprintf("G at best: %s\n", paste(round(G_best, 6), collapse=", ")))
