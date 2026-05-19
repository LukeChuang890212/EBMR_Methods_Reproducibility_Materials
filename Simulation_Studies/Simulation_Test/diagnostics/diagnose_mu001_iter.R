## Test if more outer iterations fix M3 convergence for rep 1
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

W_var <- function(g.matrix) solve(var(g.matrix))

data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]]
all_data <- NULL
for (fi in seq_along(data_file)) {
  if (file.exists(data_file[[fi]])) {
    test <- readRDS(data_file[[fi]])
    if (nrow(test) / 1000 == 2000) { all_data <- test; break }
  }
}

dat <- all_data[1:2000, ]

# Manually call WangShaoKim2014 for M3 with different max_outer settings
# We need to access the gmm function directly
# Let's test by modifying stall_limit

cat("=== Testing M3 convergence with different max_outer / stall_limit ===\n\n")

for (max_outer in c(200, 500, 1000, 2000)) {
  ps_spec_test <- list(
    formula.list = list(ps_spec_base$formula.list[[3]]),
    h_alpha.list = list(h_alpha_fn),
    inv_link = ps_spec_base$inv_link,
    outcome = ps_spec_base$outcome
  )

  # Use the package but override stall_limit via ps_specifications
  # Actually, stall_limit is hardcoded in gmm(). Let's just try more outer iters
  # by increasing stall_limit in the package call... but we can't easily.
  # Instead, call gmm directly.

  # Build ebmr to get the design matrix etc
  ebmr <- EPS$new("y", ps_spec_test, dat, W_var)
  fit <- ebmr$ps_fit.list[[1]]
  cat(sprintf("max_outer=default: grad=%.2e, obj=%.2e, iter=%d, converged=%s\n",
      fit$gmm_fit$opt$final_grad_norm, fit$gmm_fit$opt$objective,
      fit$gmm_fit$opt$iterations, fit$gmm_fit$opt$converged))
}

# The real question: is this a stall issue or a fundamental non-convergence?
# Let's look at how grad_norm evolves
cat("\n=== Tracking grad_norm across outer iterations (manual GMM) ===\n\n")

# Reconstruct M3's g function manually
formula3 <- ps_spec_base$formula.list[[3]]
inv_link <- ps_spec_base$inv_link

# Get the preprocessed pieces from the fitted object
ebmr <- EPS$new("y", list(
  formula.list = list(formula3),
  h_alpha.list = list(h_alpha_fn),
  inv_link = inv_link,
  outcome = ps_spec_base$outcome
), dat, W_var)

fit3 <- ebmr$ps_fit.list[[1]]
alpha_hat <- fit3$coefficients
g_mat <- fit3$gmm_fit$g.matrix
G_hat <- colMeans(g_mat)

cat(sprintf("Final alpha: %s\n", paste(round(alpha_hat, 4), collapse=", ")))
cat(sprintf("Final |G|: max=%.4e, norm=%.4e\n", max(abs(G_hat)), sqrt(sum(G_hat^2))))
cat(sprintf("Converged: %s, iterations: %d\n", fit3$gmm_fit$opt$converged, fit3$gmm_fit$opt$iterations))

# Check: is the issue that grad_norm=6.9e-3 is a true stationary point
# (i.e., can't go lower), or is it just stalling?
# If it's a saddle/plateau, more iterations won't help.
cat(sprintf("\nGrad norm at solution: %.6e\n", fit3$gmm_fit$opt$final_grad_norm))
cat("(tolerance is 1e-6, so 6.9e-3 means NOT converged)\n\n")

# Try with stall_limit = 200 (more patience)
# We can't easily change stall_limit without modifying the package,
# but we CAN check if running optim again from this point helps
cat("=== Running additional optim from current solution ===\n\n")

# Reconstruct objective
design_mat <- fit3$design_matrix
r_vec <- as.vector(dat$r)
n <- 2000
h_x2 <- h_alpha_fn(dat)
h_x <- cbind(1, h_x2)
h_dim <- ncol(h_x)
alpha_dim <- ncol(design_mat)

compute_pi <- function(eta) 1/(1+exp(eta))  # logistic_complement

g_fn <- function(alpha) {
  eta <- as.vector(design_mat %*% alpha)
  eta <- pmin(pmax(eta, -20), 20)
  pi_vec <- compute_pi(eta)
  as.vector(r_vec/pi_vec - 1) * h_x
}

G_fn <- function(alpha) colMeans(g_fn(alpha))

# Current W
g_cur <- g_fn(alpha_hat)
W_cur <- solve(var(g_cur))

obj_fn <- function(alpha) {
  G <- matrix(G_fn(alpha), ncol=1)
  as.numeric(t(G) %*% W_cur %*% G)
}

cat(sprintf("Obj at current: %.6e\n", obj_fn(alpha_hat)))

# Try more iterations
for (attempt in 1:5) {
  g_cur <- g_fn(alpha_hat)
  W_cur <- solve(var(g_cur))
  opt <- optim(alpha_hat, function(a) {
    G <- matrix(G_fn(a), ncol=1)
    as.numeric(t(G) %*% W_cur %*% G)
  }, method = "L-BFGS-B", control = list(maxit = 5000))
  alpha_hat <- opt$par
  G_new <- G_fn(alpha_hat)
  cat(sprintf("  Attempt %d: obj=%.6e, max|G|=%.4e, convergence=%d\n",
      attempt, opt$value, max(abs(G_new)), opt$convergence))
}
