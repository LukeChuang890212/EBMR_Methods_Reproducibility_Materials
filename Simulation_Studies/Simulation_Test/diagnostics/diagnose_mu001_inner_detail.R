## Verify inner-loop convergence more rigorously
## Check: does the inner optimizer actually reach the minimum of G'WG for fixed W?
## Verify by checking the gradient of the inner objective at the solution
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

cat("=== Detailed inner-loop convergence check (no damping) ===\n\n")
cat(sprintf("%4s | %12s | %12s | %12s | %12s\n",
    "iter", "inner_obj", "inner_grad", "inner_conv", "alpha_change"))
cat(strrep("-", 70), "\n")

alpha <- rep(0, alpha_dim)

for (t in 1:10) {
  g_mat <- g_fn(alpha)
  W_hat <- tryCatch(solve(var(g_mat)), error = function(e) diag(h_dim))

  # Inner objective
  obj_fn <- function(a) {
    G <- matrix(G_fn(a), ncol=1)
    as.numeric(t(G) %*% W_hat %*% G)
  }

  # Inner gradient (analytical): 2 * Gamma(alpha)' W G(alpha), W fixed
  grad_fn <- function(a) {
    G <- G_fn(a)
    Gamma_mat <- Gamma_fn(a)
    2 * as.vector(t(Gamma_mat) %*% W_hat %*% G)
  }

  # Run inner optimization
  opt <- optim(alpha, obj_fn, grad_fn, method = "L-BFGS-B",
               control = list(maxit = 1000))
  alpha_new <- opt$par

  # Verify: compute gradient at solution
  inner_grad_at_sol <- grad_fn(alpha_new)
  inner_grad_norm <- max(abs(inner_grad_at_sol))

  # Also verify with numerical gradient
  num_grad <- numDeriv::grad(obj_fn, alpha_new)
  num_grad_norm <- max(abs(num_grad))

  # Check objective can't be lowered further
  # Try a few random perturbations
  obj_at_sol <- obj_fn(alpha_new)
  can_improve <- FALSE
  for (trial in 1:20) {
    perturb <- alpha_new + rnorm(alpha_dim) * 0.01
    if (obj_fn(perturb) < obj_at_sol - 1e-10) {
      can_improve <- TRUE
      break
    }
  }

  alpha_change <- max(abs(alpha_new - alpha))

  cat(sprintf("%4d | %12.6e | %12.4e | %12.4e | %12.4e | can_improve=%s\n",
      t, obj_at_sol, inner_grad_norm, num_grad_norm, alpha_change,
      ifelse(can_improve, "YES!", "no")))
  cat(sprintf("       alpha = (%s)\n", paste(round(alpha_new, 4), collapse=", ")))
  cat(sprintf("       grad  = (%s)\n\n", paste(round(inner_grad_at_sol, 6), collapse=", ")))

  alpha <- alpha_new
}
