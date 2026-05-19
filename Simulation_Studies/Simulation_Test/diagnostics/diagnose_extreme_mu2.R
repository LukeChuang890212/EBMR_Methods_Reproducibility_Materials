setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
source("Simulation.r")
library(EPS)

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

# For specific extreme reps, trace what happens in each step
extreme_reps <- c(507, 255, 997, 547, 711, 989, 437)

cat(sprintf("%-6s %-12s %-12s %-12s %-12s %-12s %-12s\n",
            "Rep", "mu_final", "||alpha||", "min_pi", "max_r/pi",
            "Q_s1", "Q_s2"))

for (rep_i in extreme_reps) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]
  n <- nrow(dat)
  r_vec <- as.numeric(dat$r)
  y_vec <- dat$y

  # Use package to fit
  ebmr <- EPS$new("y", subset_ps_spec, dat, W_func)
  fit <- ebmr$ps_fit.list[[1]]
  alpha <- fit$coefficients
  design_matrix <- fit$design_matrix
  h_x <- fit$h_x
  esteq_dim <- ncol(h_x)

  g_func <- function(param) {
    eta <- design_matrix %*% param
    pi_vec <- compute_pi(eta)
    (r_vec / as.vector(pi_vec) - 1) * h_x
  }
  G_func <- function(param) matrix(colMeans(g_func(param)), esteq_dim, 1)

  # Reconstruct step 1 and step 2
  param_dim <- ncol(design_matrix)

  # Step 1: W = I
  obj_I <- function(param) {
    G_bar <- G_func(param)
    as.numeric(crossprod(G_bar))
  }
  opt1 <- optim(rep(0, param_dim), obj_I, method = "L-BFGS-B",
                lower = rep(-Inf, param_dim), upper = rep(Inf, param_dim),
                control = list(maxit = 1000))
  alpha_s1 <- opt1$par
  Q_s1 <- opt1$value

  # Step 2: W = (g'g/n)^-1 at step1
  g_s1 <- g_func(alpha_s1)
  W_s1 <- tryCatch(solve(crossprod(g_s1) / n), error = function(e) diag(esteq_dim))

  obj_s2 <- function(param) {
    G_bar <- G_func(param)
    as.numeric(crossprod(G_bar, W_s1 %*% G_bar))
  }
  opt2 <- optim(alpha_s1, obj_s2, method = "L-BFGS-B",
                lower = rep(-Inf, param_dim), upper = rep(Inf, param_dim),
                control = list(maxit = 1000))
  alpha_s2 <- opt2$par

  pi_s2 <- compute_pi(design_matrix %*% alpha_s2)
  mu_s2 <- mean(r_vec / as.vector(pi_s2) * y_vec)

  cat(sprintf("%-6d %-12.4f %-12.4f %-12.6e %-12.1f %-12.2e %-12.2e\n",
              rep_i, mu_s2, sqrt(sum(alpha_s2^2)), min(pi_s2),
              max(r_vec / as.vector(pi_s2)), Q_s1, opt2$value))
}

# Deep dive on rep 507
cat("\n=== DEEP DIVE: Rep 507 ===\n")
rep_i <- 507
dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]
n <- nrow(dat)
r_vec <- as.numeric(dat$r)
y_vec <- dat$y

ebmr <- EPS$new("y", subset_ps_spec, dat, W_func)
fit <- ebmr$ps_fit.list[[1]]
design_matrix <- fit$design_matrix
h_x <- fit$h_x
esteq_dim <- ncol(h_x)
param_dim <- ncol(design_matrix)

g_func <- function(param) {
  eta <- design_matrix %*% param
  pi_vec <- compute_pi(eta)
  (r_vec / as.vector(pi_vec) - 1) * h_x
}
G_func <- function(param) matrix(colMeans(g_func(param)), esteq_dim, 1)

# Step 1
obj_I <- function(param) {
  G_bar <- G_func(param)
  as.numeric(crossprod(G_bar))
}
opt1 <- optim(rep(0, param_dim), obj_I, method = "L-BFGS-B",
              lower = rep(-Inf, param_dim), upper = rep(Inf, param_dim),
              control = list(maxit = 1000))
alpha_s1 <- opt1$par

cat(sprintf("Step 1 alpha: [%s]\n", paste(round(alpha_s1, 4), collapse=", ")))
cat(sprintf("Step 1 ||alpha||: %.4f\n", sqrt(sum(alpha_s1^2))))
cat(sprintf("Step 1 Q(W=I): %.2e\n", opt1$value))

pi_s1 <- compute_pi(design_matrix %*% alpha_s1)
cat(sprintf("Step 1 min(pi): %.6e\n", min(pi_s1)))
cat(sprintf("Step 1 mu_ipw: %.4f\n", mean(r_vec / as.vector(pi_s1) * y_vec)))

# Check W_s1
g_s1 <- g_func(alpha_s1)
ggn <- crossprod(g_s1) / n
eig_ggn <- eigen(ggn)$values
cat(sprintf("\ng'g/n eigenvalues: [%s]\n", paste(round(eig_ggn, 6), collapse=", ")))
cat(sprintf("Condition number of g'g/n: %.1f\n", max(eig_ggn) / min(eig_ggn)))

W_s1 <- solve(ggn)
cat(sprintf("diag(W_s1): [%s]\n", paste(round(diag(W_s1), 2), collapse=", ")))

# Step 2
obj_s2 <- function(param) {
  G_bar <- G_func(param)
  as.numeric(crossprod(G_bar, W_s1 %*% G_bar))
}

# Check objective landscape along path from alpha_s1 to alpha_s2
opt2 <- optim(alpha_s1, obj_s2, method = "L-BFGS-B",
              lower = rep(-Inf, param_dim), upper = rep(Inf, param_dim),
              control = list(maxit = 1000))
alpha_s2 <- opt2$par

cat(sprintf("\nStep 2 alpha: [%s]\n", paste(round(alpha_s2, 2), collapse=", ")))
cat(sprintf("Step 2 ||alpha||: %.4f\n", sqrt(sum(alpha_s2^2))))
cat(sprintf("Step 2 Q(W_s1): %.2e\n", opt2$value))
cat(sprintf("Step 2 convergence: %d\n", opt2$convergence))

pi_s2 <- compute_pi(design_matrix %*% alpha_s2)
cat(sprintf("Step 2 min(pi): %.6e\n", min(pi_s2)))
cat(sprintf("Step 2 mu_ipw: %.4f\n", mean(r_vec / as.vector(pi_s2) * y_vec)))
cat(sprintf("# obs with eta at bound: %d\n", sum(abs(design_matrix %*% alpha_s2) >= eta_max - 0.01)))

# What is G_bar at step2?
G_s2 <- as.vector(G_func(alpha_s2))
cat(sprintf("\nG_bar at step2: [%s]\n", paste(round(G_s2, 8), collapse=", ")))
cat(sprintf("||G_bar|| at step2: %.2e\n", sqrt(sum(G_s2^2))))

G_s1 <- as.vector(G_func(alpha_s1))
cat(sprintf("G_bar at step1: [%s]\n", paste(round(G_s1, 8), collapse=", ")))
cat(sprintf("||G_bar|| at step1: %.2e\n", sqrt(sum(G_s1^2))))

# Check: does step2 objective actually decrease relative to step1 at the same W?
cat(sprintf("\nQ(alpha_s1, W_s1) = %.6e\n", obj_s2(alpha_s1)))
cat(sprintf("Q(alpha_s2, W_s1) = %.6e\n", obj_s2(alpha_s2)))

# Hessian at step2
H_s2 <- numDeriv::hessian(obj_s2, alpha_s2)
eig_H <- eigen(H_s2)$values
cat(sprintf("\nHessian eig at step2: [%s]\n", paste(round(eig_H, 6), collapse=", ")))
cat(sprintf("All positive: %s\n", all(eig_H > 0)))

cat("\nDone!\n")
