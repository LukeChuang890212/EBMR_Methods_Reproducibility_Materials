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

# Check multiple outlier AND normal reps
test_reps <- c(929, 438, 430, 221, 180, 1, 2, 3, 4, 5)

cat(sprintf("%-6s %-5s %-12s %-12s %-12s %-12s %-12s\n",
            "Rep", "Type", "eig1", "eig2", "eig3", "eig4", "cond_num"))

for (rep_i in test_reps) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]
  label <- if (rep_i %in% c(929, 438, 430, 221, 180)) "OUT" else "NORM"

  ebmr <- EBMRAlgorithmFast4$new("y", subset_ps_spec, dat, W_func)
  fit <- ebmr$ps_fit.list[[1]]
  alpha <- fit$coefficients
  h_x <- fit$h_x
  design_matrix <- fit$design_matrix
  r_vec <- as.numeric(dat$r)
  n <- nrow(dat)
  esteq_dim <- ncol(h_x)

  g_func <- function(param) {
    eta <- design_matrix %*% param
    pi_vec <- compute_pi(eta)
    (r_vec / as.vector(pi_vec) - 1) * h_x
  }
  G_func <- function(param) matrix(colMeans(g_func(param)), esteq_dim, 1)

  g_mat <- g_func(alpha)
  W_hat <- solve(crossprod(g_mat) / n)

  H <- numDeriv::hessian(function(p) {
    G_bar <- G_func(p)
    as.numeric(crossprod(G_bar, W_hat %*% G_bar))
  }, alpha)

  eig <- sort(eigen(H)$values)
  cond <- max(eig) / min(abs(eig))

  cat(sprintf("%-6d %-5s %-12.6f %-12.6f %-12.6f %-12.6f %-12.1f\n",
              rep_i, label, eig[1], eig[2], eig[3], eig[4], cond))
}

cat("\n")

# Now for rep 929: perturb alpha along the smallest eigenvector and see how Q changes
cat("=== Rep 929: Objective along smallest eigenvector ===\n")
rep_i <- 929
dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]
ebmr <- EBMRAlgorithmFast4$new("y", subset_ps_spec, dat, W_func)
fit <- ebmr$ps_fit.list[[1]]
alpha <- fit$coefficients
h_x <- fit$h_x
design_matrix <- fit$design_matrix
r_vec <- as.numeric(dat$r)
n <- nrow(dat)
esteq_dim <- ncol(h_x)

g_func <- function(param) {
  eta <- design_matrix %*% param
  pi_vec <- compute_pi(eta)
  (r_vec / as.vector(pi_vec) - 1) * h_x
}
G_func <- function(param) matrix(colMeans(g_func(param)), esteq_dim, 1)

g_mat <- g_func(alpha)
W_hat <- solve(crossprod(g_mat) / n)

H <- numDeriv::hessian(function(p) {
  G_bar <- G_func(p)
  as.numeric(crossprod(G_bar, W_hat %*% G_bar))
}, alpha)

eig_decomp <- eigen(H)
v_min <- eig_decomp$vectors[, which.min(eig_decomp$values)]
cat(sprintf("Smallest eigenvector: [%s]\n", paste(round(v_min, 4), collapse=", ")))
cat(sprintf("Smallest eigenvalue: %.6f\n\n", min(eig_decomp$values)))

# Perturb along v_min
cat(sprintf("%-10s %-12s %-12s\n", "delta", "Q(alpha+d*v)", "mu_ipw"))
for (delta in c(-2, -1, -0.5, -0.1, 0, 0.1, 0.5, 1, 2)) {
  alpha_pert <- alpha + delta * v_min
  G_bar <- G_func(alpha_pert)
  Q_val <- as.numeric(crossprod(G_bar, W_hat %*% G_bar))
  pi_pert <- as.vector(compute_pi(design_matrix %*% alpha_pert))
  mu_pert <- mean(r_vec / pi_pert * dat$y)
  cat(sprintf("%-10.1f %-12.6e %-12.4f\n", delta, Q_val, mu_pert))
}

cat("\nDone!\n")
