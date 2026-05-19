## Test convergence criterion: |alpha_{t} - alpha_{t-1}| < 1e-6
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

ps_spec_base <- get_ps_spec("7")
h_alpha_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2)

data_file <- correct_model_all_data_file.list[["setting3"]][["miss50"]]
all_data <- NULL
for (fi in seq_along(data_file)) {
  if (file.exists(data_file[[fi]])) {
    test <- readRDS(data_file[[fi]])
    if (nrow(test) / 1000 == 2000) { all_data <- test; break }
  }
}

nn <- 2000
design_formula <- ps_spec_base$formula.list[[2]]
compute_pi <- function(eta) 1/(1+exp(eta))

cat("=== Convergence criterion: |alpha_t - alpha_{t-1}| < 1e-6 ===\n\n")

for (rep_i in c(1, 3, 7, 13)) {
  dat <- all_data[((rep_i-1)*nn + 1):(rep_i*nn), ]
  design_mat <- model.matrix(design_formula, dat)
  r_vec <- as.vector(dat$r)
  n <- nn
  h_x <- cbind(1, h_alpha_fn(dat))
  h_dim <- ncol(h_x)
  alpha_dim <- ncol(design_mat)

  g_fn <- function(alpha) {
    eta <- pmin(pmax(as.vector(design_mat %*% alpha), -20), 20)
    pi_vec <- compute_pi(eta)
    as.vector(r_vec/pi_vec - 1) * h_x
  }
  G_fn <- function(alpha) colMeans(g_fn(alpha))
  Gamma_fn <- function(alpha) {
    eta <- pmin(pmax(as.vector(design_mat %*% alpha), -20), 20)
    pi_vec <- compute_pi(eta)
    common <- r_vec * (1 - pi_vec) / pi_vec
    Gamma_mat <- matrix(0, h_dim, alpha_dim)
    for (j in 1:alpha_dim) Gamma_mat[, j] <- colMeans(common * design_mat[, j] * h_x)
    Gamma_mat
  }

  alpha <- rep(0, alpha_dim)
  # Step 1: W=I
  obj_I <- function(a) { G <- G_fn(a); sum(G^2) }
  gr_I <- function(a) { G <- G_fn(a); 2 * as.vector(t(Gamma_fn(a)) %*% G) }
  opt1 <- optim(alpha, obj_I, gr_I, method = "L-BFGS-B", control = list(maxit = 1000))
  alpha <- opt1$par

  cat(sprintf("--- Rep %d ---\n", rep_i))
  for (t in 1:500) {
    alpha_old <- alpha
    g_mat <- g_fn(alpha)
    W_hat <- tryCatch(solve(t(g_mat) %*% g_mat / n), error = function(e) diag(h_dim))
    obj_fn <- function(a) { G <- matrix(G_fn(a), ncol=1); as.numeric(t(G) %*% W_hat %*% G) }
    gr_fn <- function(a) { G <- G_fn(a); 2 * as.vector(t(Gamma_fn(a)) %*% W_hat %*% G) }
    opt_t <- optim(alpha, obj_fn, gr_fn, method = "L-BFGS-B", control = list(maxit = 1000))
    alpha <- opt_t$par

    alpha_change <- max(abs(alpha - alpha_old))
    G_vec <- G_fn(alpha)
    grad_norm <- max(abs(2 * as.vector(t(Gamma_fn(alpha)) %*% W_hat %*% G_vec)))

    if (t <= 10 || t %% 20 == 0 || alpha_change < 1e-6) {
      cat(sprintf("  iter %3d: |dα|=%.4e  |2Γ'WG|=%.4e  alpha=(%s)\n",
          t, alpha_change, grad_norm, paste(round(alpha, 4), collapse=",")))
    }
    if (alpha_change < 1e-6) { cat("  *** ALPHA CONVERGED ***\n"); break }
  }
  if (t == 500) cat("  Did not converge after 500 iterations\n")
  cat("\n")
}
