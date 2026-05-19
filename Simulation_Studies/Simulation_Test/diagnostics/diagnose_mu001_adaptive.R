## Adaptive damping — the most promising strategy
## Run with more iterations and track convergence
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

inner_opt <- function(alpha, W_hat) {
  obj <- function(a) {
    G <- matrix(G_fn(a), ncol=1)
    as.numeric(t(G) %*% W_hat %*% G)
  }
  optim(alpha, obj, method = "L-BFGS-B", control = list(maxit = 1000))$par
}

cat("=== Adaptive damping with 500 iterations ===\n\n")

alpha <- rep(0, alpha_dim)
lambda <- 0.0
W_prev <- diag(h_dim)
prev_gn <- Inf

for (t in 1:500) {
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

  if (t <= 20 || t %% 20 == 0 || gn < 1e-6) {
    cat(sprintf("  iter %3d: grad=%.4e, max|G|=%.4e, lambda=%.2f\n",
        t, gn, max(abs(G_fn(alpha))), lambda))
  }
  if (gn < 1e-6) { cat("  *** CONVERGED ***\n"); break }
}

cat(sprintf("\nFinal: alpha=(%s)\n", paste(round(alpha, 4), collapse=", ")))
cat(sprintf("Final max|G|=%.4e\n", max(abs(G_fn(alpha)))))

# Now check if this solution gives reasonable SE
cat("\n=== SE at converged solution ===\n\n")
ps_spec <- list(
  formula.list = list(ps_spec_base$formula.list[[3]]),
  h_alpha.list = list(h_alpha_fn),
  inv_link = ps_spec_base$inv_link,
  outcome = ps_spec_base$outcome,
  alpha_init.list = list(alpha)
)
W_fn <- function(g.matrix) solve(var(g.matrix))
ebmr <- EPS$new("y", ps_spec, dat, W_fn)
ipw <- ebmr$EBMR_IPW(h_alpha_fn, model_indices = 1, se.fit = TRUE)
cat(sprintf("mu_ipw=%.6f, se=%.6f\n", ipw$mu_ipw, ipw$se_ipw))
cat(sprintf("GMM converged: %s, grad_norm: %.2e\n",
    ebmr$ps_fit.list[[1]]$gmm_fit$opt$converged,
    ebmr$ps_fit.list[[1]]$gmm_fit$opt$final_grad_norm))
