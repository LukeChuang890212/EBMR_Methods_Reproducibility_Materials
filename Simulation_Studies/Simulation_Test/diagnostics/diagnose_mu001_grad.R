## Check gradient of CUE objective at different solutions
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

cue_obj <- function(alpha) {
  g_mat <- g_fn(alpha)
  G <- colMeans(g_mat)
  W <- tryCatch(solve(var(g_mat)), error = function(e) diag(h_dim))
  as.numeric(t(G) %*% W %*% G)
}

# Solutions to compare
alpha_degenerate <- c(2.6178, -41.9016, 9.401, 5.8311)  # Nelder-Mead
alpha_damped <- c(2.2883, -5.1034, 1.2044, 0.5302)      # Damped CUE
alpha_iterative <- c(2.332, -5.7017, 1.2359, 0.593)     # Package iterative (stalled)

cat("=== Gradient of CUE objective at different solutions ===\n\n")

for (name in c("degenerate", "damped", "iterative")) {
  alpha <- switch(name,
    degenerate = alpha_degenerate,
    damped = alpha_damped,
    iterative = alpha_iterative)

  obj <- cue_obj(alpha)
  grad <- numDeriv::grad(cue_obj, alpha)
  G <- G_fn(alpha)
  g_mat <- g_fn(alpha)
  ps <- compute_pi(pmin(pmax(as.vector(design_mat %*% alpha), -20), 20))

  cat(sprintf("--- %s ---\n", name))
  cat(sprintf("  alpha = (%s)\n", paste(round(alpha, 4), collapse=", ")))
  cat(sprintf("  CUE obj = %.6e\n", obj))
  cat(sprintf("  max|grad(Q)| = %.6e\n", max(abs(grad))))
  cat(sprintf("  norm(grad(Q)) = %.6e\n", sqrt(sum(grad^2))))
  cat(sprintf("  max|G| = %.6e\n", max(abs(G))))
  cat(sprintf("  PS: min=%.4f, Q1=%.4f, med=%.4f, Q3=%.4f, max=%.4f\n",
      min(ps), quantile(ps,0.25), median(ps), quantile(ps,0.75), max(ps)))
  cat(sprintf("  cond(var(g)) = %.2e\n\n", kappa(var(g_mat))))
}
