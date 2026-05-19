## Check: is the full dQ/dα = 0 at L-BFGS-B solution?
## Compare 2Γ'WG vs full numerical gradient of Q(α) = G(α)'W(α)G(α)
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EPS", quiet = TRUE)
library(numDeriv)

ps_spec_base <- get_ps_spec("7")
h_alpha_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2)
W_sm <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

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

ps_spec_lb <- list(formula.list = list(design_formula), h_alpha.list = list(h_alpha_fn),
                   inv_link = ps_spec_base$inv_link, outcome = ps_spec_base$outcome)

cat("=== Full dQ/dα vs 2Γ'WG at L-BFGS-B solution ===\n\n")

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

  # Full CUE objective: Q(α) = G(α)' W(α) G(α) where W(α) = solve(g'g/n)
  Q_cue <- function(alpha) {
    g_mat <- g_fn(alpha)
    G <- colMeans(g_mat)
    W <- tryCatch(solve(t(g_mat) %*% g_mat / n), error = function(e) diag(h_dim))
    as.numeric(t(G) %*% W %*% G)
  }

  # Get L-BFGS-B solution
  ebmr <- EPS$new("y", ps_spec_lb, dat, W_sm)
  alpha_lb <- ebmr$ps_fit.list[[1]]$coefficients

  # Compute 2Γ'WG (iterative criterion)
  g_mat <- g_fn(alpha_lb)
  G_vec <- G_fn(alpha_lb)
  W_hat <- solve(t(g_mat) %*% g_mat / n)
  Gamma_hat <- Gamma_fn(alpha_lb)
  iter_grad <- 2 * as.vector(t(Gamma_hat) %*% W_hat %*% G_vec)

  # Compute full dQ/dα numerically
  full_grad <- numDeriv::grad(Q_cue, alpha_lb)

  cat(sprintf("--- Rep %d ---\n", rep_i))
  cat(sprintf("  alpha: (%s)\n", paste(round(alpha_lb, 4), collapse=", ")))
  cat(sprintf("  max|2Γ'WG|:     %.4e (iterative criterion)\n", max(abs(iter_grad))))
  cat(sprintf("  max|full dQ/dα|: %.4e (numerical)\n", max(abs(full_grad))))
  cat(sprintf("  2Γ'WG:     (%s)\n", paste(round(iter_grad, 6), collapse=", ")))
  cat(sprintf("  full dQ/dα: (%s)\n", paste(round(full_grad, 6), collapse=", ")))
  cat(sprintf("  difference: (%s)\n\n", paste(round(full_grad - iter_grad, 6), collapse=", ")))
}
