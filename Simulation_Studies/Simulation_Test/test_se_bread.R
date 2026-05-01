setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
source("Simulation.r")
library(EBMRalgorithmFast4)
library(numDeriv)

# Setting 3, model 2 alone
data_file <- misspecified_model_all_data_file.list[["setting3"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
n_val <- 2000
mu_true <- get_mu_true("setting3")

ps_spec <- get_ps_spec("9-alt1")
single_ps_spec <- list(
  formula.list = ps_spec[["formula.list"]][2],
  h_alpha.list = ps_spec[["h_alpha.list"]][2],
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]]
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

# Check a few reps: is the bread matrix (A) missing a D_W = Gamma'(dW/dalpha)G term?
test_reps <- c(1, 2, 3, 5, 10, 20)

for (rep_i in test_reps) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]

  ebmr <- EBMRAlgorithmFast4[["new"]]("y", single_ps_spec, dat, W_func)
  res <- ebmr[["EBMR_IPW"]](
    h_nu = function(dat) cbind(u1 = dat[["u1"]], u2 = dat[["u2"]], z1 = dat[["z1"]], z2 = dat[["z2"]]),
    se.fit = TRUE, type = "HT"
  )

  gmm_fit <- ebmr[["ps_fit.list"]][[1]][["gmm_fit"]]
  alpha_hat <- gmm_fit[["estimates"]]
  Gamma.hat <- gmm_fit[["Gamma.hat"]]
  W.hat <- gmm_fit[["W.hat"]]
  g.matrix <- gmm_fit[["g.matrix"]]
  eta_s <- gmm_fit[["eta_s"]]
  Q_full <- gmm_fit[["Q"]]

  n <- nrow(g.matrix)
  esteq_dim <- ncol(g.matrix)
  param_dim <- length(alpha_hat)

  # Code's bread: A_code = -(Gamma'W Gamma + M_n) = Q^{-1}
  # So Q = -A_code^{-1}, meaning A_code = -Q^{-1}
  # Q = (I - H_0n M_n)^{-1} H_0n
  # A_code = Gamma'WGamma + M_n

  GtW <- crossprod(Gamma.hat, W.hat)
  GtWG <- GtW %*% Gamma.hat
  A_code <- GtWG  # This is just the GtWG part. The full A_code = GtWG + M_n
  # Actually A_code = -Q^{-1}
  Q_inv <- solve(Q_full)
  A_code_full <- -Q_inv  # This is Gamma'WGamma + M_n

  # Compute the TRUE bread numerically: A_true = d(Gamma'WG)/dalpha
  # f(alpha) = Gamma(alpha)' W(alpha) G(alpha)
  # This is the FOC that is solved

  # Reconstruct g function from the ps_fit
  ps_fit <- ebmr[["ps_fit.list"]][[1]]
  design_matrix <- ps_fit[["design_matrix"]]
  link_type <- ps_fit[["link_type"]]
  r_vec <- dat[["r"]]

  h_alpha_names <- single_ps_spec[["h_alpha.list"]][[1]]
  h_x_raw <- as.matrix(dat[, h_alpha_names, drop = FALSE])
  h_x <- cbind(1, h_x_raw)
  h_dim <- ncol(h_x)

  inv_link_func <- single_ps_spec[["inv_link"]]

  g_func <- function(alpha) {
    eta <- as.vector(design_matrix %*% alpha)
    eta <- pmin(pmax(eta, -20), 20)
    pi_vec <- inv_link_func(eta)
    (r_vec / pi_vec - 1) * h_x
  }

  G_func <- function(alpha) {
    g_mat <- g_func(alpha)
    colMeans(g_mat)
  }

  # f(alpha) = Gamma(alpha)' W(alpha) G(alpha)
  # where Gamma = dG/dalpha, W = (g'g/n)^{-1}
  f_func <- function(alpha) {
    g_mat <- g_func(alpha)
    G_vec <- colMeans(g_mat)
    W_mat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(h_dim))
    Gamma_mat <- jacobian(G_func, alpha)  # esteq_dim x param_dim
    as.vector(t(Gamma_mat) %*% W_mat %*% G_vec)
  }

  # Numerical Jacobian of f w.r.t. alpha: this is the TRUE bread
  A_numerical <- jacobian(f_func, alpha_hat)  # param_dim x param_dim

  cat(sprintf("\n=== Rep %d ===\n", rep_i))
  cat(sprintf("||eta_s|| = %.6f\n", sqrt(sum(eta_s^2))))
  cat(sprintf("||alpha|| = %.4f\n", sqrt(sum(alpha_hat^2))))

  cat("\nA_code (Gamma'WGamma + M_n):\n")
  print(round(A_code_full, 4))

  cat("\nA_numerical (true Jacobian):\n")
  print(round(A_numerical, 4))

  cat("\nDifference (D_W = A_numerical - A_code):\n")
  D_W <- A_numerical - A_code_full
  print(round(D_W, 4))
  cat(sprintf("||D_W|| / ||A_code|| = %.4f\n", norm(D_W, "F") / norm(A_code_full, "F")))

  # Now compute SE with corrected bread
  A_corrected <- A_numerical
  Q_corrected <- -solve(A_corrected)

  # Extract H_sum = H_1n + H_2n_2 + H_2n_3
  psi_alpha <- gmm_fit[["psi"]]
  H_sum <- Q_inv %*% psi_alpha

  # Corrected psi
  psi_corrected <- Q_corrected %*% H_sum

  # Also compute psi without H_2n_3
  H_1n <- GtW %*% (t(g.matrix) - eta_s)
  W_eta <- as.vector(W.hat %*% eta_s)
  tGamma_W_eta <- as.vector(GtW %*% eta_s)
  GtW_g <- GtW %*% t(g.matrix)
  g_W_eta <- as.vector(g.matrix %*% W_eta)
  H_2n_3 <- matrix(tGamma_W_eta, param_dim, n) - GtW_g * matrix(g_W_eta, param_dim, n, byrow = TRUE)
  H_2n_2 <- H_sum - H_1n - H_2n_3

  # Corrected bread, no H_2n_3
  psi_corrected_noH2n3 <- Q_corrected %*% (H_1n + H_2n_2)

  # Compute mu SE for each variant
  pi_hat <- ps_fit[["fitted.values"]]
  y_vec <- dat[["y"]]
  base_term <- as.vector(r_vec / pi_hat * y_vec)

  ry_ps_inv2 <- as.vector(r_vec * y_vec * (pi_hat^(-2)))
  if (link_type == "logistic_complement") {
    dot_pi <- -design_matrix * (pi_hat * (1 - pi_hat))
  } else {
    dot_pi <- design_matrix * (pi_hat * (1 - pi_hat))
  }
  H_alpha.w <- colMeans(dot_pi * ry_ps_inv2)

  # Original SE
  mu_iid_orig <- base_term - as.vector(t(H_alpha.w) %*% psi_alpha)
  se_orig <- sqrt(var(mu_iid_orig) / n)

  # Corrected bread SE
  mu_iid_corrected <- base_term - as.vector(t(H_alpha.w) %*% psi_corrected)
  se_corrected <- sqrt(var(mu_iid_corrected) / n)

  # Corrected bread, no H_2n_3
  mu_iid_corr_noH2n3 <- base_term - as.vector(t(H_alpha.w) %*% psi_corrected_noH2n3)
  se_corr_noH2n3 <- sqrt(var(mu_iid_corr_noH2n3) / n)

  # No H_2n_3, original bread
  psi_noH2n3 <- Q_full %*% (H_1n + H_2n_2)
  mu_iid_noH2n3 <- base_term - as.vector(t(H_alpha.w) %*% psi_noH2n3)
  se_noH2n3 <- sqrt(var(mu_iid_noH2n3) / n)

  cat(sprintf("\n  se_orig         = %.4f\n", se_orig))
  cat(sprintf("  se_corrected    = %.4f  (corrected bread)\n", se_corrected))
  cat(sprintf("  se_noH2n3       = %.4f  (no H_2n_3, orig bread)\n", se_noH2n3))
  cat(sprintf("  se_corr_noH2n3  = %.4f  (corrected bread, no H_2n_3)\n", se_corr_noH2n3))
}
