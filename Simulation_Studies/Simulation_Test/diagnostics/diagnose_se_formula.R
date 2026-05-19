setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

n_val <- 2000
n_reps <- 200

ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
mu_true <- get_mu_true("setting4")

model_idx <- 3

single_ps <- list(
  formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
  h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]],
  alpha_init.list = list(NULL),
  optimizer = "constrained_nr"
)
W <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

alpha_dim <- 4  # from ps formula
alpha_mat <- matrix(NA_real_, n_reps, alpha_dim)
mu_vals <- rep(NA_real_, n_reps)

# Package SE for alpha
se1_alpha_mat <- matrix(NA_real_, n_reps, alpha_dim)
se2_alpha_mat <- matrix(NA_real_, n_reps, alpha_dim)

# Simple SE for alpha: diag of (Gamma'WGamma)^-1 Gamma'W Omega W Gamma (Gamma'WGamma)^-1 / n
# When W = Omega^{-1}, simplifies to (Gamma'WGamma)^{-1}/n
se_simple_alpha_mat <- matrix(NA_real_, n_reps, alpha_dim)

# Track key diagnostics
eta_norm_vals <- rep(NA_real_, n_reps)
Mn_norm_vals <- rep(NA_real_, n_reps)
Q_norm_vals <- rep(NA_real_, n_reps)
H0n_norm_vals <- rep(NA_real_, n_reps)

for (rep_i in 1:n_reps) {
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W)

    r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
    pi_hat <- ebmr$ps_fit.list[[1]]$fitted.values
    mu_hat <- mean(r_vec * y_vec / pi_hat)

    gmm_fit <- ebmr$ps_fit.list[[1]]$gmm_fit
    alpha_mat[rep_i, ] <- gmm_fit$estimates
    mu_vals[rep_i] <- mu_hat
    se1_alpha_mat[rep_i, ] <- gmm_fit$se1
    se2_alpha_mat[rep_i, ] <- gmm_fit$se2

    # Simple SE: sandwich
    Gamma_hat <- gmm_fit$Gamma.hat
    g_mat <- gmm_fit$g.matrix
    W_hat <- gmm_fit$W.hat
    GtW <- crossprod(Gamma_hat, W_hat)
    GtWG <- GtW %*% Gamma_hat
    H_inv <- tryCatch(solve(GtWG), error = function(e) NULL)
    if (!is.null(H_inv)) {
      # Full sandwich: H_inv %*% GtW %*% Omega %*% t(GtW) %*% H_inv / n
      Omega <- crossprod(g_mat) / n_val
      sandwich <- H_inv %*% GtW %*% Omega %*% t(GtW) %*% H_inv / n_val
      se_simple_alpha_mat[rep_i, ] <- sqrt(diag(sandwich))
    }

    # Track diagnostics from gmm_fit internals
    eta_norm_vals[rep_i] <- sqrt(sum(gmm_fit$eta_s^2))

  }, error = function(e) NULL)
}

valid <- !is.na(mu_vals)

cat("=== Alpha parameter SE comparison (setting4, model3, constrained_nr) ===\n\n")
cat(sprintf("Valid reps: %d/%d\n\n", sum(valid), n_reps))

# ESD for each alpha parameter
alpha_esd <- apply(alpha_mat[valid, ], 2, sd)
cat("Alpha ESD across reps:\n")
cat(sprintf("  %s\n", paste(sprintf("%.6f", alpha_esd), collapse=", ")))

# ESE1/ESD for each alpha
ese1_alpha <- colMeans(se1_alpha_mat[valid, ])
cat("\nAlpha SE1 (package): ESE/ESD per param:\n")
cat(sprintf("  %s\n", paste(sprintf("%.3f", ese1_alpha/alpha_esd), collapse=", ")))

# ESE2/ESD for each alpha
v2 <- valid & apply(!is.na(se2_alpha_mat), 1, all)
if (sum(v2) > 0) {
  ese2_alpha <- colMeans(se2_alpha_mat[v2, ])
  alpha_esd2 <- apply(alpha_mat[v2, ], 2, sd)
  cat(sprintf("\nAlpha SE2 (package): ESE/ESD per param (valid=%d):\n", sum(v2)))
  cat(sprintf("  %s\n", paste(sprintf("%.3f", ese2_alpha/alpha_esd2), collapse=", ")))
}

# Simple SE
vs <- valid & apply(!is.na(se_simple_alpha_mat), 1, all)
if (sum(vs) > 0) {
  ese_simple <- colMeans(se_simple_alpha_mat[vs, ])
  alpha_esds <- apply(alpha_mat[vs, ], 2, sd)
  cat(sprintf("\nAlpha SE (simple sandwich): ESE/ESD per param (valid=%d):\n", sum(vs)))
  cat(sprintf("  %s\n", paste(sprintf("%.3f", ese_simple/alpha_esds), collapse=", ")))
}

# Eta norm distribution
cat(sprintf("\n||eta_s|| distribution: median=%.4e, mean=%.4e, max=%.4e\n",
    median(eta_norm_vals[valid]), mean(eta_norm_vals[valid]), max(eta_norm_vals[valid])))

# Now do the same for model 2
cat("\n\n=== Now checking model 2 with L-BFGS-B ===\n\n")

single_ps2 <- list(
  formula.list = list(ps_spec[["formula.list"]][[2]]),
  h_alpha.list = list(ps_spec[["h_alpha.list"]][[2]]),
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]],
  alpha_init.list = list(NULL),
  optimizer = "L-BFGS-B"
)

alpha_mat2 <- matrix(NA_real_, n_reps, alpha_dim)
se1_alpha_mat2 <- matrix(NA_real_, n_reps, alpha_dim)
se2_alpha_mat2 <- matrix(NA_real_, n_reps, alpha_dim)
se_simple_alpha_mat2 <- matrix(NA_real_, n_reps, alpha_dim)
eta_norm_vals2 <- rep(NA_real_, n_reps)

for (rep_i in 1:n_reps) {
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", single_ps2, dat, W)
    gmm_fit <- ebmr$ps_fit.list[[1]]$gmm_fit
    alpha_mat2[rep_i, ] <- gmm_fit$estimates
    se1_alpha_mat2[rep_i, ] <- gmm_fit$se1
    se2_alpha_mat2[rep_i, ] <- gmm_fit$se2
    eta_norm_vals2[rep_i] <- sqrt(sum(gmm_fit$eta_s^2))

    Gamma_hat <- gmm_fit$Gamma.hat
    g_mat <- gmm_fit$g.matrix
    W_hat <- gmm_fit$W.hat
    GtW <- crossprod(Gamma_hat, W_hat)
    GtWG <- GtW %*% Gamma_hat
    H_inv <- tryCatch(solve(GtWG), error = function(e) NULL)
    if (!is.null(H_inv)) {
      Omega <- crossprod(g_mat) / n_val
      sandwich <- H_inv %*% GtW %*% Omega %*% t(GtW) %*% H_inv / n_val
      se_simple_alpha_mat2[rep_i, ] <- sqrt(diag(sandwich))
    }
  }, error = function(e) NULL)
}

valid2 <- !is.na(alpha_mat2[,1])
alpha_esd2 <- apply(alpha_mat2[valid2, ], 2, sd)
cat(sprintf("Valid reps: %d/%d\n\n", sum(valid2), n_reps))
cat("Alpha ESD across reps:\n")
cat(sprintf("  %s\n", paste(sprintf("%.6f", alpha_esd2), collapse=", ")))

ese1_a2 <- colMeans(se1_alpha_mat2[valid2, ])
cat("\nAlpha SE1 (package): ESE/ESD per param:\n")
cat(sprintf("  %s\n", paste(sprintf("%.3f", ese1_a2/alpha_esd2), collapse=", ")))

v2b <- valid2 & apply(!is.na(se2_alpha_mat2), 1, all)
if (sum(v2b) > 0) {
  ese2_a2 <- colMeans(se2_alpha_mat2[v2b, ])
  alpha_esd2b <- apply(alpha_mat2[v2b, ], 2, sd)
  cat(sprintf("\nAlpha SE2 (package): ESE/ESD per param (valid=%d):\n", sum(v2b)))
  cat(sprintf("  %s\n", paste(sprintf("%.3f", ese2_a2/alpha_esd2b), collapse=", ")))
}

vs2 <- valid2 & apply(!is.na(se_simple_alpha_mat2), 1, all)
if (sum(vs2) > 0) {
  ese_s2 <- colMeans(se_simple_alpha_mat2[vs2, ])
  alpha_esds2 <- apply(alpha_mat2[vs2, ], 2, sd)
  cat(sprintf("\nAlpha SE (simple sandwich): ESE/ESD per param (valid=%d):\n", sum(vs2)))
  cat(sprintf("  %s\n", paste(sprintf("%.3f", ese_s2/alpha_esds2), collapse=", ")))
}

cat(sprintf("\n||eta_s|| distribution: median=%.4e, mean=%.4e, max=%.4e\n",
    median(eta_norm_vals2[valid2]), mean(eta_norm_vals2[valid2]), max(eta_norm_vals2[valid2])))
