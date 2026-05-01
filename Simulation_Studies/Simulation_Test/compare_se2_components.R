setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

n_val <- 2000
ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

model_idx <- 3

analyze_rep <- function(rep_i) {
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
  r_vec <- dat[["r"]]; y_vec <- dat[["y"]]; n <- n_val

  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL),
    optimizer = "constrained_nr"
  )
  ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
  gmm_fit <- ebmr$ps_fit.list[[1]]$gmm_fit
  pi_hat <- ebmr$ps_fit.list[[1]]$fitted.values
  design_mat <- ebmr$ps_fit.list[[1]]$design_matrix
  link_type <- ebmr$ps_fit.list[[1]]$link_type
  psi2_mat <- gmm_fit$psi2

  mu_hat <- mean(r_vec * y_vec / pi_hat)
  sol_cond <- gmm_fit$opt$solution_cond

  cat(sprintf("\n--- Rep %d ---\n", rep_i))
  cat(sprintf("  mu_hat=%.4f, cond=%.2e, grad=%.2e\n",
      mu_hat, sol_cond, gmm_fit$opt$final_grad_norm))
  cat(sprintf("  alpha = %s\n", paste(round(gmm_fit$estimates, 4), collapse=", ")))
  cat(sprintf("  pi: min=%.6f, max=%.6f, mean=%.4f\n",
      min(pi_hat), max(pi_hat), mean(pi_hat)))
  cat(sprintf("  # pi > 0.99: %d, # pi > 0.999: %d\n",
      sum(pi_hat > 0.99), sum(pi_hat > 0.999)))

  # Raw IPW term: r/pi * y
  raw_ipw <- r_vec * y_vec / pi_hat
  cat(sprintf("  var(r/pi*y) = %.4f, mean(r/pi*y) = %.4f\n",
      var(raw_ipw), mean(raw_ipw)))

  if (!is.null(psi2_mat) && is.matrix(psi2_mat) && !any(is.na(psi2_mat))) {
    # H_alpha
    if (!is.null(link_type) && link_type == "logistic_complement") {
      dot_pi <- -design_mat * (pi_hat * (1 - pi_hat))
    } else {
      dot_pi <- design_mat * (pi_hat * (1 - pi_hat))
    }
    ry_ps_inv2 <- as.vector(r_vec * y_vec * (pi_hat^(-2)))
    H_alpha <- colMeans(dot_pi * ry_ps_inv2)
    cat(sprintf("  H_alpha = %s\n", paste(round(H_alpha, 4), collapse=", ")))
    cat(sprintf("  ||H_alpha|| = %.4f\n", sqrt(sum(H_alpha^2))))

    # Adjustment term: H_alpha' * psi2
    adjust <- as.vector(t(H_alpha) %*% psi2_mat)
    cat(sprintf("  var(adjustment) = %.6f\n", var(adjust)))
    cat(sprintf("  mean(adjustment) = %.6f\n", mean(adjust)))

    # Full mu_iid
    mu_iid <- raw_ipw - adjust
    se2 <- sqrt(var(mu_iid) / n)
    cat(sprintf("  var(mu_iid) = %.6f\n", var(mu_iid)))
    cat(sprintf("  SE2 = %.4f\n", se2))

    # Correlation between raw_ipw and adjustment
    cat(sprintf("  cor(raw_ipw, adjust) = %.4f\n", cor(raw_ipw, adjust)))

    # Variance decomposition
    cat(sprintf("  var(raw) = %.4f, var(adj) = %.6f, 2*cov = %.6f\n",
        var(raw_ipw), var(adjust), -2*cov(raw_ipw, adjust)))
    cat(sprintf("  var(mu_iid) = var(raw) + var(adj) - 2*cov = %.6f\n",
        var(raw_ipw) + var(adjust) - 2*cov(raw_ipw, adjust)))

    # psi2 stats
    cat(sprintf("  var(psi2) per row: %s\n",
        paste(round(apply(psi2_mat, 1, var), 4), collapse=", ")))
    cat(sprintf("  ||psi2||_F = %.4f\n", sqrt(sum(psi2_mat^2))))

    # H2 eigenvalues
    Gamma_hat <- gmm_fit$Gamma.hat
    W_hat <- gmm_fit$W.hat
    GtWG <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
    eig_GtWG <- eigen(GtWG, symmetric = TRUE)$values
    cat(sprintf("  eig(GtWG) = %s\n", paste(round(eig_GtWG, 6), collapse=", ")))
  } else {
    cat("  psi2 is NA\n")
  }
}

cat("=== Normal reps ===")
for (r in c(1, 10, 50)) analyze_rep(r)

cat("\n\n=== Degenerate reps ===")
for (r in c(35, 68, 166)) analyze_rep(r)
