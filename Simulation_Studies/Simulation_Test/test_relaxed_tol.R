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
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

# Run constrained_nr (package default tol=1e-6), then externally decide
# whether to accept the CNR result or fall back to L-BFGS-B based on
# a relaxed tolerance threshold.
run_with_relaxed_tol <- function(model_idx, accept_tol) {
  mu_vals <- se2_mu <- sol_cond <- foc_vals <- rep(NA_real_, n_reps)
  method_used <- rep(NA_character_, n_reps)

  for (rep_i in 1:n_reps) {
    dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
    tryCatch({
      # Run constrained_nr
      single_ps <- list(
        formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
        h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
        inv_link = ps_spec[["inv_link"]],
        outcome = ps_spec[["outcome"]],
        alpha_init.list = list(NULL),
        optimizer = "constrained_nr"
      )
      ebmr_cnr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
      gmm_cnr <- ebmr_cnr$ps_fit.list[[1]]$gmm_fit
      cnr_grad <- gmm_cnr$opt$final_grad_norm
      cnr_cond <- gmm_cnr$opt$solution_cond

      # Decision: if CNR converged (grad < 1e-6) OR
      # if CNR is near-converged (grad < accept_tol) with non-degenerate cond,
      # use the CNR result. Otherwise use L-BFGS-B fallback result (which the
      # package already computed internally).
      # But we can't separate CNR-only vs fallback from the package output.
      # So: run L-BFGS-B separately and pick externally.

      if (cnr_grad < accept_tol && cnr_cond < 1e8) {
        # Accept CNR boundary solution
        ebmr <- ebmr_cnr; gmm_fit <- gmm_cnr
        method_used[rep_i] <- "cnr_accepted"
      } else if (cnr_grad < 1e-6) {
        # CNR converged (possibly via fallback to degenerate solution)
        ebmr <- ebmr_cnr; gmm_fit <- gmm_cnr
        method_used[rep_i] <- "cnr_converged"
      } else {
        # Neither: run L-BFGS-B
        single_ps$optimizer <- "L-BFGS-B"
        ebmr_lb <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
        gmm_lb <- ebmr_lb$ps_fit.list[[1]]$gmm_fit
        ebmr <- ebmr_lb; gmm_fit <- gmm_lb
        method_used[rep_i] <- "lbfgsb"
      }

      r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
      pi_hat <- ebmr$ps_fit.list[[1]]$fitted.values
      mu_vals[rep_i] <- mean(r_vec * y_vec / pi_hat)
      sol_cond[rep_i] <- gmm_fit$opt$solution_cond

      # Compute FOC at the solution
      GtW <- crossprod(gmm_fit$Gamma.hat, gmm_fit$W.hat)
      foc_vals[rep_i] <- max(abs(GtW %*% gmm_fit$eta_s))

      psi2_mat <- gmm_fit$psi2
      design_mat <- ebmr$ps_fit.list[[1]]$design_matrix
      link_type <- ebmr$ps_fit.list[[1]]$link_type

      if (!is.null(psi2_mat) && is.matrix(psi2_mat) && !any(is.na(psi2_mat))) {
        if (!is.null(link_type) && link_type == "logistic_complement") {
          dot_pi <- -design_mat * (pi_hat * (1 - pi_hat))
        } else {
          dot_pi <- design_mat * (pi_hat * (1 - pi_hat))
        }
        ry_ps_inv2 <- as.vector(r_vec * y_vec * (pi_hat^(-2)))
        H_alpha <- colMeans(dot_pi * ry_ps_inv2)
        mu_iid <- as.vector(t(r_vec/pi_hat*y_vec) - t(H_alpha) %*% psi2_mat)
        se2_mu[rep_i] <- sqrt(var(mu_iid)/n_val)
      }
    }, error = function(e) NULL)
  }

  valid <- !is.na(mu_vals) & !is.na(se2_mu)

  cat(sprintf("\n  M%d accept_tol=%.0e:\n", model_idx, accept_tol))
  cat(sprintf("    All valid:  n=%3d, Bias=%.4f, ESD=%.4f, ESE2/ESD=%.3f\n",
      sum(valid), mean(mu_vals[valid]) - mu_true, sd(mu_vals[valid]),
      mean(se2_mu[valid]) / sd(mu_vals[valid])))

  for (m in sort(unique(method_used[valid]))) {
    mask <- valid & method_used == m
    if (sum(mask) > 0) {
      cat(sprintf("    %-20s n=%3d, median_cond=%.1e, median_foc=%.1e",
          paste0(m, ":"), sum(mask), median(sol_cond[mask]), median(foc_vals[mask])))
      if (sum(mask) > 5) {
        cat(sprintf(", ESE2/ESD=%.3f", mean(se2_mu[mask]) / sd(mu_vals[mask])))
      }
      cat("\n")
    }
  }
}

cat("=== Relaxed tolerance test (M3 only) ===\n")
for (tol_val in c(1e-6, 1e-4, 1e-3, 5e-3, 1e-2)) {
  run_with_relaxed_tol(3, tol_val)
}

cat("\n=== Verify M1 and M2 unaffected (accept_tol=1e-2) ===\n")
run_with_relaxed_tol(1, 1e-2)
run_with_relaxed_tol(2, 1e-2)
