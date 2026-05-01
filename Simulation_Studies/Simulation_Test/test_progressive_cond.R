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

# Simulate progressive relaxation externally:
# For each rep, try constrained_nr with cond_threshold = 1e8, 1e10, 1e12, then L-BFGS-B
run_progressive <- function(model_idx) {
  cond_levels <- c(1e8, 1e10, 1e12)

  mu_vals <- se2_mu <- sol_cond <- rep(NA_real_, n_reps)
  converged_flag <- rep(FALSE, n_reps)
  method_used <- rep(NA_character_, n_reps)

  for (rep_i in 1:n_reps) {
    dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]

    found <- FALSE
    for (cl in cond_levels) {
      # We can't pass cond_threshold to the package directly,
      # but we CAN run constrained_nr and check solution_cond
      # The package always uses cond_threshold=1e8 internally.
      # So for cl > 1e8, we need to simulate by running L-BFGS-B
      # and checking if solution_cond < cl.
      # Actually, we need a different approach:
      # Run constrained_nr (uses 1e8 guard) -> if converged, done
      # If not, run L-BFGS-B and check if its solution_cond < cl
      # This is a simplification but tests the concept.
      next
    }

    # Simpler approach: run both optimizers, then pick based on cond preference
    tryCatch({
      # Run constrained_nr (package default cond=1e8)
      single_ps_cnr <- list(
        formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
        h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
        inv_link = ps_spec[["inv_link"]],
        outcome = ps_spec[["outcome"]],
        alpha_init.list = list(NULL),
        optimizer = "constrained_nr"
      )
      ebmr_cnr <- EBMRAlgorithmFast4$new("y", single_ps_cnr, dat, W_fn)
      gmm_cnr <- ebmr_cnr$ps_fit.list[[1]]$gmm_fit
      cnr_cond <- gmm_cnr$opt$solution_cond
      cnr_grad <- gmm_cnr$opt$final_grad_norm
      cnr_conv <- gmm_cnr$opt$converged

      # If constrained_nr converged with low cond, use it
      if (cnr_conv && cnr_cond < 1e8) {
        ebmr <- ebmr_cnr; gmm_fit <- gmm_cnr
        method_used[rep_i] <- "cnr_1e8"
      } else {
        # Run L-BFGS-B to get the unconstrained solution
        single_ps_lb <- single_ps_cnr
        single_ps_lb$optimizer <- "L-BFGS-B"
        ebmr_lb <- EBMRAlgorithmFast4$new("y", single_ps_lb, dat, W_fn)
        gmm_lb <- ebmr_lb$ps_fit.list[[1]]$gmm_fit
        lb_cond <- gmm_lb$opt$solution_cond
        lb_conv <- gmm_lb$opt$converged

        # Progressive: prefer lowest-cond converged solution
        if (cnr_conv && cnr_cond < 1e10) {
          # CNR converged at boundary ~1e8, prefer it over degenerate L-BFGS-B
          ebmr <- ebmr_cnr; gmm_fit <- gmm_cnr
          method_used[rep_i] <- "cnr_1e10"
        } else if (lb_conv && lb_cond < 1e10) {
          ebmr <- ebmr_lb; gmm_fit <- gmm_lb
          method_used[rep_i] <- "lb_1e10"
        } else if (cnr_conv && cnr_cond < 1e12) {
          ebmr <- ebmr_cnr; gmm_fit <- gmm_cnr
          method_used[rep_i] <- "cnr_1e12"
        } else if (lb_conv && lb_cond < 1e12) {
          ebmr <- ebmr_lb; gmm_fit <- gmm_lb
          method_used[rep_i] <- "lb_1e12"
        } else if (lb_conv) {
          # Final fallback: accept whatever L-BFGS-B found
          ebmr <- ebmr_lb; gmm_fit <- gmm_lb
          method_used[rep_i] <- "lb_unconstrained"
        } else if (cnr_conv) {
          ebmr <- ebmr_cnr; gmm_fit <- gmm_cnr
          method_used[rep_i] <- "cnr_unconstrained"
        } else {
          # Neither converged, pick best grad
          if (cnr_grad <= gmm_lb$opt$final_grad_norm) {
            ebmr <- ebmr_cnr; gmm_fit <- gmm_cnr
          } else {
            ebmr <- ebmr_lb; gmm_fit <- gmm_lb
          }
          method_used[rep_i] <- "best_effort"
        }
      }

      r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
      pi_hat <- ebmr$ps_fit.list[[1]]$fitted.values
      mu_vals[rep_i] <- mean(r_vec * y_vec / pi_hat)
      converged_flag[rep_i] <- gmm_fit$opt$converged
      sol_cond[rep_i] <- gmm_fit$opt$solution_cond

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

  cat(sprintf("\n  M%d progressive constrained_nr:\n", model_idx))
  cat(sprintf("    All valid:  n=%3d, Bias=%.4f, ESD=%.4f, ESE2/ESD=%.3f\n",
      sum(valid), mean(mu_vals[valid]) - mu_true, sd(mu_vals[valid]),
      mean(se2_mu[valid]) / sd(mu_vals[valid])))

  # Breakdown by method used
  for (m in unique(method_used[valid])) {
    mask <- valid & method_used == m
    if (sum(mask) > 0) {
      cat(sprintf("    %-22s n=%3d, median_cond=%.1e",
          paste0(m, ":"), sum(mask), median(sol_cond[mask])))
      if (sum(mask) > 5) {
        cat(sprintf(", Bias=%.4f, ESD=%.4f, ESE2/ESD=%.3f",
            mean(mu_vals[mask]) - mu_true, sd(mu_vals[mask]),
            mean(se2_mu[mask]) / sd(mu_vals[mask])))
      }
      cat("\n")
    }
  }

  # By cond bucket
  for (thresh in c(1e8, 1e10, 1e12)) {
    mask <- valid & sol_cond < thresh
    if (sum(mask) > 10) {
      cat(sprintf("    cond < %.0e:          n=%3d, Bias=%.4f, ESD=%.4f, ESE2/ESD=%.3f\n",
          thresh, sum(mask), mean(mu_vals[mask]) - mu_true, sd(mu_vals[mask]),
          mean(se2_mu[mask]) / sd(mu_vals[mask])))
    }
  }
}

cat("=== Progressive cond relaxation test ===\n")
for (m in 1:3) {
  run_progressive(m)
}
