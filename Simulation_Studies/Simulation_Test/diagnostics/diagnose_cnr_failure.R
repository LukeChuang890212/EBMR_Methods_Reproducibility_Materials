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

# Test a single rep for M2 and M3 with constrained_nr, with verbose diagnostics
diagnose_single <- function(model_idx, rep_i, cond_thresh = 1e8, trust_r = 2.0) {
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
  n <- n_val
  r_vec <- dat[["r"]]

  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL),
    optimizer = "constrained_nr"
  )

  # Run with default settings
  ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
  gmm_fit <- ebmr$ps_fit.list[[1]]$gmm_fit
  cat(sprintf("  Default: converged=%s, grad=%.2e, obj=%.6f, iter=%d\n",
      gmm_fit$opt$converged, gmm_fit$opt$final_grad_norm,
      gmm_fit$opt$objective, gmm_fit$opt$iterations))
  cat(sprintf("  alpha = %s\n", paste(round(gmm_fit$estimates, 4), collapse=", ")))

  # Check condition number at solution
  Gamma <- gmm_fit$Gamma.hat
  W_hat <- gmm_fit$W.hat
  H <- crossprod(Gamma, W_hat %*% Gamma)
  eigs <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  cond <- max(eigs) / max(min(eigs), 1e-15)
  cat(sprintf("  cond(GtWG) at solution = %.2e\n", cond))

  # Check pi range
  pi_hat <- ebmr$ps_fit.list[[1]]$fitted.values
  cat(sprintf("  pi range: [%.6f, %.6f], mean=%.4f\n", min(pi_hat), max(pi_hat), mean(pi_hat)))
  cat(sprintf("  r mean=%.4f, n_respond=%d\n", mean(r_vec), sum(r_vec)))

  return(list(converged = gmm_fit$opt$converged, cond = cond))
}

# M2: diagnose failure
cat("=== M2 constrained_nr (almost completely fails) ===\n")
for (rep_i in 1:5) {
  cat(sprintf("\nRep %d:\n", rep_i))
  diagnose_single(2, rep_i)
}

# M3: find a failing rep
cat("\n\n=== M3 constrained_nr (22/200 fail) ===\n")
# Run through reps to find non-converged ones
fail_reps <- c()
for (rep_i in 1:200) {
  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[3]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[3]]),
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL),
    optimizer = "constrained_nr"
  )
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
    if (!ebmr$ps_fit.list[[1]]$gmm_fit$opt$converged) {
      fail_reps <- c(fail_reps, rep_i)
    }
  }, error = function(e) { fail_reps <<- c(fail_reps, rep_i) })
  if (length(fail_reps) >= 5) break
}

cat(sprintf("First failed reps: %s\n", paste(fail_reps, collapse=", ")))
for (rep_i in fail_reps[1:min(3, length(fail_reps))]) {
  cat(sprintf("\nFailed Rep %d:\n", rep_i))
  diagnose_single(3, rep_i)
}

# Also show a successful rep for comparison
cat("\n\nSuccessful Rep 1:\n")
diagnose_single(3, 1)

# ============================================================
# Now test with relaxed cond_threshold
# ============================================================
cat("\n\n=== Testing relaxed cond_threshold for M2 ===\n")

# We need to call gmm directly with different cond_threshold
# Rebuild the model with different thresholds
for (ct in c(1e8, 1e10, 1e12, Inf)) {
  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[2]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[2]]),
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL),
    optimizer = "constrained_nr"
  )
  dat <- all_data[1:n_val, ]

  # To pass cond_threshold, we need to modify EBMRAlgorithm to forward it
  # For now, directly call the internal gmm function
  # Let's temporarily patch the package

  # Actually, the gmm function accepts cond_threshold as parameter
  # But EBMRAlgorithm.r might not pass it through
  # Let me check what happens with L-BFGS-B first as baseline
  single_ps_lb <- single_ps
  single_ps_lb$optimizer <- "L-BFGS-B"
  ebmr_lb <- EBMRAlgorithmFast4$new("y", single_ps_lb, dat, W_fn)
  gmm_lb <- ebmr_lb$ps_fit.list[[1]]$gmm_fit

  cat(sprintf("\n  cond_threshold=%s:\n", format(ct, scientific=TRUE)))
  cat(sprintf("    L-BFGS-B baseline: conv=%s, grad=%.2e, obj=%.6f\n",
      gmm_lb$opt$converged, gmm_lb$opt$final_grad_norm, gmm_lb$opt$objective))
  cat(sprintf("    L-BFGS-B alpha = %s\n", paste(round(gmm_lb$estimates, 4), collapse=", ")))

  # Check cond at L-BFGS-B solution
  Gamma_lb <- gmm_lb$Gamma.hat
  W_lb <- gmm_lb$W.hat
  H_lb <- crossprod(Gamma_lb, W_lb %*% Gamma_lb)
  eigs_lb <- eigen(H_lb, symmetric = TRUE, only.values = TRUE)$values
  cond_lb <- max(eigs_lb) / max(min(eigs_lb), 1e-15)
  cat(sprintf("    L-BFGS-B cond(GtWG) = %.2e\n", cond_lb))
}
