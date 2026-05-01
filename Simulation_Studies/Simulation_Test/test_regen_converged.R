setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

n_val <- 2000
n_target <- 1000
ps_spec <- get_ps_spec("9-alt1")
mu_true <- get_mu_true("setting4")
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
model_idx <- 3

# Storage
mu_vals <- se2_vals <- cond_vals <- numeric(n_target)
n_tried <- 0
n_done <- 0

set.seed(2026)

while (n_done < n_target) {
  n_tried <- n_tried + 1

  # Generate fresh data
  dat <- setting4.B1(n_val)

  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL),
    optimizer = "constrained_nr"
  )

  ok <- tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
    gmm_fit <- ebmr$ps_fit.list[[1]]$gmm_fit
    sol_cond <- gmm_fit$opt$solution_cond

    if (sol_cond >= 1e8) {
      FALSE  # degenerate, skip
    } else {
      pi_hat <- ebmr$ps_fit.list[[1]]$fitted.values
      design_mat <- ebmr$ps_fit.list[[1]]$design_matrix
      link_type <- ebmr$ps_fit.list[[1]]$link_type
      r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
      mu <- mean(r_vec * y_vec / pi_hat)

      psi2_mat <- gmm_fit$psi2
      se2 <- NA
      if (!is.null(psi2_mat) && is.matrix(psi2_mat) && !any(is.na(psi2_mat))) {
        if (!is.null(link_type) && link_type == "logistic_complement") {
          dot_pi <- -design_mat * (pi_hat * (1 - pi_hat))
        } else {
          dot_pi <- design_mat * (pi_hat * (1 - pi_hat))
        }
        ry_ps_inv2 <- as.vector(r_vec * y_vec * (pi_hat^(-2)))
        H_alpha <- colMeans(dot_pi * ry_ps_inv2)
        mu_iid <- as.vector(t(r_vec/pi_hat*y_vec) - t(H_alpha) %*% psi2_mat)
        se2 <- sqrt(var(mu_iid) / n_val)
      }

      if (!is.na(se2)) {
        n_done <- n_done + 1
        mu_vals[n_done] <- mu
        se2_vals[n_done] <- se2
        cond_vals[n_done] <- sol_cond
        TRUE
      } else {
        FALSE
      }
    }
  }, error = function(e) FALSE)

  if (n_tried %% 50 == 0) {
    cat(sprintf("  tried=%d, done=%d\n", n_tried, n_done))
  }
}

cat(sprintf("\n=== Results: 200 converged reps (regenerated) ===\n"))
cat(sprintf("  Total tried: %d (skipped %d = %.1f%%)\n",
    n_tried, n_tried - n_target, 100*(n_tried - n_target)/n_tried))
cat(sprintf("  mu_true = %.6f\n", mu_true))
cat(sprintf("  Bias    = %.4f\n", mean(mu_vals) - mu_true))
cat(sprintf("  ESD     = %.4f\n", sd(mu_vals)))
cat(sprintf("  ESE2    = %.4f\n", mean(se2_vals)))
cat(sprintf("  ESE2/ESD = %.3f\n", mean(se2_vals) / sd(mu_vals)))
cat(sprintf("  CP95    = %.3f\n",
    mean(abs(mu_vals - mu_true) < 1.96 * se2_vals)))
