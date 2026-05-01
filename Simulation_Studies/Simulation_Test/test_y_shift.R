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
cond_limit <- 1e8
m_idx <- 3

# ===== Test 1: Stubborn degenerate reps with y+1 shift =====
stubborn_reps <- c(35, 38, 143, 155, 156, 190)
cat("=== Stubborn reps: y shifted by +1 ===\n")
for (rep_i in stubborn_reps) {
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
  dat_shifted <- dat
  dat_shifted$y <- dat$y + 1  # Shift y from {0,1} to {1,2}

  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[m_idx]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[m_idx]]),
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL),
    optimizer = "constrained_nr"
  )
  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat_shifted, W_fn)
    gmm_fit <- ebmr$ps_fit.list[[1]]$gmm_fit
    sol_cond <- gmm_fit$opt$solution_cond
    grad <- gmm_fit$opt$final_grad_norm
    pi_hat <- ebmr$ps_fit.list[[1]]$fitted.values
    r_vec <- dat$r
    y_shifted <- dat_shifted$y
    mu_shifted <- mean(r_vec * y_shifted / pi_hat)
    mu <- mu_shifted - 1  # Subtract shift
    alpha <- gmm_fit$alpha_hat
    cat(sprintf("  Rep %3d: grad=%.2e, cond=%.2e, mu=%.4f, alpha[2]=%.3f, #pi>0.99=%d\n",
        rep_i, grad, sol_cond, mu, alpha[2], sum(pi_hat > 0.99)))
  }, error = function(e) {
    cat(sprintf("  Rep %3d: ERROR %s\n", rep_i, e$message))
  })
}

# ===== Test 2: All 200 reps with y+1 shift =====
cat("\n=== All 200 reps: M3 with y shifted by +1 ===\n")
mu_vals <- se_vals <- cond_vals <- grad_vals <- rep(NA_real_, n_reps)

for (rep_i in 1:n_reps) {
  if (rep_i %% 50 == 0) cat(sprintf("  rep %d...\n", rep_i))
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
  dat_shifted <- dat
  dat_shifted$y <- dat$y + 1

  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[m_idx]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[m_idx]]),
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL),
    optimizer = "constrained_nr"
  )
  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat_shifted, W_fn)
    gmm_fit <- ebmr$ps_fit.list[[1]]$gmm_fit
    cond_vals[rep_i] <- gmm_fit$opt$solution_cond
    grad_vals[rep_i] <- gmm_fit$opt$final_grad_norm

    pi_hat <- ebmr$ps_fit.list[[1]]$fitted.values
    design_mat <- ebmr$ps_fit.list[[1]]$design_matrix
    link_type <- ebmr$ps_fit.list[[1]]$link_type
    r_vec <- dat$r
    y_shifted <- dat_shifted$y

    # SE using psi from package (now SE2/CUE)
    psi_mat <- gmm_fit$psi
    if (!is.null(psi_mat) && is.matrix(psi_mat) && !any(is.na(psi_mat))) {
      if (!is.null(link_type) && link_type == "logistic_complement") {
        dot_pi <- -design_mat * (pi_hat * (1 - pi_hat))
      } else {
        dot_pi <- design_mat * (pi_hat * (1 - pi_hat))
      }
      # Note: use y_shifted for the SE computation since that's what the model sees
      ry_ps_inv2 <- as.vector(r_vec * y_shifted * (pi_hat^(-2)))
      H_alpha <- colMeans(dot_pi * ry_ps_inv2)
      mu_iid <- as.vector(t(r_vec/pi_hat*y_shifted) - t(H_alpha) %*% psi_mat)
      mu_vals[rep_i] <- mean(r_vec * y_shifted / pi_hat) - 1  # Subtract shift
      # SE is the same whether shifted or not (shift is constant, doesn't affect variance)
      se_vals[rep_i] <- sqrt(var(mu_iid) / n_val)
    }
  }, error = function(e) {
    cat(sprintf("  Rep %d ERROR: %s\n", rep_i, e$message))
  })
}

valid <- !is.na(mu_vals) & !is.na(se_vals)
n_degen <- sum(cond_vals[valid] >= cond_limit, na.rm = TRUE)
n_conv <- sum(grad_vals[valid] < 1e-6, na.rm = TRUE)
cat(sprintf("\n  Valid: %d/%d, Degenerate: %d, Converged: %d\n",
    sum(valid), n_reps, n_degen, n_conv))
cat(sprintf("  Bias     = %.4f\n", mean(mu_vals[valid]) - mu_true))
cat(sprintf("  ESD      = %.4f\n", sd(mu_vals[valid])))
cat(sprintf("  ESE      = %.4f\n", mean(se_vals[valid])))
cat(sprintf("  ESE/ESD  = %.3f\n", mean(se_vals[valid]) / sd(mu_vals[valid])))

ci_lo <- mu_vals[valid] - 1.96 * se_vals[valid]
ci_hi <- mu_vals[valid] + 1.96 * se_vals[valid]
cp <- mean(mu_true >= ci_lo & mu_true <= ci_hi)
cat(sprintf("  CP       = %.3f\n", cp))

well_cond <- valid & cond_vals < cond_limit
if (sum(well_cond) > 1 && sum(well_cond) < sum(valid)) {
  cat(sprintf("\n  [Cond<1e8] n=%d, Bias=%.4f, ESD=%.4f, ESE=%.4f, ESE/ESD=%.3f\n",
      sum(well_cond), mean(mu_vals[well_cond]) - mu_true, sd(mu_vals[well_cond]),
      mean(se_vals[well_cond]), mean(se_vals[well_cond]) / sd(mu_vals[well_cond])))
}

# Show degenerate reps if any
degen_idx <- which(valid & cond_vals >= cond_limit)
if (length(degen_idx) > 0) {
  cat(sprintf("\n  Degenerate reps: %s\n", paste(degen_idx, collapse=", ")))
}
