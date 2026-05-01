setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

n_val <- 2000
n_reps <- 200
ps_spec_orig <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
mu_true <- get_mu_true("setting4")
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
cond_limit <- 1e8
m_idx <- 3

# u1, z1 are binary -> squares are redundant
# Only add: u2^2, z2^2, and cross interactions involving at least one continuous
h_alpha_enrich <- function(dat) {
  cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2,
        u2_sq = dat$u2^2, z2_sq = dat$z2^2,
        u1_u2 = dat$u1*dat$u2, u1_z2 = dat$u1*dat$z2,
        z1_u2 = dat$z1*dat$u2, z1_z2 = dat$z1*dat$z2,
        u2_z2 = dat$u2*dat$z2)
}

stubborn_reps <- c(35, 38, 143, 155, 156, 190)

# ===== Test 1: Stubborn reps =====
cat("=== Stubborn reps with enriched h_alpha (11 instruments) ===\n")
for (rep_i in stubborn_reps) {
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
  single_ps <- list(
    formula.list = list(ps_spec_orig[["formula.list"]][[m_idx]]),
    h_alpha.list = list(h_alpha_enrich),
    inv_link = ps_spec_orig[["inv_link"]],
    outcome = ps_spec_orig[["outcome"]],
    alpha_init.list = list(NULL),
    optimizer = "constrained_nr"
  )
  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
    gmm_fit <- ebmr$ps_fit.list[[1]]$gmm_fit
    sol_cond <- gmm_fit$opt$solution_cond
    grad <- gmm_fit$opt$final_grad_norm
    pi_hat <- ebmr$ps_fit.list[[1]]$fitted.values
    r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
    mu <- mean(r_vec * y_vec / pi_hat)
    alpha <- gmm_fit$alpha_hat
    cat(sprintf("  Rep %3d: grad=%.2e, cond=%.2e, mu=%.4f, alpha[2]=%.3f, #pi>0.99=%d\n",
        rep_i, grad, sol_cond, mu, alpha[2], sum(pi_hat > 0.99)))
  }, error = function(e) {
    cat(sprintf("  Rep %3d: ERROR %s\n", rep_i, e$message))
  })
}

# ===== Test 2: All 200 reps =====
cat("\n=== All 200 reps: M3 with enriched h_alpha (11 instruments) ===\n")
mu_vals <- se2_vals <- cond_vals <- grad_vals <- rep(NA_real_, n_reps)

for (rep_i in 1:n_reps) {
  if (rep_i %% 50 == 0) cat(sprintf("  rep %d...\n", rep_i))
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
  single_ps <- list(
    formula.list = list(ps_spec_orig[["formula.list"]][[m_idx]]),
    h_alpha.list = list(h_alpha_enrich),
    inv_link = ps_spec_orig[["inv_link"]],
    outcome = ps_spec_orig[["outcome"]],
    alpha_init.list = list(NULL),
    optimizer = "constrained_nr"
  )
  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
    gmm_fit <- ebmr$ps_fit.list[[1]]$gmm_fit
    sol_cond <- gmm_fit$opt$solution_cond
    cond_vals[rep_i] <- sol_cond
    grad_vals[rep_i] <- gmm_fit$opt$final_grad_norm

    pi_hat <- ebmr$ps_fit.list[[1]]$fitted.values
    design_mat <- ebmr$ps_fit.list[[1]]$design_matrix
    link_type <- ebmr$ps_fit.list[[1]]$link_type
    r_vec <- dat[["r"]]; y_vec <- dat[["y"]]

    psi2_mat <- gmm_fit$psi2
    if (!is.null(psi2_mat) && is.matrix(psi2_mat) && !any(is.na(psi2_mat))) {
      if (!is.null(link_type) && link_type == "logistic_complement") {
        dot_pi <- -design_mat * (pi_hat * (1 - pi_hat))
      } else {
        dot_pi <- design_mat * (pi_hat * (1 - pi_hat))
      }
      ry_ps_inv2 <- as.vector(r_vec * y_vec * (pi_hat^(-2)))
      H_alpha <- colMeans(dot_pi * ry_ps_inv2)
      mu_iid <- as.vector(t(r_vec/pi_hat*y_vec) - t(H_alpha) %*% psi2_mat)
      mu_vals[rep_i] <- mean(r_vec * y_vec / pi_hat)
      se2_vals[rep_i] <- sqrt(var(mu_iid) / n_val)
    }
  }, error = function(e) {
    cat(sprintf("  Rep %d ERROR: %s\n", rep_i, e$message))
  })
}

valid <- !is.na(mu_vals) & !is.na(se2_vals)
n_degen <- sum(cond_vals[valid] >= cond_limit, na.rm = TRUE)
n_conv <- sum(grad_vals[valid] < 1e-6, na.rm = TRUE)
cat(sprintf("\n  Valid: %d/%d, Degenerate: %d, Converged: %d\n",
    sum(valid), n_reps, n_degen, n_conv))
cat(sprintf("  Bias  = %.4f\n", mean(mu_vals[valid]) - mu_true))
cat(sprintf("  ESD   = %.4f\n", sd(mu_vals[valid])))
cat(sprintf("  ESE2  = %.4f\n", mean(se2_vals[valid])))
cat(sprintf("  ESE2/ESD = %.3f\n", mean(se2_vals[valid]) / sd(mu_vals[valid])))

well_cond <- valid & cond_vals < cond_limit
if (sum(well_cond) > 1 && sum(well_cond) < sum(valid)) {
  cat(sprintf("\n  [Cond<1e8] n=%d, Bias=%.4f, ESD=%.4f, ESE2=%.4f, ESE2/ESD=%.3f\n",
      sum(well_cond), mean(mu_vals[well_cond]) - mu_true, sd(mu_vals[well_cond]),
      mean(se2_vals[well_cond]), mean(se2_vals[well_cond]) / sd(mu_vals[well_cond])))
}

degen_idx <- which(valid & cond_vals >= cond_limit)
if (length(degen_idx) > 0) {
  cat(sprintf("\n  Degenerate reps: %s\n", paste(degen_idx, collapse=", ")))
}
