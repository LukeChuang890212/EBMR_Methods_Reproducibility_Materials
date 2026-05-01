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
h_nu_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1_u2=dat$u1*dat$u2)

# Use models 2+3
model_set <- c(2, 3)

mu_vals <- se_vals <- rep(NA_real_, n_reps)
cond_m2 <- cond_m3 <- rep(NA_real_, n_reps)

for (rep_i in 1:n_reps) {
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]

  ps_sub <- list(
    formula.list = ps_spec[["formula.list"]][model_set],
    h_alpha.list = ps_spec[["h_alpha.list"]][model_set],
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL, NULL),
    optimizer = "constrained_nr"
  )

  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", ps_sub, dat, W_fn)
    res <- ebmr$EBMR_IPW(h_nu = h_nu_fn, type = "HT", se.fit = TRUE)
    mu_vals[rep_i] <- res$mu_ipw
    se_vals[rep_i] <- res$se_ipw
    cond_m2[rep_i] <- ebmr$ps_fit.list[[1]]$gmm_fit$opt$solution_cond
    cond_m3[rep_i] <- ebmr$ps_fit.list[[2]]$gmm_fit$opt$solution_cond
  }, error = function(e) {
    cat(sprintf("  Rep %d ERROR: %s\n", rep_i, e$message))
  })
}

valid <- !is.na(mu_vals) & !is.na(se_vals)
cat(sprintf("\n=== Ensemble M2+M3, setting4, n=2000 ===\n"))
cat(sprintf("  Valid: %d/%d\n", sum(valid), n_reps))
cat(sprintf("  Bias  = %.4f\n", mean(mu_vals[valid]) - mu_true))
cat(sprintf("  ESD   = %.4f\n", sd(mu_vals[valid])))
cat(sprintf("  ESE   = %.4f\n", mean(se_vals[valid])))
cat(sprintf("  ESE/ESD = %.3f\n", mean(se_vals[valid]) / sd(mu_vals[valid])))

# Check degeneracy per model
cat(sprintf("\n  M2 degenerate (cond>=1e8): %d\n", sum(cond_m2[valid] >= 1e8, na.rm = TRUE)))
cat(sprintf("  M3 degenerate (cond>=1e8): %d\n", sum(cond_m3[valid] >= 1e8, na.rm = TRUE)))
cat(sprintf("  Either degenerate: %d\n",
    sum(cond_m2[valid] >= 1e8 | cond_m3[valid] >= 1e8, na.rm = TRUE)))

# Report conditional on both converged
both_ok <- valid & cond_m2 < 1e8 & cond_m3 < 1e8
cat(sprintf("\n=== Both M2 & M3 converged (cond<1e8) ===\n"))
cat(sprintf("  n=%d\n", sum(both_ok)))
if (sum(both_ok) > 10) {
  cat(sprintf("  Bias  = %.4f\n", mean(mu_vals[both_ok]) - mu_true))
  cat(sprintf("  ESD   = %.4f\n", sd(mu_vals[both_ok])))
  cat(sprintf("  ESE   = %.4f\n", mean(se_vals[both_ok])))
  cat(sprintf("  ESE/ESD = %.3f\n", mean(se_vals[both_ok]) / sd(mu_vals[both_ok])))
}

# Also run individual models for comparison
cat("\n=== Individual models for comparison ===\n")
for (m_idx in c(2, 3)) {
  mu_m <- se_m <- rep(NA_real_, n_reps)
  for (rep_i in 1:n_reps) {
    dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
    ps_sub <- list(
      formula.list = list(ps_spec[["formula.list"]][[m_idx]]),
      h_alpha.list = list(ps_spec[["h_alpha.list"]][[m_idx]]),
      inv_link = ps_spec[["inv_link"]],
      outcome = ps_spec[["outcome"]],
      alpha_init.list = list(NULL),
      optimizer = "constrained_nr"
    )
    tryCatch({
      ebmr <- EBMRAlgorithmFast4$new("y", ps_sub, dat, W_fn)
      res <- ebmr$EBMR_IPW(h_nu = h_nu_fn, type = "HT", se.fit = TRUE)
      mu_m[rep_i] <- res$mu_ipw
      se_m[rep_i] <- res$se_ipw
    }, error = function(e) NULL)
  }
  v <- !is.na(mu_m) & !is.na(se_m)
  cat(sprintf("  M%d: Valid=%d, Bias=%.4f, ESD=%.4f, ESE=%.4f, ESE/ESD=%.3f\n",
      m_idx, sum(v), mean(mu_m[v]) - mu_true, sd(mu_m[v]),
      mean(se_m[v]), mean(se_m[v]) / sd(mu_m[v])))
}
