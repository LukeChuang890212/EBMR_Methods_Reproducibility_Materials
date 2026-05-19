setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
  library(EBMRalgorithmFast4)
})

W_func  <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_nu_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1_u2=dat$u1*dat$u2)

# Model 3 only from scenario 9-3 (ps_spec = 9-alt1): formula = u2_z2, h_alpha = full
ps_spec <- get_ps_spec("9-alt1")
ps_sub  <- list(formula.list = ps_spec[["formula.list"]][3],
                h_alpha.list = ps_spec[["h_alpha.list"]][3],
                inv_link     = ps_spec[["inv_link"]],
                outcome      = ps_spec[["outcome"]])
mu_true  <- 0.75
all_data <- readRDS("Simulation_Data/setting4.B1_n2000_replicate1000.RDS")
n_val    <- 2000
n_reps   <- 500

mu_vals  <- se_vals <- numeric(n_reps)
hess_pd  <- logical(n_reps)
obj_vals <- numeric(n_reps)
alpha_mat <- matrix(NA, n_reps, 4)  # model 3 has 4 params: intercept, u2, z2 + h_alpha full=4

for (i in 1:n_reps) {
  set.seed(12345 + i)
  dat  <- all_data[((i-1)*n_val + 1):(i*n_val), ]
  ebmr <- EBMRAlgorithmFast4$new("y", ps_sub, dat, W_func)
  res  <- ebmr$EBMR_IPW(h_nu=h_nu_fn, type="HT", se.fit=TRUE)
  mu_vals[i]  <- res$mu_ipw
  se_vals[i]  <- res$se_ipw
  fit          <- ebmr$ps_fit.list[[1]]$gmm_fit
  hess_pd[i]  <- isTRUE(fit$opt$hessian_pd)
  obj_vals[i] <- fit$opt$objective
  alpha_mat[i, seq_along(fit$estimates)] <- fit$estimates
}

# Outlier detection via clean_sim_result
sim_result <- rbind(mu_ipw      = mu_vals,
                    mu_ipw.true = rep(mu_true, n_reps),
                    se_ipw      = se_vals,
                    se_ipw.true = se_vals)
colnames(sim_result) <- as.character(seq_len(n_reps))
cleaned  <- clean_sim_result(sim_result, multiplier = 3, max_pct = 0.01, verbose = FALSE)
kept_idx <- as.integer(colnames(cleaned$result))
is_out   <- !(seq_len(n_reps) %in% kept_idx)
norm_idx <- !is_out

cat(sprintf("n_reps=%d  mu_true=%.4f  ESD=%.4f  ESE=%.4f  ESE/ESD=%.3f\n",
    n_reps, mu_true, sd(mu_vals), mean(se_vals), mean(se_vals)/sd(mu_vals)))
cat(sprintf("Outliers (IQR method, cap 1%%): %d/%d\n\n", sum(is_out), n_reps))
cat(sprintf("Excluding outliers: ESD=%.4f  ESE=%.4f  ESE/ESD=%.3f  (n=%d)\n\n",
    sd(mu_vals[norm_idx]), mean(se_vals[norm_idx]),
    mean(se_vals[norm_idx])/sd(mu_vals[norm_idx]), sum(norm_idx)))

# GMM objective distribution
cat("=== GMM objective distribution ===\n")
cat(sprintf("mean=%.6f  sd=%.6f  min=%.6f  max=%.6f\n",
    mean(obj_vals), sd(obj_vals), min(obj_vals), max(obj_vals)))
cat(sprintf("obj < 1e-6:  %d/%d (%.1f%%)\n",
    sum(obj_vals < 1e-6), n_reps, 100*mean(obj_vals < 1e-6)))
cat(sprintf("obj < 1e-4:  %d/%d (%.1f%%)\n\n",
    sum(obj_vals < 1e-4), n_reps, 100*mean(obj_vals < 1e-4)))

# Alpha distribution
alpha_dim <- sum(!is.na(alpha_mat[1,]))
alpha_norm <- sqrt(rowSums(alpha_mat[, 1:alpha_dim, drop=FALSE]^2, na.rm=TRUE))
cat("=== Alpha norm ===\n")
cat(sprintf("mean=%.3f  sd=%.3f  min=%.3f  max=%.3f\n\n",
    mean(alpha_norm), sd(alpha_norm), min(alpha_norm), max(alpha_norm)))

# Check for near-zero obj with extreme alpha
extreme_alpha <- alpha_norm > 10
near_zero_obj <- obj_vals < 1e-4
cat(sprintf("Near-zero obj (< 1e-4) AND extreme alpha (norm > 10): %d/%d\n",
    sum(near_zero_obj & extreme_alpha), n_reps))
cat(sprintf("Near-zero obj only: %d/%d\n",
    sum(near_zero_obj & !extreme_alpha), n_reps))
cat(sprintf("Extreme alpha only: %d/%d\n\n",
    sum(!near_zero_obj & extreme_alpha), n_reps))

# Per-component alpha stats
cat("=== Alpha per-component (all reps) ===\n")
for (j in 1:alpha_dim) {
  cat(sprintf("  alpha[%d]: mean=%7.3f  sd=%6.3f  min=%8.3f  max=%8.3f\n",
      j, mean(alpha_mat[,j], na.rm=TRUE), sd(alpha_mat[,j], na.rm=TRUE),
      min(alpha_mat[,j], na.rm=TRUE), max(alpha_mat[,j], na.rm=TRUE)))
}

# Show extreme-alpha reps
if (sum(extreme_alpha) > 0) {
  cat("\n=== Reps with alpha norm > 10 ===\n")
  cat(sprintf("%-5s %-8s %-8s %-8s %s\n", "Rep", "mu", "obj", "alpha_norm", "alpha"))
  for (i in which(extreme_alpha)) {
    cat(sprintf("%-5d %8.4f %8.5f %10.3f  [%s]\n",
        i, mu_vals[i], obj_vals[i], alpha_norm[i],
        paste(round(alpha_mat[i, 1:alpha_dim], 3), collapse=", ")))
  }
}
