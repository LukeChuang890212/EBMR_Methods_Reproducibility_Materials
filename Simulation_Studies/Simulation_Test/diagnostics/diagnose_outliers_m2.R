setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
  library(EBMRalgorithmFast4)
})

W_func  <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_nu_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1_u2=dat$u1*dat$u2)

ps_spec <- get_ps_spec("9-alt1")
ps_sub  <- list(formula.list = ps_spec[["formula.list"]][2],
                h_alpha.list = ps_spec[["h_alpha.list"]][2],
                inv_link     = ps_spec[["inv_link"]],
                outcome      = ps_spec[["outcome"]])
mu_true <- get_mu_true("setting3")
all_data <- readRDS("Simulation_Data/Setting3.B1_n2000_replicate1000.RDS")
n_val <- 2000
n_reps <- 500

mu_vals <- se_vals <- numeric(n_reps)
hess_pd <- logical(n_reps)
alpha_mat <- matrix(NA, n_reps, 4)
obj_vals  <- numeric(n_reps)

for (i in 1:n_reps) {
  set.seed(12345 + i)
  dat  <- all_data[((i-1)*n_val + 1):(i*n_val), ]
  ebmr <- EBMRAlgorithmFast4$new("y", ps_sub, dat, W_func)
  res  <- ebmr$EBMR_IPW(h_nu=h_nu_fn, type="HT", se.fit=TRUE)
  mu_vals[i]   <- res$mu_ipw
  se_vals[i]   <- res$se_ipw
  fit          <- ebmr$ps_fit.list[[1]]$gmm_fit
  hess_pd[i]   <- isTRUE(fit$opt$hessian_pd)
  alpha_mat[i,] <- fit$estimates
  obj_vals[i]  <- fit$opt$objective
}

# Detect outliers using same IQR logic as clean_sim_result
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

# Summary by hess_pd and outlier status
cat("         hess_pd=T  hess_pd=F\n")
cat(sprintf("Normal:    %4d       %4d\n", sum(!is_out & hess_pd), sum(!is_out & !hess_pd)))
cat(sprintf("Outlier:   %4d       %4d\n", sum( is_out & hess_pd), sum( is_out & !hess_pd)))

# ESE/ESD after excluding outliers
cat(sprintf("\nExcluding outliers: ESD=%.4f  ESE=%.4f  ESE/ESD=%.3f  (n=%d)\n",
    sd(mu_vals[norm_idx]), mean(se_vals[norm_idx]),
    mean(se_vals[norm_idx])/sd(mu_vals[norm_idx]), sum(norm_idx)))

# Show outlier details
cat("\n=== Outlier details ===\n")
cat(sprintf("%-5s %-7s %-8s %-8s %-8s %s\n",
    "Rep", "hpd", "mu", "se", "obj", "alpha"))
for (i in which(is_out)) {
  cat(sprintf("%-5d %-7s %8.3f %8.4f %8.5f [%s]\n",
      i, hess_pd[i],
      mu_vals[i], se_vals[i], obj_vals[i],
      paste(round(alpha_mat[i,], 3), collapse=", ")))
}

# Alpha magnitude analysis
alpha_norm <- sqrt(rowSums(alpha_mat^2))
cat(sprintf("\n=== Alpha norm: outliers vs normal ===\n"))
cat(sprintf("Normal:  mean=%.3f  max=%.3f\n",
    mean(alpha_norm[norm_idx]), max(alpha_norm[norm_idx])))
cat(sprintf("Outlier: mean=%.3f  max=%.3f\n",
    mean(alpha_norm[is_out]), max(alpha_norm[is_out])))

# Per-component alpha stats for outliers vs normal
cat("\nPer-component alpha (intercept, y, u1, z2):\n")
comp_names <- c("intercept", "y", "u1", "z2")
for (j in 1:4) {
  cat(sprintf("  %s: normal=[%.3f, %.3f]  outlier=[%.3f, %.3f]\n",
      comp_names[j],
      min(alpha_mat[norm_idx, j]), max(alpha_mat[norm_idx, j]),
      min(alpha_mat[is_out, j]),  max(alpha_mat[is_out, j])))
}
