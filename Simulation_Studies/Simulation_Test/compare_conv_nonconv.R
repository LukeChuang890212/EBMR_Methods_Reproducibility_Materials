setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
  devtools::load_all("../EBMRalgorithmFast4")
})

W_func  <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_nu_fn <- function(dat) cbind(u1=dat[["u1"]], u2=dat[["u2"]], z1=dat[["z1"]], z2=dat[["z2"]], u1_u2=dat[["u1"]]*dat[["u2"]])
n_val   <- 2000
n_reps  <- 200

ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)

ps_sub3 <- list(
  formula.list = ps_spec[["formula.list"]][3],
  h_alpha.list = ps_spec[["h_alpha.list"]][3],
  inv_link     = ps_spec[["inv_link"]],
  outcome      = ps_spec[["outcome"]]
)

alpha_mat  <- matrix(NA, n_reps, 4)
conv_vals  <- rep(NA, n_reps)
grad_vals  <- rep(NA_real_, n_reps)
obj_vals   <- rep(NA_real_, n_reps)
iter_vals  <- rep(NA_integer_, n_reps)
mu_vals    <- rep(NA_real_, n_reps)
pi_mean_vals <- rep(NA_real_, n_reps)

for (i in 1:n_reps) {
  dat <- all_data[((i-1)*n_val + 1):(i*n_val), ]
  tryCatch({
    ebmr <- EBMRAlgorithmFast4[["new"]]("y", ps_sub3, dat, W_func)
    res  <- ebmr[["EBMR_IPW"]](h_nu=h_nu_fn, type="HT", se.fit=FALSE)
    gf <- ebmr[["ps_fit.list"]][[1]][["gmm_fit"]]
    alpha_mat[i,]  <- gf[["estimates"]]
    conv_vals[i]   <- gf[["opt"]][["converged"]]
    grad_vals[i]   <- gf[["opt"]][["final_grad_norm"]]
    obj_vals[i]    <- gf[["opt"]][["objective"]]
    iter_vals[i]   <- gf[["opt"]][["iterations"]]
    mu_vals[i]     <- res[["mu_ipw"]]
    pi_mean_vals[i] <- mean(ebmr[["ps_fit.list"]][[1]][["fitted.values"]])
  }, error = function(e) NULL)
}

valid <- !is.na(mu_vals)
conv <- conv_vals[valid]
colnames(alpha_mat) <- c("intercept", "y", "u2", "z2")

cat("=== Alpha estimates: Converged vs Not converged ===\n\n")
for (j in 1:4) {
  a_c  <- alpha_mat[valid, j][conv]
  a_nc <- alpha_mat[valid, j][!conv]
  cat(sprintf("alpha[%d] (%s):\n", j, colnames(alpha_mat)[j]))
  cat(sprintf("  Converged (n=%d):     mean=%7.4f, sd=%7.4f, median=%7.4f, range=[%7.4f, %7.4f]\n",
      sum(conv), mean(a_c), sd(a_c), median(a_c), min(a_c), max(a_c)))
  cat(sprintf("  Not converged (n=%d): mean=%7.4f, sd=%7.4f, median=%7.4f, range=[%7.4f, %7.4f]\n\n",
      sum(!conv), mean(a_nc), sd(a_nc), median(a_nc), min(a_nc), max(a_nc)))
}

cat("=== Objective ===\n")
cat(sprintf("  Converged:     mean=%.6f, sd=%.6f, median=%.6f, range=[%.6f, %.6f]\n",
    mean(obj_vals[valid][conv]), sd(obj_vals[valid][conv]), median(obj_vals[valid][conv]),
    min(obj_vals[valid][conv]), max(obj_vals[valid][conv])))
cat(sprintf("  Not converged: mean=%.6f, sd=%.6f, median=%.6f, range=[%.6f, %.6f]\n\n",
    mean(obj_vals[valid][!conv]), sd(obj_vals[valid][!conv]), median(obj_vals[valid][!conv]),
    min(obj_vals[valid][!conv]), max(obj_vals[valid][!conv])))

cat("=== Gradient norm ===\n")
cat(sprintf("  Converged:     mean=%.2e, sd=%.2e, median=%.2e, range=[%.2e, %.2e]\n",
    mean(grad_vals[valid][conv]), sd(grad_vals[valid][conv]), median(grad_vals[valid][conv]),
    min(grad_vals[valid][conv]), max(grad_vals[valid][conv])))
cat(sprintf("  Not converged: mean=%.2e, sd=%.2e, median=%.2e, range=[%.2e, %.2e]\n\n",
    mean(grad_vals[valid][!conv]), sd(grad_vals[valid][!conv]), median(grad_vals[valid][!conv]),
    min(grad_vals[valid][!conv]), max(grad_vals[valid][!conv])))

cat("=== Iterations ===\n")
cat(sprintf("  Converged:     mean=%.1f, median=%.0f, range=[%d, %d]\n",
    mean(iter_vals[valid][conv]), median(iter_vals[valid][conv]),
    min(iter_vals[valid][conv]), max(iter_vals[valid][conv])))
cat(sprintf("  Not converged: mean=%.1f, median=%.0f, range=[%d, %d]\n\n",
    mean(iter_vals[valid][!conv]), median(iter_vals[valid][!conv]),
    min(iter_vals[valid][!conv]), max(iter_vals[valid][!conv])))

cat("=== Mean pi(alpha) ===\n")
cat(sprintf("  Converged:     mean=%.4f, sd=%.4f, range=[%.4f, %.4f]\n",
    mean(pi_mean_vals[valid][conv]), sd(pi_mean_vals[valid][conv]),
    min(pi_mean_vals[valid][conv]), max(pi_mean_vals[valid][conv])))
cat(sprintf("  Not converged: mean=%.4f, sd=%.4f, range=[%.4f, %.4f]\n\n",
    mean(pi_mean_vals[valid][!conv]), sd(pi_mean_vals[valid][!conv]),
    min(pi_mean_vals[valid][!conv]), max(pi_mean_vals[valid][!conv])))

cat("=== mu_ipw ===\n")
cat(sprintf("  Converged:     mean=%.4f, sd=%.4f, range=[%.4f, %.4f]\n",
    mean(mu_vals[valid][conv]), sd(mu_vals[valid][conv]),
    min(mu_vals[valid][conv]), max(mu_vals[valid][conv])))
cat(sprintf("  Not converged: mean=%.4f, sd=%.4f, range=[%.4f, %.4f]\n\n",
    mean(mu_vals[valid][!conv]), sd(mu_vals[valid][!conv]),
    min(mu_vals[valid][!conv]), max(mu_vals[valid][!conv])))

# Print first 10 converged and first 10 not-converged
cat("=== Sample of converged reps ===\n")
c_idx <- which(valid)[conv]
for (ii in head(c_idx, 10)) {
  cat(sprintf("  Rep %3d: alpha=(%7.3f,%7.3f,%7.3f,%7.3f), obj=%.6f, grad=%.2e, iter=%3d, mu=%.4f, mean_pi=%.4f\n",
      ii, alpha_mat[ii,1], alpha_mat[ii,2], alpha_mat[ii,3], alpha_mat[ii,4],
      obj_vals[ii], grad_vals[ii], iter_vals[ii], mu_vals[ii], pi_mean_vals[ii]))
}

cat("\n=== Sample of NOT converged reps ===\n")
nc_idx <- which(valid)[!conv]
for (ii in head(nc_idx, 10)) {
  cat(sprintf("  Rep %3d: alpha=(%7.3f,%7.3f,%7.3f,%7.3f), obj=%.6f, grad=%.2e, iter=%3d, mu=%.4f, mean_pi=%.4f\n",
      ii, alpha_mat[ii,1], alpha_mat[ii,2], alpha_mat[ii,3], alpha_mat[ii,4],
      obj_vals[ii], grad_vals[ii], iter_vals[ii], mu_vals[ii], pi_mean_vals[ii]))
}
