setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
  devtools::load_all("../EBMRalgorithmFast4")
})
W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
n_val <- 2000
n_reps <- 200

ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)

# Check model 3 alpha convergence across reps
alpha_mat <- matrix(NA, n_reps, 4)
obj_vals <- rep(NA, n_reps)
grad_vals <- rep(NA, n_reps)
conv_vals <- rep(NA, n_reps)
iter_vals <- rep(NA, n_reps)

ps_sub3 <- list(
  formula.list = ps_spec[["formula.list"]][3],
  h_alpha.list = ps_spec[["h_alpha.list"]][3],
  inv_link     = ps_spec[["inv_link"]],
  outcome      = ps_spec[["outcome"]]
)

for (i in 1:n_reps) {
  dat <- all_data[((i-1)*n_val + 1):(i*n_val), ]
  tryCatch({
    ebmr <- EBMRAlgorithmFast4[["new"]]("y", ps_sub3, dat, W_func)
    gf <- ebmr[["ps_fit.list"]][[1]][["gmm_fit"]]
    alpha_mat[i,] <- gf[["estimates"]]
    obj_vals[i] <- gf[["opt"]][["objective"]]
    grad_vals[i] <- gf[["opt"]][["final_grad_norm"]]
    conv_vals[i] <- gf[["opt"]][["converged"]]
    iter_vals[i] <- gf[["opt"]][["iterations"]]
  }, error = function(e) {
    cat(sprintf("Error rep %d: %s\n", i, conditionMessage(e)))
  })
}

cat("=== Model 3 alpha convergence diagnostics ===\n\n")
cat(sprintf("Converged: %d / %d\n", sum(conv_vals, na.rm=TRUE), n_reps))
cat(sprintf("NOT converged: %d\n", sum(!conv_vals, na.rm=TRUE)))

cat(sprintf("\nObjective: mean=%.6f, median=%.6f, max=%.6f\n",
    mean(obj_vals, na.rm=TRUE), median(obj_vals, na.rm=TRUE), max(obj_vals, na.rm=TRUE)))
cat(sprintf("Grad norm: mean=%.2e, median=%.2e, max=%.2e\n",
    mean(grad_vals, na.rm=TRUE), median(grad_vals, na.rm=TRUE), max(grad_vals, na.rm=TRUE)))
cat(sprintf("Iterations: mean=%.1f, median=%.0f, max=%d\n",
    mean(iter_vals, na.rm=TRUE), median(iter_vals, na.rm=TRUE), max(iter_vals, na.rm=TRUE)))

cat("\nAlpha estimates (intercept, y, u2, z2):\n")
for (j in 1:4) {
  cat(sprintf("  alpha[%d]: mean=%.4f, sd=%.4f, range=[%.4f, %.4f]\n",
      j, mean(alpha_mat[,j], na.rm=TRUE), sd(alpha_mat[,j], na.rm=TRUE),
      min(alpha_mat[,j], na.rm=TRUE), max(alpha_mat[,j], na.rm=TRUE)))
}

# Check for outlier alphas
cat("\n--- Reps with large grad norm (>1e-4) ---\n")
bad <- which(grad_vals > 1e-4)
if (length(bad) > 0) {
  for (idx in bad) {
    cat(sprintf("  Rep %d: grad=%.2e, obj=%.6f, iter=%d, alpha=%s\n",
        idx, grad_vals[idx], obj_vals[idx], iter_vals[idx],
        paste(round(alpha_mat[idx,], 4), collapse=", ")))
  }
} else {
  cat("  None\n")
}

cat("\n--- Reps with large objective (>0.01) ---\n")
bad2 <- which(obj_vals > 0.01)
if (length(bad2) > 0) {
  for (idx in head(bad2, 20)) {
    cat(sprintf("  Rep %d: obj=%.6f, grad=%.2e, iter=%d, alpha=%s\n",
        idx, obj_vals[idx], grad_vals[idx], iter_vals[idx],
        paste(round(alpha_mat[idx,], 4), collapse=", ")))
  }
} else {
  cat("  None\n")
}

cat("\n--- Reps that hit max iterations (200) ---\n")
bad3 <- which(iter_vals >= 200)
if (length(bad3) > 0) {
  for (idx in head(bad3, 20)) {
    cat(sprintf("  Rep %d: iter=%d, grad=%.2e, obj=%.6f\n",
        idx, iter_vals[idx], grad_vals[idx], obj_vals[idx]))
  }
} else {
  cat("  None\n")
}

# Also compare with GLM init
cat("\n\n=== Compare GMM alpha3 vs GLM init (first 10 reps) ===\n")
for (i in 1:10) {
  dat <- all_data[((i-1)*n_val + 1):(i*n_val), ]
  glm_fit <- glm(r ~ y + u2 + z2, data=dat, family=binomial(link="logit"))
  cat(sprintf("Rep %d: GLM=(%s)  GMM=(%s)\n", i,
      paste(round(coef(glm_fit), 4), collapse=", "),
      paste(round(alpha_mat[i,], 4), collapse=", ")))
}
