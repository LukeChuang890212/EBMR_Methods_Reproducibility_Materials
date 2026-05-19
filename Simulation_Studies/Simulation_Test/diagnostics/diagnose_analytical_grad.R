setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
source("Data_Generation.r")
source("config/scenarios.R")

old_wd <- getwd()
setwd("../EBMRalgorithmFast4")
source("R/EBMRAlgorithm.r")
setwd(old_wd)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
inv_link_fn <- function(eta) 1 / (1 + exp(eta))

ps_spec_m2 <- list(
  formula.list = list(FORMULAS$u1_z2),
  h_alpha.list = list(H_ALPHA$full),
  inv_link = inv_link_fn,
  outcome = "y",
  alpha_init.list = list(NULL)
)

set.seed(999)
big <- setting3.A1(1e6)
mu_true <- mean(big$y)
rm(big); gc()
cat(sprintf("mu_true = %.6f\n\n", mu_true))

set.seed(42)
n_reps <- 200
n_val <- 2000

mu_hat <- numeric(n_reps)
alpha_mat <- list()
converged <- logical(n_reps)
iterations <- integer(n_reps)
objective <- numeric(n_reps)
min_ps <- numeric(n_reps)
max_ipw <- numeric(n_reps)

for (i in 1:n_reps) {
  dat <- setting3.B1(n_val)

  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", ps_spec_m2, dat, W_func, bounds = 5)

    fit <- ebmr$ps_fit.list[[1]]
    alpha_hat <- fit$coefficients
    ps_vals <- fit$fitted.values
    gmm_fit <- fit$gmm_fit

    mu_hat[i] <- mean(dat$r / ps_vals * dat$y)
    alpha_mat[[i]] <- alpha_hat
    converged[i] <- gmm_fit$opt$converged
    iterations[i] <- gmm_fit$opt$iterations
    objective[i] <- gmm_fit$opt$objective
    min_ps[i] <- min(ps_vals)
    max_ipw[i] <- max(dat$r / ps_vals)

  }, error = function(e) {
    mu_hat[i] <<- NA
    alpha_mat[[i]] <<- rep(NA, 4)
    converged[i] <<- NA
    cat(sprintf("Rep %d: ERROR - %s\n", i, conditionMessage(e)))
  })

  if (i %% 50 == 0) cat(sprintf("  %d/%d done\n", i, n_reps))
}

valid <- !is.na(mu_hat)
cat(sprintf("\n=== Analytical Gradient with W' term (bounds=5): %d/%d valid reps ===\n", sum(valid), n_reps))

mu_v <- mu_hat[valid]
q1 <- quantile(mu_v, 0.25); q3 <- quantile(mu_v, 0.75); iqr <- q3 - q1
is_outlier <- mu_v < (q1 - 3*iqr) | mu_v > (q3 + 3*iqr)
cat(sprintf("Outliers: %d/%d (%.1f%%)\n", sum(is_outlier), sum(valid), 100*mean(is_outlier)))
cat(sprintf("Bias: %.4f, RMSE: %.4f, SD: %.4f\n",
            mean(mu_v) - mu_true, sqrt(mean((mu_v - mu_true)^2)), sd(mu_v)))

cat(sprintf("\nBaseline (numerical grad, no W' term):\n"))
cat(sprintf("  Outliers: 3/200 (1.5%%)\n"))
cat(sprintf("  Bias: 0.2779, RMSE: 0.3596, SD: 0.2288\n"))
cat(sprintf("  Convergence: 57/200 (28%%)\n"))
cat(sprintf("  Iterations: mean=40.7\n"))

valid_idx <- which(valid)
outlier_orig_idx <- valid_idx[is_outlier]
normal_orig_idx <- valid_idx[!is_outlier]

# Convergence
cat(sprintf("\n=== CONVERGENCE ===\n"))
cat(sprintf("  Overall converged: %d/%d (%.0f%%)\n", sum(converged[valid]), sum(valid), 100*mean(converged[valid])))

# Iterations
cat(sprintf("\n=== ITERATIONS ===\n"))
cat(sprintf("  Mean=%.1f, median=%.0f, range=[%d, %d]\n",
            mean(iterations[valid]), median(iterations[valid]),
            min(iterations[valid]), max(iterations[valid])))

# GMM objective
cat(sprintf("\n=== GMM OBJECTIVE ===\n"))
cat(sprintf("  Mean=%.6f, median=%.6f, range=[%.6f, %.6f]\n",
            mean(objective[valid]), median(objective[valid]),
            min(objective[valid]), max(objective[valid])))

# Alpha
cat(sprintf("\n=== ALPHA COEFFICIENTS ===\n"))
alpha_all <- do.call(rbind, alpha_mat[valid])
anames <- names(alpha_mat[[valid_idx[1]]])
if (is.null(anames)) anames <- paste0("alpha", 1:ncol(alpha_all))

for (k in 1:ncol(alpha_all)) {
  cat(sprintf("  %s: mean=%7.3f, sd=%.3f\n", anames[k],
              mean(alpha_all[,k]), sd(alpha_all[,k])))
}

# PS and IPW
cat(sprintf("\n=== PROPENSITY SCORES ===\n"))
cat(sprintf("  min(PS): mean=%.6f, min=%.6f\n",
            mean(min_ps[valid]), min(min_ps[valid])))
cat(sprintf("  max(r/PS): mean=%.1f, max=%.1f\n",
            mean(max_ipw[valid]), max(max_ipw[valid])))

# Correlations
cat(sprintf("\n=== CORRELATIONS ===\n"))
cat(sprintf("  Cor(max_ipw, |mu-median|) = %.3f\n",
            cor(max_ipw[valid], abs(mu_v - median(mu_v)))))

# Outlier details
if (sum(is_outlier) > 0) {
  cat(sprintf("\n=== OUTLIER REPS ===\n"))
  for (ii in seq_along(outlier_orig_idx)) {
    k <- outlier_orig_idx[ii]
    cat(sprintf("  Rep %3d: mu=%8.4f, obj=%.6f, conv=%s, iter=%d, min_ps=%.6f, max_ipw=%6.1f, alpha=[%s]\n",
                k, mu_hat[k], objective[k], converged[k], iterations[k],
                min_ps[k], max_ipw[k],
                paste(round(alpha_mat[[k]], 3), collapse=", ")))
  }
}

cat("\nDone!\n")
