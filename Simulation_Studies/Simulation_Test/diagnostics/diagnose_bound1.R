setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
source("Data_Generation.r")
source("config/scenarios.R")

old_wd <- getwd()
setwd("../EPS")
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
n_respond <- integer(n_reps)

for (i in 1:n_reps) {
  dat <- setting3.B1(n_val)
  n_respond[i] <- sum(dat$r)

  tryCatch({
    ebmr <- EPS$new("y", ps_spec_m2, dat, W_func, bounds = 3)

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
cat(sprintf("\n=== BOUNDS=1: %d/%d valid reps ===\n", sum(valid), n_reps))

mu_v <- mu_hat[valid]
q1 <- quantile(mu_v, 0.25); q3 <- quantile(mu_v, 0.75); iqr <- q3 - q1
is_outlier <- mu_v < (q1 - 3*iqr) | mu_v > (q3 + 3*iqr)
cat(sprintf("Outliers: %d/%d (%.1f%%)\n", sum(is_outlier), sum(valid), 100*mean(is_outlier)))
cat(sprintf("Bias: %.4f, RMSE: %.4f, SD: %.4f\n",
            mean(mu_v) - mu_true, sqrt(mean((mu_v - mu_true)^2)), sd(mu_v)))

valid_idx <- which(valid)
outlier_orig_idx <- valid_idx[is_outlier]
normal_orig_idx <- valid_idx[!is_outlier]

# Convergence
cat(sprintf("\n=== CONVERGENCE ===\n"))
cat(sprintf("  Overall converged: %d/%d (%.0f%%)\n", sum(converged[valid]), sum(valid), 100*mean(converged[valid])))
if (sum(is_outlier) > 0) {
  cat(sprintf("  Outlier converged:     %d/%d (%.0f%%)\n", sum(converged[outlier_orig_idx]), length(outlier_orig_idx), 100*mean(converged[outlier_orig_idx])))
}
cat(sprintf("  Non-outlier converged: %d/%d (%.0f%%)\n", sum(converged[normal_orig_idx]), length(normal_orig_idx), 100*mean(converged[normal_orig_idx])))

# Iterations
cat(sprintf("\n=== ITERATIONS ===\n"))
if (sum(is_outlier) > 0) {
  cat(sprintf("  Outlier:     mean=%.1f, median=%.0f, range=[%d, %d]\n",
              mean(iterations[outlier_orig_idx]), median(iterations[outlier_orig_idx]),
              min(iterations[outlier_orig_idx]), max(iterations[outlier_orig_idx])))
}
cat(sprintf("  Non-outlier: mean=%.1f, median=%.0f, range=[%d, %d]\n",
            mean(iterations[normal_orig_idx]), median(iterations[normal_orig_idx]),
            min(iterations[normal_orig_idx]), max(iterations[normal_orig_idx])))

# GMM objective
cat(sprintf("\n=== GMM OBJECTIVE ===\n"))
if (sum(is_outlier) > 0) {
  cat(sprintf("  Outlier:     mean=%.6f, median=%.6f, range=[%.6f, %.6f]\n",
              mean(objective[outlier_orig_idx]), median(objective[outlier_orig_idx]),
              min(objective[outlier_orig_idx]), max(objective[outlier_orig_idx])))
}
cat(sprintf("  Non-outlier: mean=%.6f, median=%.6f, range=[%.6f, %.6f]\n",
            mean(objective[normal_orig_idx]), median(objective[normal_orig_idx]),
            min(objective[normal_orig_idx]), max(objective[normal_orig_idx])))

# Alpha coefficients
cat(sprintf("\n=== ALPHA COEFFICIENTS ===\n"))
alpha_all <- do.call(rbind, alpha_mat[valid])
anames <- names(alpha_mat[[valid_idx[1]]])
if (is.null(anames)) anames <- paste0("alpha", 1:ncol(alpha_all))

for (k in 1:ncol(alpha_all)) {
  cat(sprintf("  %s:\n", anames[k]))
  if (sum(is_outlier) > 0) {
    cat(sprintf("    Outlier:     mean=%7.3f, sd=%.3f, range=[%7.3f, %7.3f]\n",
                mean(alpha_all[is_outlier, k]), sd(alpha_all[is_outlier, k]),
                min(alpha_all[is_outlier, k]), max(alpha_all[is_outlier, k])))
  }
  cat(sprintf("    Non-outlier: mean=%7.3f, sd=%.3f, range=[%7.3f, %7.3f]\n",
              mean(alpha_all[!is_outlier, k]), sd(alpha_all[!is_outlier, k]),
              min(alpha_all[!is_outlier, k]), max(alpha_all[!is_outlier, k])))
}

max_alpha <- apply(abs(alpha_all), 1, max)
cat(sprintf("\n  Max|alpha|:\n"))
if (sum(is_outlier) > 0) {
  cat(sprintf("    Outlier:     mean=%.3f, median=%.3f, range=[%.3f, %.3f]\n",
              mean(max_alpha[is_outlier]), median(max_alpha[is_outlier]),
              min(max_alpha[is_outlier]), max(max_alpha[is_outlier])))
}
cat(sprintf("    Non-outlier: mean=%.3f, median=%.3f, range=[%.3f, %.3f]\n",
            mean(max_alpha[!is_outlier]), median(max_alpha[!is_outlier]),
            min(max_alpha[!is_outlier]), max(max_alpha[!is_outlier])))

at_bound <- max_alpha >= 2.99
cat(sprintf("\n  At bound (|alpha|>=2.99):\n"))
if (sum(is_outlier) > 0) {
  cat(sprintf("    Outlier: %d/%d\n", sum(at_bound[is_outlier]), sum(is_outlier)))
}
cat(sprintf("    Non-outlier: %d/%d\n", sum(at_bound[!is_outlier]), sum(!is_outlier)))

# PS and IPW
cat(sprintf("\n=== PROPENSITY SCORES ===\n"))
if (sum(is_outlier) > 0) {
  cat(sprintf("  min(PS) outlier:     mean=%.6f, median=%.6f, min=%.6f\n",
              mean(min_ps[outlier_orig_idx]), median(min_ps[outlier_orig_idx]), min(min_ps[outlier_orig_idx])))
}
cat(sprintf("  min(PS) non-outlier: mean=%.6f, median=%.6f, min=%.6f\n",
            mean(min_ps[normal_orig_idx]), median(min_ps[normal_orig_idx]), min(min_ps[normal_orig_idx])))

cat(sprintf("\n  max(r/PS) (max IPW weight):\n"))
if (sum(is_outlier) > 0) {
  cat(sprintf("    Outlier:     mean=%.1f, median=%.1f, max=%.1f\n",
              mean(max_ipw[outlier_orig_idx]), median(max_ipw[outlier_orig_idx]), max(max_ipw[outlier_orig_idx])))
}
cat(sprintf("    Non-outlier: mean=%.1f, median=%.1f, max=%.1f\n",
            mean(max_ipw[normal_orig_idx]), median(max_ipw[normal_orig_idx]), max(max_ipw[normal_orig_idx])))

# Correlations
cat(sprintf("\n=== CORRELATIONS ===\n"))
cat(sprintf("  Cor(max_ipw, |mu-median|) = %.3f\n", cor(max_ipw[valid], abs(mu_v - median(mu_v)))))
cat(sprintf("  Cor(objective, |mu-median|) = %.3f\n", cor(objective[valid], abs(mu_v - median(mu_v)))))
cat(sprintf("  Cor(iterations, |mu-median|) = %.3f\n", cor(iterations[valid], abs(mu_v - median(mu_v)))))

# Converged vs not
cat(sprintf("\n=== CONVERGED vs NOT ===\n"))
conv_v <- converged[valid]
cat(sprintf("  Converged (%d): mu mean=%.4f, sd=%.4f, max_ipw mean=%.1f\n",
            sum(conv_v), mean(mu_v[conv_v]), sd(mu_v[conv_v]),
            mean(max_ipw[valid][conv_v])))
cat(sprintf("  Not converged (%d): mu mean=%.4f, sd=%.4f, max_ipw mean=%.1f\n",
            sum(!conv_v), mean(mu_v[!conv_v]), sd(mu_v[!conv_v]),
            mean(max_ipw[valid][!conv_v])))

# All outlier reps
if (sum(is_outlier) > 0) {
  cat(sprintf("\n=== ALL OUTLIER REPS ===\n"))
  for (ii in seq_along(outlier_orig_idx)) {
    k <- outlier_orig_idx[ii]
    cat(sprintf("  Rep %3d: mu=%8.4f, obj=%.6f, conv=%s, iter=%d, min_ps=%.6f, max_ipw=%6.1f, alpha=[%s]\n",
                k, mu_hat[k], objective[k], converged[k], iterations[k],
                min_ps[k], max_ipw[k],
                paste(round(alpha_mat[[k]], 3), collapse=", ")))
  }
}

# Iteration distribution
cat(sprintf("\n=== ITERATION DISTRIBUTION ===\n"))
iter_v <- iterations[valid]
breaks <- c(0, 5, 10, 20, 50, 100, 200, 1000)
for (b in 2:length(breaks)) {
  mask <- iter_v >= breaks[b-1] & iter_v < breaks[b]
  if (sum(mask) > 0) {
    cat(sprintf("  iter [%3d,%3d): n=%3d, mu_sd=%.4f, max_ipw_mean=%.1f\n",
                breaks[b-1], breaks[b], sum(mask),
                sd(mu_v[mask]), mean(max_ipw[valid][mask])))
  }
}

cat("\nDone!\n")
