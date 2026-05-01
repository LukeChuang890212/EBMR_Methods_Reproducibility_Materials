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

h_nu_func <- function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1,
                                  z2 = dat$z2, u1_u2 = dat$u1 * dat$u2)

set.seed(999)
big <- setting3.A1(1e6)
mu_true <- mean(big$y)
rm(big); gc()
cat(sprintf("mu_true = %.6f\n\n", mu_true))

set.seed(42)
n_reps <- 200
n_val <- 2000
datasets <- vector("list", n_reps)
for (i in 1:n_reps) datasets[[i]] <- setting3.B1(n_val)

bounds_list <- c(5, Inf)

for (b in bounds_list) {
  cat(sprintf("\n========== GD with bounds = %s ==========\n", ifelse(is.finite(b), b, "Inf (none)")))
  mu_hat <- numeric(n_reps)
  se_hat <- numeric(n_reps)
  alpha_mat <- vector("list", n_reps)
  converged <- logical(n_reps)
  iterations <- integer(n_reps)
  objective <- numeric(n_reps)
  grad_norm_vec <- numeric(n_reps)
  min_ps <- numeric(n_reps)
  max_ipw <- numeric(n_reps)

  t_start <- proc.time()
  for (i in 1:n_reps) {
    dat <- datasets[[i]]
    tryCatch({
      ebmr <- EBMRAlgorithmFast4$new("y", ps_spec_m2, dat, W_func, bounds = b, gmm_method = "gd")
      ipw_result <- ebmr$EBMR_IPW(h_nu = h_nu_func)
      fit <- ebmr$ps_fit.list[[1]]
      ps_vals <- fit$fitted.values
      gmm_fit <- fit$gmm_fit

      mu_hat[i] <- ipw_result$mu_ipw
      se_hat[i] <- ipw_result$se_ipw
      alpha_mat[[i]] <- fit$coefficients
      converged[i] <- gmm_fit$opt$converged
      iterations[i] <- gmm_fit$opt$iterations
      objective[i] <- gmm_fit$opt$objective
      grad_norm_vec[i] <- ifelse(is.null(gmm_fit$opt$final_grad_norm), NA, gmm_fit$opt$final_grad_norm)
      min_ps[i] <- min(ps_vals)
      max_ipw[i] <- max(dat$r / ps_vals)
    }, error = function(e) {
      mu_hat[i] <<- NA; se_hat[i] <<- NA
      alpha_mat[[i]] <<- rep(NA, 4); converged[i] <<- NA
      cat(sprintf("Rep %d: ERROR - %s\n", i, conditionMessage(e)))
    })
    if (i %% 50 == 0) cat(sprintf("  %d/%d done\n", i, n_reps))
  }
  elapsed <- (proc.time() - t_start)[3]
  cat(sprintf("  Elapsed: %.1f sec\n", elapsed))

  # Report
  valid <- !is.na(mu_hat) & !is.na(se_hat)
  cat(sprintf("Valid reps: %d/%d\n", sum(valid), n_reps))
  if (sum(valid) == 0) next

  mu_v <- mu_hat[valid]; se_v <- se_hat[valid]
  q1 <- quantile(mu_v, 0.25); q3 <- quantile(mu_v, 0.75); iqr <- q3 - q1
  is_outlier <- mu_v < (q1 - 3*iqr) | mu_v > (q3 + 3*iqr)
  cat(sprintf("Outliers removed: %d/%d (%.1f%%)\n", sum(is_outlier), sum(valid), 100*mean(is_outlier)))

  mu_v <- mu_v[!is_outlier]; se_v <- se_v[!is_outlier]; n_kept <- length(mu_v)
  empirical_sd <- sd(mu_v)
  cat(sprintf("Bias: %.4f, RMSE: %.4f, SD: %.4f\n",
              mean(mu_v) - mu_true, sqrt(mean((mu_v - mu_true)^2)), empirical_sd))
  cat(sprintf("  mean(SE)/SD: %.4f, median(SE)/SD: %.4f\n", mean(se_v)/empirical_sd, median(se_v)/empirical_sd))
  ci_lower <- mu_v - 1.96*se_v; ci_upper <- mu_v + 1.96*se_v
  covers <- (ci_lower <= mu_true) & (mu_true <= ci_upper)
  cat(sprintf("  95%% CI coverage: %d/%d (%.1f%%)\n", sum(covers), n_kept, 100*mean(covers)))

  valid_idx <- which(valid); kept_idx <- valid_idx[!is_outlier]
  cat(sprintf("  Convergence: %d/%d (%.0f%%)\n", sum(converged[kept_idx]), n_kept, 100*mean(converged[kept_idx])))
  cat(sprintf("  Iterations: mean=%.1f, range=[%d, %d]\n",
              mean(iterations[kept_idx]), min(iterations[kept_idx]), max(iterations[kept_idx])))
  fgn <- grad_norm_vec[kept_idx]; fgn <- fgn[!is.na(fgn)]
  if (length(fgn) > 0) cat(sprintf("  Final grad norm: mean=%.6f, median=%.6f, max=%.6f\n", mean(fgn), median(fgn), max(fgn)))

  alpha_all <- do.call(rbind, alpha_mat[kept_idx])
  anames <- names(alpha_mat[[kept_idx[1]]])
  if (is.null(anames)) anames <- paste0("alpha", 1:ncol(alpha_all))
  cat(sprintf("  Alpha coefficients:\n"))
  for (k in 1:ncol(alpha_all)) {
    cat(sprintf("    %s: mean=%7.3f, sd=%.3f, range=[%.3f, %.3f]\n",
                anames[k], mean(alpha_all[,k]), sd(alpha_all[,k]), min(alpha_all[,k]), max(alpha_all[,k])))
  }
  cat(sprintf("  min(PS): mean=%.6f, min=%.6f\n", mean(min_ps[kept_idx]), min(min_ps[kept_idx])))
  cat(sprintf("  max(r/PS): mean=%.1f, max=%.1f\n", mean(max_ipw[kept_idx]), max(max_ipw[kept_idx])))
}

cat("\nDone!\n")
