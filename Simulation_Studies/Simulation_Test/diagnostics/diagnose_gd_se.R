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

# Generate datasets
set.seed(42)
n_reps <- 200
n_val <- 2000
datasets <- vector("list", n_reps)
for (i in 1:n_reps) datasets[[i]] <- setting3.B1(n_val)

methods <- c("iterated", "gd")

results <- list()
for (m in methods) {
  results[[m]] <- list(
    mu_hat = numeric(n_reps), se_hat = numeric(n_reps),
    alpha_mat = vector("list", n_reps),
    converged = logical(n_reps), iterations = integer(n_reps),
    newton_iters = integer(n_reps),
    objective = numeric(n_reps), final_grad_norm = numeric(n_reps),
    min_ps = numeric(n_reps), max_ipw = numeric(n_reps)
  )
}

for (m in methods) {
  cat(sprintf("\n========== METHOD: %s ==========\n", m))
  t_start <- proc.time()

  for (i in 1:n_reps) {
    dat <- datasets[[i]]
    tryCatch({
      ebmr <- EBMRAlgorithmFast4$new("y", ps_spec_m2, dat, W_func, bounds = 5, gmm_method = m)

      # Call EBMR_IPW to get both mu and SE
      ipw_result <- ebmr$EBMR_IPW(h_nu = h_nu_func)

      fit <- ebmr$ps_fit.list[[1]]
      ps_vals <- fit$fitted.values
      gmm_fit <- fit$gmm_fit

      results[[m]]$mu_hat[i] <- ipw_result$mu_ipw
      results[[m]]$se_hat[i] <- ipw_result$se_ipw
      results[[m]]$alpha_mat[[i]] <- fit$coefficients
      results[[m]]$converged[i] <- gmm_fit$opt$converged
      results[[m]]$iterations[i] <- gmm_fit$opt$iterations
      results[[m]]$objective[i] <- gmm_fit$opt$objective
      results[[m]]$final_grad_norm[i] <- ifelse(is.null(gmm_fit$opt$final_grad_norm), NA, gmm_fit$opt$final_grad_norm)
      results[[m]]$newton_iters[i] <- ifelse(is.null(gmm_fit$opt$newton_iters), 0L, gmm_fit$opt$newton_iters)
      results[[m]]$min_ps[i] <- min(ps_vals)
      results[[m]]$max_ipw[i] <- max(dat$r / ps_vals)
    }, error = function(e) {
      results[[m]]$mu_hat[i] <<- NA
      results[[m]]$se_hat[i] <<- NA
      results[[m]]$alpha_mat[[i]] <<- rep(NA, 4)
      results[[m]]$converged[i] <<- NA
      cat(sprintf("Rep %d: ERROR - %s\n", i, conditionMessage(e)))
    })
    if (i %% 50 == 0) cat(sprintf("  %d/%d done\n", i, n_reps))
  }
  elapsed <- (proc.time() - t_start)[3]
  cat(sprintf("  Elapsed: %.1f sec\n", elapsed))
}

# Report
for (m in methods) {
  cat(sprintf("\n\n########## METHOD: %s ##########\n", toupper(m)))
  res <- results[[m]]
  valid <- !is.na(res$mu_hat) & !is.na(res$se_hat)
  cat(sprintf("Valid reps: %d/%d\n", sum(valid), n_reps))
  if (sum(valid) == 0) next

  mu_v <- res$mu_hat[valid]
  se_v <- res$se_hat[valid]

  q1 <- quantile(mu_v, 0.25); q3 <- quantile(mu_v, 0.75); iqr <- q3 - q1
  is_outlier <- mu_v < (q1 - 3*iqr) | mu_v > (q3 + 3*iqr)
  n_outlier <- sum(is_outlier)
  cat(sprintf("Outliers removed: %d/%d (%.1f%%)\n", n_outlier, sum(valid), 100*mean(is_outlier)))

  # Remove outliers
  mu_v <- mu_v[!is_outlier]
  se_v <- se_v[!is_outlier]
  n_kept <- length(mu_v)

  # Bias, RMSE, SD
  empirical_sd <- sd(mu_v)
  cat(sprintf("Bias: %.4f, RMSE: %.4f, SD(empirical): %.4f\n",
              mean(mu_v) - mu_true, sqrt(mean((mu_v - mu_true)^2)), empirical_sd))

  # SE summary
  cat(sprintf("\n  SE summary:\n"))
  cat(sprintf("    mean(SE): %.4f, median(SE): %.4f, sd(SE): %.4f\n",
              mean(se_v), median(se_v), sd(se_v)))
  cat(sprintf("    range(SE): [%.4f, %.4f]\n", min(se_v), max(se_v)))
  cat(sprintf("    Empirical SD:  %.4f\n", empirical_sd))
  cat(sprintf("    mean(SE)/SD:   %.4f  (ideal = 1.0)\n", mean(se_v) / empirical_sd))
  cat(sprintf("    median(SE)/SD: %.4f  (ideal = 1.0)\n", median(se_v) / empirical_sd))

  # Coverage: 95% CI
  ci_lower <- mu_v - 1.96 * se_v
  ci_upper <- mu_v + 1.96 * se_v
  covers <- (ci_lower <= mu_true) & (mu_true <= ci_upper)
  cat(sprintf("    95%% CI coverage: %d/%d (%.1f%%)\n", sum(covers), n_kept, 100*mean(covers)))

  # Use non-outlier indices for remaining summaries
  valid_idx <- which(valid)
  kept_idx <- valid_idx[!is_outlier]

  # Convergence & iterations
  cat(sprintf("\n  Convergence: %d/%d (%.0f%%)\n",
              sum(res$converged[kept_idx]), n_kept, 100*mean(res$converged[kept_idx])))
  cat(sprintf("  Iterations: mean=%.1f, range=[%d, %d]\n",
              mean(res$iterations[kept_idx]), min(res$iterations[kept_idx]), max(res$iterations[kept_idx])))
  cat(sprintf("  Objective: mean=%.6f, median=%.6f\n",
              mean(res$objective[kept_idx]), median(res$objective[kept_idx])))
  fgn <- res$final_grad_norm[kept_idx]
  fgn <- fgn[!is.na(fgn)]
  if (length(fgn) > 0) {
    cat(sprintf("  Final grad norm: mean=%.6f, median=%.6f, max=%.6f\n",
                mean(fgn), median(fgn), max(fgn)))
    ni <- res$newton_iters[kept_idx]
    cat(sprintf("  Newton iters: mean=%.1f, range=[%d, %d]\n",
                mean(ni), min(ni), max(ni)))
  }

  # Alpha
  alpha_all <- do.call(rbind, res$alpha_mat[kept_idx])
  anames <- names(res$alpha_mat[[kept_idx[1]]])
  if (is.null(anames)) anames <- paste0("alpha", 1:ncol(alpha_all))
  cat(sprintf("\n  Alpha coefficients:\n"))
  for (k in 1:ncol(alpha_all)) {
    cat(sprintf("    %s: mean=%7.3f, sd=%.3f\n", anames[k], mean(alpha_all[,k]), sd(alpha_all[,k])))
  }

  cat(sprintf("\n  min(PS): mean=%.6f, min=%.6f\n", mean(res$min_ps[kept_idx]), min(res$min_ps[kept_idx])))
  cat(sprintf("  max(r/PS): mean=%.1f, max=%.1f\n", mean(res$max_ipw[kept_idx]), max(res$max_ipw[kept_idx])))
}

cat("\nDone!\n")
