setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
source("Data_Generation.r")
source("config/scenarios.R")

old_wd <- getwd()
setwd("../EBMRalgorithmFast4")
source("R/EBMRAlgorithm.r")
setwd(old_wd)

library(numDeriv)

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
methods <- c("iterated", "gd", "newton")

# Generate all datasets first
set.seed(42)
datasets <- vector("list", n_reps)
for (i in 1:n_reps) {
  datasets[[i]] <- setting3.B1(n_val)
}

# Storage
results <- list()
for (m in methods) {
  results[[m]] <- list(
    mu_hat = numeric(n_reps),
    alpha_mat = vector("list", n_reps),
    converged = logical(n_reps),
    iterations = integer(n_reps),
    objective = numeric(n_reps),
    min_ps = numeric(n_reps),
    max_ipw = numeric(n_reps),
    phase2_info = vector("list", n_reps)
  )
}

for (m in methods) {
  cat(sprintf("\n========== METHOD: %s ==========\n", m))
  t_start <- proc.time()

  for (i in 1:n_reps) {
    dat <- datasets[[i]]

    tryCatch({
      ebmr <- EBMRAlgorithmFast4$new("y", ps_spec_m2, dat, W_func, bounds = 5, gmm_method = m)

      fit <- ebmr$ps_fit.list[[1]]
      alpha_hat <- fit$coefficients
      ps_vals <- fit$fitted.values
      gmm_fit <- fit$gmm_fit

      results[[m]]$mu_hat[i] <- mean(dat$r / ps_vals * dat$y)
      results[[m]]$alpha_mat[[i]] <- alpha_hat
      results[[m]]$converged[i] <- gmm_fit$opt$converged
      results[[m]]$iterations[i] <- gmm_fit$opt$iterations
      results[[m]]$objective[i] <- gmm_fit$opt$objective
      results[[m]]$min_ps[i] <- min(ps_vals)
      results[[m]]$max_ipw[i] <- max(dat$r / ps_vals)
      results[[m]]$phase2_info[[i]] <- gmm_fit$opt$phase2_info

    }, error = function(e) {
      results[[m]]$mu_hat[i] <<- NA
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

  valid <- !is.na(res$mu_hat)
  cat(sprintf("Valid reps: %d/%d\n", sum(valid), n_reps))
  if (sum(valid) == 0) next

  mu_v <- res$mu_hat[valid]
  q1 <- quantile(mu_v, 0.25); q3 <- quantile(mu_v, 0.75); iqr <- q3 - q1
  is_outlier <- mu_v < (q1 - 3*iqr) | mu_v > (q3 + 3*iqr)
  cat(sprintf("Outliers: %d/%d (%.1f%%)\n", sum(is_outlier), sum(valid), 100*mean(is_outlier)))
  cat(sprintf("Bias: %.4f, RMSE: %.4f, SD: %.4f\n",
              mean(mu_v) - mu_true, sqrt(mean((mu_v - mu_true)^2)), sd(mu_v)))

  valid_idx <- which(valid)

  cat(sprintf("\n  Convergence: %d/%d (%.0f%%)\n",
              sum(res$converged[valid]), sum(valid), 100*mean(res$converged[valid])))
  cat(sprintf("  Iterations: mean=%.1f, range=[%d, %d]\n",
              mean(res$iterations[valid]), min(res$iterations[valid]), max(res$iterations[valid])))
  cat(sprintf("  Objective: mean=%.6f, median=%.6f, range=[%.6f, %.6f]\n",
              mean(res$objective[valid]), median(res$objective[valid]),
              min(res$objective[valid]), max(res$objective[valid])))

  # Alpha
  alpha_all <- do.call(rbind, res$alpha_mat[valid])
  anames <- names(res$alpha_mat[[valid_idx[1]]])
  if (is.null(anames)) anames <- paste0("alpha", 1:ncol(alpha_all))
  cat(sprintf("\n  Alpha coefficients:\n"))
  for (k in 1:ncol(alpha_all)) {
    cat(sprintf("    %s: mean=%7.3f, sd=%.3f\n", anames[k], mean(alpha_all[,k]), sd(alpha_all[,k])))
  }

  cat(sprintf("\n  min(PS): mean=%.6f, min=%.6f\n", mean(res$min_ps[valid]), min(res$min_ps[valid])))
  cat(sprintf("  max(r/PS): mean=%.1f, max=%.1f\n", mean(res$max_ipw[valid]), max(res$max_ipw[valid])))
  cat(sprintf("  Cor(max_ipw, |mu-median|) = %.3f\n", cor(res$max_ipw[valid], abs(mu_v - median(mu_v)))))

  # Phase 2 info for gd/newton
  if (m %in% c("gd", "newton")) {
    p2_valid <- res$phase2_info[valid]
    p2_valid <- p2_valid[!sapply(p2_valid, is.null)]
    if (length(p2_valid) > 0) {
      shifts <- sapply(p2_valid, function(x) x$shift)
      focs <- sapply(p2_valid, function(x) x$foc_max)
      p2iters <- sapply(p2_valid, function(x) x$phase2_iterations)
      Q_before <- sapply(p2_valid, function(x) x$Q_twostep)
      Q_after <- sapply(p2_valid, function(x) x$Q_final)
      cat(sprintf("\n  Phase 2 details:\n"))
      cat(sprintf("    Phase2 iterations: mean=%.1f, range=[%d, %d]\n",
                  mean(p2iters), min(p2iters), max(p2iters)))
      cat(sprintf("    Shift from 2-step: mean=%.6f, max=%.6f\n", mean(shifts), max(shifts)))
      cat(sprintf("    FOC max|grad|:     mean=%.6f, max=%.6f\n", mean(focs), max(focs)))
      cat(sprintf("    Q(2-step):   mean=%.6f, range=[%.6f, %.6f]\n",
                  mean(Q_before), min(Q_before), max(Q_before)))
      cat(sprintf("    Q(final):    mean=%.6f, range=[%.6f, %.6f]\n",
                  mean(Q_after), min(Q_after), max(Q_after)))
      cat(sprintf("    Q decreased: %d/%d\n", sum(Q_after < Q_before), length(Q_after)))
    }
  }

  # Outlier details
  outlier_orig_idx <- valid_idx[is_outlier]
  if (sum(is_outlier) > 0) {
    cat(sprintf("\n  Outlier reps:\n"))
    for (ii in seq_along(outlier_orig_idx)) {
      k <- outlier_orig_idx[ii]
      cat(sprintf("    Rep %3d: mu=%8.4f, obj=%.6f, iter=%d, max_ipw=%6.1f\n",
                  k, res$mu_hat[k], res$objective[k], res$iterations[k], res$max_ipw[k]))
    }
  }
}

# Cross-method comparison
cat("\n\n########## CROSS-METHOD COMPARISON ##########\n")
for (i in 1:n_reps) {
  mus <- sapply(methods, function(m) results[[m]]$mu_hat[i])
  if (any(is.na(mus))) next
  spread <- max(mus) - min(mus)
  if (spread > 0.1) {
    cat(sprintf("Rep %3d: iterated=%.4f, gd=%.4f, newton=%.4f (spread=%.4f)\n",
                i, mus[1], mus[2], mus[3], spread))
  }
}

cat("\nDone!\n")
