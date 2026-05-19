## Test L-BFGS-B with different bounds: 20, 10, 5, 3
## Which bound best balances outlier reduction vs estimation quality?

setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
source("Data_Generation.r")
source("config/scenarios.R")

old_wd <- getwd()
setwd("../EBMRalgorithmFast4")
source("R/EBMRAlgorithm.r")
setwd(old_wd)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

ps_spec <- PS_SPECS[["9"]]
ps_spec$alpha_init.list <- list(NULL, NULL, NULL)
ps_spec$outcome <- "y"

h_nu_func <- function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1,
                                  z2 = dat$z2, u1_u2 = dat$u1 * dat$u2)

# Compute true mean
set.seed(999)
big <- setting4.A1(1e6)
mu_true <- mean(big$y)
rm(big); gc()
cat(sprintf("mu_true = %.4f\n\n", mu_true))

n_reps <- 50
n_val <- 500

for (bound in c(20, 10, 5, 3)) {
  cat(sprintf("\n========== bounds = [-%d, %d] ==========\n", bound, bound))

  set.seed(123)
  mu_hat <- numeric(n_reps)
  has_extreme_coef <- logical(n_reps)
  max_coef_all <- numeric(n_reps)
  min_ens_ps <- numeric(n_reps)
  max_ipw_wt <- numeric(n_reps)
  gmm_obj_ps <- matrix(NA, n_reps, 3)

  for (i in 1:n_reps) {
    dat <- setting4.B1(n_val)

    tryCatch({
      ebmr <- EBMRAlgorithmFast4$new("y", ps_spec, dat, W_func,
                                      method = "L-BFGS-B", bounds = bound)
      result <- ebmr$EBMR_IPW(h_nu = h_nu_func, method = "L-BFGS-B")

      mu_hat[i] <- result$mu_ipw
      nu_hat <- result$nu.hat
      w_hat <- nu_hat^2 / sum(nu_hat^2)

      max_c <- 0
      extreme_this <- FALSE
      for (j in 1:3) {
        fit <- ebmr$ps_fit.list[[j]]
        alpha_hat <- fit$coefficients
        if (max(abs(alpha_hat)) > 5) extreme_this <- TRUE
        max_c <- max(max_c, max(abs(alpha_hat)))
        gmm_obj_ps[i, j] <- fit$opt$objective
      }

      has_extreme_coef[i] <- extreme_this
      max_coef_all[i] <- max_c

      ps_vals <- sapply(1:3, function(j) ebmr$ps_fit.list[[j]]$fitted.values)
      ens_ps <- ps_vals %*% w_hat
      min_ens_ps[i] <- min(ens_ps)
      ipw_wts <- dat$r / ens_ps
      max_ipw_wt[i] <- max(ipw_wts[dat$r == 1])

    }, error = function(e) {
      mu_hat[i] <<- NA
      has_extreme_coef[i] <<- NA
    })

    if (i %% 10 == 0) cat(sprintf("  %d/%d done\n", i, n_reps))
  }

  valid <- !is.na(mu_hat)
  resid <- mu_hat[valid] - mu_true
  q1 <- quantile(resid, 0.25)
  q3 <- quantile(resid, 0.75)
  iqr <- q3 - q1
  is_outlier <- abs(resid) > 3 * iqr

  cat(sprintf("\nResults (%d/%d valid):\n", sum(valid), n_reps))
  cat(sprintf("  Outliers: %d/%d (%.1f%%)\n", sum(is_outlier), sum(valid), 100*mean(is_outlier)))
  cat(sprintf("  Bias: %.4f, RMSE: %.4f\n", mean(resid), sqrt(mean(resid^2))))
  cat(sprintf("  Extreme coefs: %d/%d (%.1f%%)\n",
              sum(has_extreme_coef[valid]), sum(valid), 100*mean(has_extreme_coef[valid])))
  cat(sprintf("  max|coef|: mean=%.2f, med=%.2f, max=%.2f\n",
              mean(max_coef_all[valid]), median(max_coef_all[valid]), max(max_coef_all[valid])))
  cat(sprintf("  max_ipw: mean=%.2f, med=%.2f, max=%.2f\n",
              mean(max_ipw_wt[valid]), median(max_ipw_wt[valid]), max(max_ipw_wt[valid])))
  cat(sprintf("  min_ens_ps: mean=%.4f, med=%.4f, min=%.6f\n",
              mean(min_ens_ps[valid]), median(min_ens_ps[valid]), min(min_ens_ps[valid])))
  cat(sprintf("  PS obj (mean): M1=%.4f, M2=%.4f, M3=%.4f\n",
              mean(gmm_obj_ps[valid,1], na.rm=T),
              mean(gmm_obj_ps[valid,2], na.rm=T),
              mean(gmm_obj_ps[valid,3], na.rm=T)))
}

cat("\nDone!\n")
