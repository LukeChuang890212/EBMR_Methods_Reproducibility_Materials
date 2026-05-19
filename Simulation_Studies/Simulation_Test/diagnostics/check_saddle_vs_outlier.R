## Does a PS saddle point cause an outlier in mu_ipw?
## Run many reps, flag saddle points and outliers, check correlation.

setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
source("Data_Generation.r")
source("config/scenarios.R")

old_wd <- getwd()
setwd("../EPS")
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

set.seed(123)
n_reps <- 50
n_val <- 500
cat("Starting simulation loop...\n")

# Storage
mu_hat <- numeric(n_reps)
has_saddle <- logical(n_reps)
has_extreme_coef <- logical(n_reps)
min_ens_ps <- numeric(n_reps)
max_ipw_wt <- numeric(n_reps)
max_coef_all <- numeric(n_reps)

for (i in 1:n_reps) {
  dat <- setting4.B1(n_val)

  tryCatch({
    ebmr <- EPS$new("y", ps_spec, dat, W_func, method = "GN")
    result <- ebmr$EBMR_IPW(h_nu = h_nu_func, method = "GN")

    mu_hat[i] <- result$mu_ipw
    nu_hat <- result$nu.hat
    w_hat <- nu_hat^2 / sum(nu_hat^2)

    # Check each PS model for saddle points
    saddle_this <- FALSE
    extreme_this <- FALSE
    max_c <- 0
    for (j in 1:3) {
      fit <- ebmr$ps_fit.list[[j]]
      alpha_hat <- fit$coefficients
      Gamma_hat <- fit$gmm_fit$Gamma.hat
      W_hat <- fit$gmm_fit$W.hat

      if (!any(is.na(Gamma_hat)) && !any(is.na(W_hat))) {
        GtWG <- t(Gamma_hat) %*% W_hat %*% Gamma_hat
        eig_vals <- eigen(GtWG, symmetric = TRUE, only.values = TRUE)$values
        if (min(eig_vals) <= 0) saddle_this <- TRUE
      }
      if (max(abs(alpha_hat)) > 5) extreme_this <- TRUE
      max_c <- max(max_c, max(abs(alpha_hat)))
    }

    has_saddle[i] <- saddle_this
    has_extreme_coef[i] <- extreme_this
    max_coef_all[i] <- max_c

    # Ensemble PS stats
    ps_vals <- sapply(1:3, function(j) ebmr$ps_fit.list[[j]]$fitted.values)
    ens_ps <- ps_vals %*% w_hat
    min_ens_ps[i] <- min(ens_ps)
    ipw_wts <- dat$r / ens_ps
    max_ipw_wt[i] <- max(ipw_wts[dat$r == 1])

  }, error = function(e) {
    mu_hat[i] <<- NA
    has_saddle[i] <<- NA
    has_extreme_coef[i] <<- NA
  })

  if (i %% 10 == 0) cat(sprintf("  %d/%d done\n", i, n_reps))
}

valid <- !is.na(mu_hat)
cat(sprintf("\n=== Results (%d/%d valid) ===\n", sum(valid), n_reps))

# Define outlier: |mu_hat - mu_true| > 3*IQR from median
resid <- mu_hat[valid] - mu_true
q1 <- quantile(resid, 0.25)
q3 <- quantile(resid, 0.75)
iqr <- q3 - q1
is_outlier <- abs(resid) > 3 * iqr
outlier_thresh <- 3 * iqr
cat(sprintf("Outlier threshold: |mu - mu_true| > %.4f (3*IQR)\n", outlier_thresh))
cat(sprintf("Total outliers: %d/%d (%.1f%%)\n\n", sum(is_outlier), sum(valid), 100*mean(is_outlier)))

# Cross-tabulate saddle points vs outliers
saddle_v <- has_saddle[valid]
extreme_v <- has_extreme_coef[valid]

cat("--- Saddle Point vs Outlier ---\n")
cat(sprintf("  Saddle & Outlier:     %d\n", sum(saddle_v & is_outlier)))
cat(sprintf("  Saddle & No outlier:  %d\n", sum(saddle_v & !is_outlier)))
cat(sprintf("  No saddle & Outlier:  %d\n", sum(!saddle_v & is_outlier)))
cat(sprintf("  No saddle & No outlier: %d\n", sum(!saddle_v & !is_outlier)))
cat(sprintf("  P(outlier | saddle) = %.1f%%\n", 100 * mean(is_outlier[saddle_v])))
cat(sprintf("  P(outlier | no saddle) = %.1f%%\n\n", 100 * mean(is_outlier[!saddle_v])))

cat("--- Extreme Coef (|alpha|>5) vs Outlier ---\n")
cat(sprintf("  Extreme & Outlier:     %d\n", sum(extreme_v & is_outlier)))
cat(sprintf("  Extreme & No outlier:  %d\n", sum(extreme_v & !is_outlier)))
cat(sprintf("  No extreme & Outlier:  %d\n", sum(!extreme_v & is_outlier)))
cat(sprintf("  No extreme & No outlier: %d\n", sum(!extreme_v & !is_outlier)))
cat(sprintf("  P(outlier | extreme) = %.1f%%\n", 100 * mean(is_outlier[extreme_v])))
cat(sprintf("  P(outlier | no extreme) = %.1f%%\n\n", 100 * mean(is_outlier[!extreme_v])))

cat("--- Min ensemble PS vs Outlier ---\n")
cat(sprintf("  Outlier reps:    min_ens_ps mean=%.6f, median=%.6f\n",
            mean(min_ens_ps[valid][is_outlier]), median(min_ens_ps[valid][is_outlier])))
cat(sprintf("  Non-outlier reps: min_ens_ps mean=%.6f, median=%.6f\n\n",
            mean(min_ens_ps[valid][!is_outlier]), median(min_ens_ps[valid][!is_outlier])))

cat("--- Max IPW weight vs Outlier ---\n")
cat(sprintf("  Outlier reps:    max_ipw mean=%.2f, median=%.2f\n",
            mean(max_ipw_wt[valid][is_outlier]), median(max_ipw_wt[valid][is_outlier])))
cat(sprintf("  Non-outlier reps: max_ipw mean=%.2f, median=%.2f\n\n",
            mean(max_ipw_wt[valid][!is_outlier]), median(max_ipw_wt[valid][!is_outlier])))

cat("--- Max |coef| vs Outlier ---\n")
cat(sprintf("  Outlier reps:    max|coef| mean=%.2f, median=%.2f\n",
            mean(max_coef_all[valid][is_outlier]), median(max_coef_all[valid][is_outlier])))
cat(sprintf("  Non-outlier reps: max|coef| mean=%.2f, median=%.2f\n",
            mean(max_coef_all[valid][!is_outlier]), median(max_coef_all[valid][!is_outlier])))

cat("\nDone!\n")
