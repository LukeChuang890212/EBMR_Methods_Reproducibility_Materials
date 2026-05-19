## Check second-order conditions at GN solutions
## Q(alpha) = G(alpha)' W G(alpha)
## FOC: dQ/dalpha = 2 Gamma' W G = 0
## SOC: d2Q/dalpha2 = 2(Gamma'W Gamma + higher-order terms) is positive definite
##
## Also check: do misspecified models produce extreme coefficients?

setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
source("Data_Generation.r")
source("config/scenarios.R")

old_wd <- getwd()
setwd("../EBMRalgorithmFast4")
source("R/EBMRAlgorithm.r")
setwd(old_wd)

library(numDeriv)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

ps_spec <- PS_SPECS[["9"]]
ps_spec$alpha_init.list <- list(NULL, NULL, NULL)
ps_spec$outcome <- "y"

h_nu_func <- function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1,
                                  z2 = dat$z2, u1_u2 = dat$u1 * dat$u2)

cat("=== Second-Order Conditions Check ===\n\n")

set.seed(123)
n_reps <- 20
n_saddle_ps <- 0
n_saddle_ens <- 0
n_extreme_coef <- 0

for (i in 1:n_reps) {
  dat <- setting4.B1(500)
  n <- nrow(dat)

  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", ps_spec, dat, W_func, method = "GN")

    cat(sprintf("--- Rep %d ---\n", i))

    # Check each PS model
    for (j in 1:3) {
      fit <- ebmr$ps_fit.list[[j]]
      alpha_hat <- fit$coefficients
      g_mat <- fit$gmm_fit$g.matrix
      W_hat <- fit$gmm_fit$W.hat
      Gamma_hat <- fit$gmm_fit$Gamma.hat

      # Check FOC: Gamma' W G should be ~0
      G_hat <- colMeans(g_mat)
      foc <- as.vector(t(Gamma_hat) %*% W_hat %*% G_hat)
      max_foc <- max(abs(foc))

      # Check SOC: Gamma' W Gamma should be positive definite
      GtWG <- t(Gamma_hat) %*% W_hat %*% Gamma_hat
      eig_vals <- eigen(GtWG, symmetric = TRUE, only.values = TRUE)$values
      min_eig <- min(eig_vals)
      is_pd <- min_eig > 0

      # Check for extreme coefficients
      max_coef <- max(abs(alpha_hat))

      cat(sprintf("  PS Model %d: max|FOC|=%.2e, min_eig(GtWG)=%.4f, PD=%s, max|coef|=%.2f, obj=%.2e\n",
                  j, max_foc, min_eig, is_pd, max_coef, fit$opt$objective))

      if (!is_pd) {
        n_saddle_ps <- n_saddle_ps + 1
        cat(sprintf("    ** SADDLE POINT! Eigenvalues: %s\n",
                    paste(round(eig_vals, 4), collapse = ", ")))
      }
      if (max_coef > 5) {
        n_extreme_coef <- n_extreme_coef + 1
        cat(sprintf("    ** EXTREME COEF: %s\n",
                    paste(round(alpha_hat, 4), collapse = ", ")))
      }
    }

    # Check ensemble
    result <- ebmr$EBMR_IPW(h_nu = h_nu_func, method = "GN")
    ens_fit <- result$ensemble_fit$gmm_fit
    nu_hat <- result$nu.hat

    if (!all(is.na(ens_fit$Gamma.hat))) {
      Gamma_ens <- ens_fit$Gamma.hat
      W_ens <- ens_fit$W.hat
      g_ens <- ens_fit$g.matrix
      G_ens <- colMeans(g_ens)

      foc_ens <- as.vector(t(Gamma_ens) %*% W_ens %*% G_ens)
      GtWG_ens <- t(Gamma_ens) %*% W_ens %*% Gamma_ens
      eig_ens <- eigen(GtWG_ens, symmetric = TRUE, only.values = TRUE)$values
      min_eig_ens <- min(eig_ens)

      cat(sprintf("  Ensemble: max|FOC|=%.2e, min_eig=%.4f, PD=%s, obj=%.2e, mu=%.4f\n",
                  max(abs(foc_ens)), min_eig_ens, min_eig_ens > 0,
                  ens_fit$opt$objective, result$mu_ipw))
      cat(sprintf("    nu=%s, w=%s\n",
                  paste(round(nu_hat, 4), collapse=", "),
                  paste(round(nu_hat^2/sum(nu_hat^2), 4), collapse=", ")))

      if (min_eig_ens <= 0) {
        n_saddle_ens <- n_saddle_ens + 1
        cat(sprintf("    ** ENSEMBLE SADDLE POINT! Eigenvalues: %s\n",
                    paste(round(eig_ens, 4), collapse = ", ")))
      }
    }

    # Check: what are the individual PS values at the solution?
    ps_vals <- sapply(1:3, function(j) ebmr$ps_fit.list[[j]]$fitted.values)
    cat(sprintf("  PS ranges: M1=[%.4f,%.4f], M2=[%.4f,%.4f], M3=[%.4f,%.4f]\n",
                min(ps_vals[,1]), max(ps_vals[,1]),
                min(ps_vals[,2]), max(ps_vals[,2]),
                min(ps_vals[,3]), max(ps_vals[,3])))

    # Ensemble PS
    w_hat <- nu_hat^2 / sum(nu_hat^2)
    ens_ps <- ps_vals %*% w_hat
    cat(sprintf("  Ensemble PS: [%.6f, %.4f], mean=%.4f\n",
                min(ens_ps), max(ens_ps), mean(ens_ps)))

    # IPW weights
    ipw_weights <- dat$r / ens_ps
    cat(sprintf("  IPW weights (r/pi): max=%.2f, P(>10)=%.3f, P(>50)=%.3f\n\n",
                max(ipw_weights[dat$r == 1]),
                mean(ipw_weights[dat$r == 1] > 10),
                mean(ipw_weights[dat$r == 1] > 50)))

  }, error = function(e) {
    cat(sprintf("Rep %d: ERROR - %s\n\n", i, conditionMessage(e)))
  })
}

cat(sprintf("\n=== Summary ===\n"))
cat(sprintf("PS saddle points: %d/%d\n", n_saddle_ps, n_reps * 3))
cat(sprintf("Ensemble saddle points: %d/%d\n", n_saddle_ens, n_reps))
cat(sprintf("Extreme PS coefficients (|coef|>5): %d/%d\n", n_extreme_coef, n_reps * 3))
cat("\nDone!\n")
