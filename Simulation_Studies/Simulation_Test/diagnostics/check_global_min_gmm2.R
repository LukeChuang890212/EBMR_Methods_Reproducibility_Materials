## Check ensemble GMM: compare objectives with multi-start
## Also test the original multi-start approach (commented out in ensemble)

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

# ---- Test 1: Compare GMM objectives ----
cat("=== Ensemble: Compare GMM Objectives (default vs random init) ===\n\n")

set.seed(123)
n_reps <- 10

for (i in 1:n_reps) {
  dat <- setting4.B1(500)

  tryCatch({
    ebmr <- EPS$new("y", ps_spec, dat, W_func, method = "GN")
    J <- length(ebmr$ps_fit.list)

    # Check dimensions on first rep
    if (i == 1) {
      r_tmp <- ebmr$EBMR_IPW(h_nu = h_nu_func, method = "GN")
      h_dim <- ncol(r_tmp$ensemble_fit$h_x)
      cat(sprintf("Ensemble: h_dim=%d (esteq), J=%d (params)\n", h_dim, J))
      if (h_dim > J) cat("  -> OVER-IDENTIFIED\n\n")
      else cat("  -> JUST-IDENTIFIED\n\n")
    }

    # Default init (zeros)
    r_def <- ebmr$EBMR_IPW(h_nu = h_nu_func, method = "GN")
    obj_def <- r_def$ensemble_fit$gmm_fit$opt$objective

    # Try J identity-like inits (the commented-out approach)
    best_obj <- obj_def
    best_mu <- r_def$mu_ipw
    best_nu <- r_def$nu.hat
    best_src <- "default(0)"

    init_mat <- diag(J) * 0.95
    for (j in 1:J) {
      r_j <- ebmr$EBMR_IPW(h_nu = h_nu_func, nu_init = init_mat[j, ], method = "GN")
      obj_j <- r_j$ensemble_fit$gmm_fit$opt$objective
      if (obj_j < best_obj - 1e-10) {
        best_obj <- obj_j
        best_mu <- r_j$mu_ipw
        best_nu <- r_j$nu.hat
        best_src <- sprintf("init_%d", j)
      }
    }

    # Also try random inits
    for (s in 1:5) {
      r_s <- ebmr$EBMR_IPW(h_nu = h_nu_func, nu_init = rnorm(J, sd = 0.5), method = "GN")
      obj_s <- r_s$ensemble_fit$gmm_fit$opt$objective
      if (obj_s < best_obj - 1e-10) {
        best_obj <- obj_s
        best_mu <- r_s$mu_ipw
        best_nu <- r_s$nu.hat
        best_src <- sprintf("random_%d", s)
      }
    }

    if (best_src == "default(0)") {
      cat(sprintf("Rep %2d: default is best. obj=%.6e, mu=%.4f, nu=%s\n",
                  i, obj_def, r_def$mu_ipw, paste(round(r_def$nu.hat, 4), collapse=", ")))
    } else {
      cat(sprintf("Rep %2d: BETTER found by %s!\n", i, best_src))
      cat(sprintf("  default: obj=%.6e, mu=%.4f, nu=%s\n",
                  obj_def, r_def$mu_ipw, paste(round(r_def$nu.hat, 4), collapse=", ")))
      cat(sprintf("  best:    obj=%.6e, mu=%.4f, nu=%s\n",
                  best_obj, best_mu, paste(round(best_nu, 4), collapse=", ")))
    }
  }, error = function(e) {
    cat(sprintf("Rep %2d: ERROR - %s\n", i, conditionMessage(e)))
  })
}

cat("\nDone!\n")
