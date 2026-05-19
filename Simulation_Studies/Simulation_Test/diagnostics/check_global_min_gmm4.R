## Compare ensemble: GN multi-start vs L-BFGS-B with [-20,20] bounds
## Does L-BFGS-B find better solutions?

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

cat("=== GN multi-start vs L-BFGS-B comparison ===\n\n")

set.seed(123)
n_reps <- 20
n_gn_better <- 0
n_lbfgs_better <- 0
n_same <- 0
n_err <- 0

for (i in 1:n_reps) {
  dat <- setting4.B1(500)

  tryCatch({
    # GN multi-start (current default)
    ebmr_gn <- EPS$new("y", ps_spec, dat, W_func, method = "GN")
    r_gn <- ebmr_gn$EBMR_IPW(h_nu = h_nu_func, method = "GN")
    obj_gn <- r_gn$ensemble_fit$gmm_fit$opt$objective

    # L-BFGS-B (fallback method)
    ebmr_lb <- EPS$new("y", ps_spec, dat, W_func, method = "L-BFGS-B")
    r_lb <- ebmr_lb$EBMR_IPW(h_nu = h_nu_func, method = "L-BFGS-B")
    obj_lb <- r_lb$ensemble_fit$gmm_fit$opt$objective

    mu_diff <- abs(r_gn$mu_ipw - r_lb$mu_ipw)
    obj_ratio <- obj_gn / max(obj_lb, 1e-15)

    if (mu_diff < 0.01) {
      n_same <- n_same + 1
      cat(sprintf("Rep %2d: SAME mu=%.4f (GN obj=%.2e, LB obj=%.2e)\n",
                  i, r_gn$mu_ipw, obj_gn, obj_lb))
    } else if (obj_gn < obj_lb - 1e-10) {
      n_gn_better <- n_gn_better + 1
      cat(sprintf("Rep %2d: GN BETTER  mu_gn=%.4f(obj=%.2e) vs mu_lb=%.4f(obj=%.2e)\n",
                  i, r_gn$mu_ipw, obj_gn, r_lb$mu_ipw, obj_lb))
    } else {
      n_lbfgs_better <- n_lbfgs_better + 1
      cat(sprintf("Rep %2d: LB BETTER  mu_gn=%.4f(obj=%.2e) vs mu_lb=%.4f(obj=%.2e)\n",
                  i, r_gn$mu_ipw, obj_gn, r_lb$mu_ipw, obj_lb))
    }
    cat(sprintf("  nu_gn=%s\n", paste(round(r_gn$nu.hat, 4), collapse=", ")))
    cat(sprintf("  nu_lb=%s\n", paste(round(r_lb$nu.hat, 4), collapse=", ")))

  }, error = function(e) {
    n_err <<- n_err + 1
    cat(sprintf("Rep %2d: ERROR - %s\n", i, conditionMessage(e)))
  })
}

cat(sprintf("\n=== Summary ===\n"))
cat(sprintf("Same: %d, GN better: %d, L-BFGS-B better: %d, Errors: %d\n",
            n_same, n_gn_better, n_lbfgs_better, n_err))
cat("\nDone!\n")
