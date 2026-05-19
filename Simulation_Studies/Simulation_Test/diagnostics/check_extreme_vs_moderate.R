## For reps where GN finds extreme alpha, is it really the global min?
## Compare: unconstrained GN vs L-BFGS-B with [-5,5] bounds on each PS model

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

set.seed(123)
n_reps <- 50
n_val <- 500

cat("Comparing unconstrained GN vs constrained L-BFGS-B for each PS model\n\n")

for (i in 1:n_reps) {
  dat <- setting4.B1(n_val)

  tryCatch({
    # Fit with GN (unconstrained)
    ebmr_gn <- EPS$new("y", ps_spec, dat, W_func, method = "GN")

    any_extreme <- FALSE
    for (j in 1:3) {
      fit_gn <- ebmr_gn$ps_fit.list[[j]]
      alpha_gn <- fit_gn$coefficients
      max_coef <- max(abs(alpha_gn))

      if (max_coef > 5) {
        any_extreme <- TRUE
        obj_gn <- fit_gn$opt$objective

        # Now fit same PS model with L-BFGS-B [-5,5]
        ebmr_lb <- EPS$new("y", ps_spec, dat, W_func, method = "L-BFGS-B")
        fit_lb <- ebmr_lb$ps_fit.list[[j]]
        alpha_lb <- fit_lb$coefficients
        obj_lb <- fit_lb$opt$objective

        cat(sprintf("Rep %2d, PS %d: EXTREME\n", i, j))
        cat(sprintf("  GN:     alpha=[%s], max|a|=%.2f, obj=%.6f\n",
                    paste(round(alpha_gn, 3), collapse=", "), max_coef, obj_gn))
        cat(sprintf("  LBFGS:  alpha=[%s], max|a|=%.2f, obj=%.6f\n",
                    paste(round(alpha_lb, 3), collapse=", "), max(abs(alpha_lb)), obj_lb))
        cat(sprintf("  obj_GN %s obj_LB (diff=%.6f)\n\n",
                    ifelse(obj_gn < obj_lb, "<", ifelse(obj_gn > obj_lb, ">", "=")),
                    obj_gn - obj_lb))
      }
    }
  }, error = function(e) {
    cat(sprintf("Rep %2d: ERROR - %s\n", i, conditionMessage(e)))
  })
}

cat("Done!\n")
