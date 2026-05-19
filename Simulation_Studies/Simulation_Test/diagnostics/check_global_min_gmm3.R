## Test multi-start ensemble: does it now always find the best solution?

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

cat("=== Multi-start ensemble: stability test ===\n\n")

set.seed(123)
n_reps <- 20
n_diff <- 0

for (i in 1:n_reps) {
  dat <- setting4.B1(500)

  tryCatch({
    ebmr <- EPS$new("y", ps_spec, dat, W_func, method = "GN")
    J <- length(ebmr$ps_fit.list)

    # Default (multi-start picks best)
    r_def <- ebmr$EBMR_IPW(h_nu = h_nu_func, method = "GN")

    # Extra random init (will be added to multi-start candidates)
    r_rnd <- ebmr$EBMR_IPW(h_nu = h_nu_func, nu_init = rnorm(J, sd = 0.5), method = "GN")

    obj_def <- r_def$ensemble_fit$gmm_fit$opt$objective
    obj_rnd <- r_rnd$ensemble_fit$gmm_fit$opt$objective

    if (abs(r_def$mu_ipw - r_rnd$mu_ipw) > 0.01) {
      n_diff <- n_diff + 1
      winner <- ifelse(obj_def <= obj_rnd, "default", "random")
      cat(sprintf("Rep %2d: DIFFERENT! def=%.4f(obj=%.2e), rnd=%.4f(obj=%.2e) [%s wins]\n",
                  i, r_def$mu_ipw, obj_def, r_rnd$mu_ipw, obj_rnd, winner))
      cat(sprintf("  nu_def=%s\n", paste(round(r_def$nu.hat, 4), collapse=", ")))
      cat(sprintf("  nu_rnd=%s\n", paste(round(r_rnd$nu.hat, 4), collapse=", ")))
    } else {
      cat(sprintf("Rep %2d: SAME mu=%.4f (obj_def=%.2e, obj_rnd=%.2e)\n",
                  i, r_def$mu_ipw, obj_def, obj_rnd))
    }
  }, error = function(e) {
    cat(sprintf("Rep %2d: ERROR - %s\n", i, conditionMessage(e)))
  })
}

cat(sprintf("\n%d/%d replicates had different results.\n", n_diff, n_reps))
cat("\nDone!\n")
