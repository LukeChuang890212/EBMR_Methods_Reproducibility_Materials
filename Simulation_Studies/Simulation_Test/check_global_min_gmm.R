## Check whether GN finds the global minimum of the GMM objective
## by comparing with multi-start L-BFGS-B optimization.
##
## Key question: is the GMM problem just-identified or over-identified?
## Just-identified: Q(alpha*) = 0, unique global min.
## Over-identified: Q(alpha*) > 0, local minima possible.

setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
source("Data_Generation.r")
source("config/scenarios.R")

old_wd <- getwd()
setwd("../EBMRalgorithmFast4")
source("R/EBMRAlgorithm.r")
setwd(old_wd)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

# ---- Test 1: Check dimensions (over-identified?) ----
cat("=== Dimension Check ===\n\n")

ps_spec <- PS_SPECS[["9"]]
ps_spec$alpha_init.list <- list(NULL, NULL, NULL)
ps_spec$outcome <- "y"

set.seed(42)
dat <- setting3.B1(500)

ebmr <- EBMRAlgorithmFast4$new("y", ps_spec, dat, W_func, method = "GN")

for (j in 1:3) {
  fit <- ebmr$ps_fit.list[[j]]
  alpha_dim <- length(fit$coefficients)
  h_dim <- ncol(fit$h_x)
  cat(sprintf("Model %d: alpha_dim=%d (params), h_dim=%d (esteq)\n",
              j, alpha_dim, h_dim))
  if (h_dim > alpha_dim) {
    cat("  -> OVER-IDENTIFIED (h_dim > alpha_dim)\n")
  } else if (h_dim == alpha_dim) {
    cat("  -> JUST-IDENTIFIED (h_dim == alpha_dim)\n")
  } else {
    cat("  -> UNDER-IDENTIFIED (h_dim < alpha_dim)\n")
  }
  cat(sprintf("  GMM objective: %.6e\n", fit$opt$objective))
  cat(sprintf("  Converged: %s\n\n", fit$opt$converged))
}

# ---- Test 2: Multi-start comparison for PS estimation ----
cat("\n=== Multi-Start Comparison for PS Estimation ===\n\n")

# For each PS model, try multiple initial values and compare
set.seed(42)
dat <- setting4.B1(500)

for (j in 1:3) {
  cat(sprintf("--- Model %d ---\n", j))

  # GN with default init (zero)
  ebmr_default <- EBMRAlgorithmFast4$new("y", ps_spec, dat, W_func, method = "GN")
  fit_default <- ebmr_default$ps_fit.list[[j]]
  obj_default <- fit_default$opt$objective
  alpha_default <- fit_default$coefficients

  cat(sprintf("  Default init (0): obj=%.6e, alpha=%s\n",
              obj_default, paste(round(alpha_default, 4), collapse = ", ")))

  # Try several random initial values
  n_starts <- 10
  alpha_dim <- length(alpha_default)
  best_obj <- obj_default
  best_alpha <- alpha_default

  for (s in 1:n_starts) {
    init_s <- rnorm(alpha_dim, sd = 0.5)
    ps_spec_s <- ps_spec
    ps_spec_s$alpha_init.list[[j]] <- init_s

    tryCatch({
      ebmr_s <- EBMRAlgorithmFast4$new("y", ps_spec_s, dat, W_func, method = "GN")
      fit_s <- ebmr_s$ps_fit.list[[j]]
      obj_s <- fit_s$opt$objective

      if (obj_s < best_obj - 1e-10) {
        cat(sprintf("  Start %d: obj=%.6e ** BETTER than default **\n", s, obj_s))
        best_obj <- obj_s
        best_alpha <- fit_s$coefficients
      }
    }, error = function(e) {
      cat(sprintf("  Start %d: ERROR - %s\n", s, conditionMessage(e)))
    })
  }

  if (abs(best_obj - obj_default) < 1e-10) {
    cat(sprintf("  All %d starts converged to same objective. GN finds global min.\n\n", n_starts))
  } else {
    cat(sprintf("  DIFFERENT minima found! Default=%.6e, Best=%.6e\n", obj_default, best_obj))
    cat(sprintf("  Best alpha: %s\n\n", paste(round(best_alpha, 4), collapse = ", ")))
  }
}

# ---- Test 3: Multi-start for ensemble ----
cat("\n=== Multi-Start Comparison for Ensemble ===\n\n")

h_nu_func <- function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1,
                                  z2 = dat$z2, u1_u2 = dat$u1 * dat$u2)

# Run 20 replicates, check if different inits give different ensemble results
set.seed(123)
n_reps <- 20
n_diff <- 0

for (i in 1:n_reps) {
  dat <- setting4.B1(500)

  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", ps_spec, dat, W_func, method = "GN")

    # Default init
    result_default <- ebmr$EBMR_IPW(h_nu = h_nu_func, method = "GN")
    mu_default <- result_default$mu_ipw
    nu_default <- result_default$nu.hat

    # Random init
    J <- length(ebmr$ps_fit.list)
    result_random <- ebmr$EBMR_IPW(h_nu = h_nu_func,
                                    nu_init = rnorm(J, sd = 0.5),
                                    method = "GN")
    mu_random <- result_random$mu_ipw
    nu_random <- result_random$nu.hat

    if (abs(mu_default - mu_random) > 0.01) {
      n_diff <- n_diff + 1
      cat(sprintf("  Rep %d: DIFFERENT! default=%.4f, random=%.4f (diff=%.4f)\n",
                  i, mu_default, mu_random, mu_default - mu_random))
      cat(sprintf("    nu_default: %s\n", paste(round(nu_default, 4), collapse = ", ")))
      cat(sprintf("    nu_random:  %s\n", paste(round(nu_random, 4), collapse = ", ")))
    }
  }, error = function(e) {
    cat(sprintf("  Rep %d: ERROR - %s\n", i, conditionMessage(e)))
  })
}

cat(sprintf("\n  %d/%d replicates had different results with different inits.\n", n_diff, n_reps))
if (n_diff == 0) {
  cat("  GN ensemble is robust to initial values.\n")
} else {
  cat("  GN ensemble is SENSITIVE to initial values - multiple local minima exist!\n")
}

cat("\nDone!\n")
