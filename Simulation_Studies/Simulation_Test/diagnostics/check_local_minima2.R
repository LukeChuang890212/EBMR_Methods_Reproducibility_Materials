setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Data_Generation.r")
source("config/scenarios.R")
library(EPS)
library(numDeriv)
source("Basic_setup.r")

data_file <- misspecified_model_all_data_file.list$setting4$miss50[[1]]
all_data <- readRDS(data_file)
n_val <- 2000

ps_spec <- get_ps_spec("9-alt1")
subset_ps_spec <- list(
  formula.list = ps_spec$formula.list[c(1, 3)],
  h_alpha.list = ps_spec$h_alpha.list[c(1, 3)],
  inv_link = ps_spec$inv_link,
  outcome = ps_spec$outcome
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_nu_func <- function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2, u1_u2 = dat$u1*dat$u2)

# Test outlier rep 27 (stored w1=0.0008)
rep_i <- 27
dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]

cat("=== Rep 27 [OUTLIER] - Package result ===\n")
ebmr <- EPS$new("y", subset_ps_spec, dat, W_func)

# Call EBMR_IPW to trigger ensemble
result <- ebmr$EBMR_IPW(h_nu = h_nu_func, se.fit = FALSE)
cat(sprintf("Package: nu=[%.4f,%.4f], w=[%.4f,%.4f], mu=%.4f\n",
            result$nu.hat[1], result$nu.hat[2],
            result$w.hat[1], result$w.hat[2],
            result$mu_ipw))

# Now check: does the package's ensemble_fit have convergence info?
# Access the ensemble fit from the private method
# We can't easily access private, but let's check what EBMR_IPW returns
cat("\nResult fields:", paste(names(result), collapse=", "), "\n")

# Let's also check the alpha estimates
for (j in 1:2) {
  fit <- ebmr$ps_fit.list[[j]]
  cat(sprintf("Model %d: alpha=[%s], iter=%d, conv=%s, grad=%.2e\n",
              j, paste(round(fit$coefficients, 4), collapse=","),
              fit$gmm_fit$opt$iterations,
              fit$gmm_fit$opt$converged,
              fit$gmm_fit$opt$final_grad_norm))
}

# Now test: what happens if we disable analytical gradients?
# We can't easily do that, but let's compare the Phi_nu objective
# at the package's solution vs the manual solution
ps.matrix <- do.call(cbind, lapply(ebmr$ps_fit.list, function(f) f$fitted.values))
h_x2 <- h_nu_func(dat)
d <- ps.matrix[, -1, drop = FALSE] - ps.matrix[, 1]
h_x <- cbind(d, h_x2)
h_dim <- ncol(h_x)
r_vec <- as.numeric(dat$r)
n <- nrow(dat)

Phi_nu <- function(nu) {
  ps_nu <- as.vector(ps.matrix %*% nu)
  (r_vec / ps_nu - 1) * h_x
}

gmm_obj <- function(nu) {
  g_mat <- Phi_nu(nu)
  G_val <- colMeans(g_mat)
  W_mat <- tryCatch(solve(t(g_mat) %*% g_mat / n), error = function(e) diag(h_dim))
  as.numeric(t(G_val) %*% W_mat %*% G_val)
}

# Package solution
nu_pkg <- result$nu.hat
cat(sprintf("\nObjective at package solution nu=[%.4f,%.4f]: %.8f\n",
            nu_pkg[1], nu_pkg[2], gmm_obj(nu_pkg)))

# Manual "correct" solution (from the local minima check)
nu_manual <- c(1.0626, -0.0494)
cat(sprintf("Objective at manual solution  nu=[%.4f,%.4f]: %.8f\n",
            nu_manual[1], nu_manual[2], gmm_obj(nu_manual)))

# Also check: what nu does the stored result correspond to?
nu_stored <- res <- readRDS("Simulation_Results/EBMR_IPW_setting4-miss50-scenario9-2_13_n2000_replicate1000_test57.RDS")
nu1_s <- res["nu.hat1", rep_i]; nu2_s <- res["nu.hat2", rep_i]
cat(sprintf("Stored nu=[%.4f,%.4f], obj=%.8f\n", nu1_s, nu2_s, gmm_obj(c(nu1_s, nu2_s))))

cat("\nDone!\n")
