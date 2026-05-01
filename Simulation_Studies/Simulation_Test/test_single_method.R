setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
source("Data_Generation.r")
source("config/scenarios.R")

old_wd <- getwd()
setwd("../EBMRalgorithmFast4")
source("R/EBMRAlgorithm.r")
setwd(old_wd)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
inv_link_fn <- function(eta) 1 / (1 + exp(eta))

ps_spec <- list(
  formula.list = list(FORMULAS$u1_z2),
  h_alpha.list = list(H_ALPHA$full),
  inv_link = inv_link_fn,
  outcome = "y",
  alpha_init.list = list(NULL)
)

h_nu_func <- function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1,
                                  z2 = dat$z2, u1_u2 = dat$u1 * dat$u2)

# Quick test: one dataset
set.seed(42)
dat <- setting3.B1(2000)
ebmr <- EBMRAlgorithmFast4$new("y", ps_spec, dat, W_func, bounds = Inf)
ipw <- ebmr$EBMR_IPW(h_nu = h_nu_func)
fit <- ebmr$ps_fit.list[[1]]

cat(sprintf("mu_ipw = %.4f, se_ipw = %.4f\n", ipw$mu_ipw, ipw$se_ipw))
cat(sprintf("alpha = [%s]\n", paste(round(fit$coefficients, 4), collapse=", ")))
cat(sprintf("converged = %s, iterations = %d, grad_norm = %.2e\n",
            fit$gmm_fit$opt$converged, fit$gmm_fit$opt$iterations, fit$gmm_fit$opt$final_grad_norm))
cat("OK\n")
