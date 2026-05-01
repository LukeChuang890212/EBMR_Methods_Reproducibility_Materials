setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
source("Simulation.r")
library(EBMRalgorithmFast4)

ps_spec <- get_ps_spec("9-alt1")
ps_spec_2 <- list(
  formula.list = ps_spec[["formula.list"]][2],
  h_alpha.list = ps_spec[["h_alpha.list"]][2],
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]]
)
W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

# Use a representative "bounded at 5" rep
set.seed(12346)
dat <- setting4.B1(2000)

cat("=== CUE obj and hessian_pd at different bound levels ===\n")
cat(sprintf("%-8s  %12s  %10s  %10s  %8s\n",
    "bounds", "alpha_norm", "obj", "hess_pd", "pi_bnd"))

for (B in c(Inf, 50, 20, 15, 10, 7, 5, 3, 2, 1)) {
  lo <- if (is.infinite(B)) -Inf else -B
  hi <- if (is.infinite(B)) Inf  else  B

  ebmr <- EBMRAlgorithmFast4$new("y", ps_spec_2, dat, W_func,
                                  lower = lo, upper = hi)
  invisible(ebmr$EBMR_IPW(
    h_nu = function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1_u2=dat$u1*dat$u2),
    type="HT", se.fit=FALSE
  ))
  ps_fit <- ebmr$ps_fit.list[[1]]
  alpha  <- ps_fit$coefficients
  gfit   <- ps_fit$gmm_fit$opt
  pi_hat <- ps_fit$fitted.values
  pi_bnd <- sum(pi_hat < 0.05 | pi_hat > 0.95)
  cat(sprintf("%-8s  %12.4f  %10.6f  %8s  %8d\n",
      if (is.infinite(B)) "none" else paste0("±", B),
      sqrt(sum(alpha^2)), gfit$objective,
      as.character(gfit$hessian_pd), pi_bnd))
}
