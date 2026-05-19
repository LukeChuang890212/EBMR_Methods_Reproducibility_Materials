## End-to-end sanity check: package supports nu_optimizer/nu_cond_threshold and
## scenarios.R correctly threads it from optimizer_overrides.txt to simulate().
## Runs 3 reps of s2.B1 9-3 _23 manually with both default and cnr/1e2 nu fits.
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({ source("Basic_setup.r"); source("Data_Generation.r"); source("Simulation.r") })
library(EBMRalgorithmFast4)

W_sm   <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_full <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2)
h_nu_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1_u2=dat$u1*dat$u2)
inv_link_fn <- function(eta) 1/(1+exp(eta))
ps_spec_23 <- list(
  formula.list = list(r ~ y + u1 + z2, r ~ y + u2 + z2),
  h_alpha.list = list(h_full, h_full),
  inv_link = inv_link_fn, outcome = "y",
  optimizer = list("L-BFGS-B", "constrained_nr"),
  cond_threshold = list(1e8, 1e4)
)
alpha.true <- misspecified_model_alpha.true.list$setting2$miss50[[1]]
ps_true_fn <- function(dat) {
  X <- cbind(rep(1, nrow(dat)), dat$y, dat$u1, dat$u2)
  1 / (1 + exp(X %*% alpha.true))
}
all_data <- readRDS("Simulation_Data/setting2.B1_n2000_replicate1000.RDS")
nn <- 2000

cat("Default nu (L-BFGS-B/1e8) vs cnr/1e2 across 3 reps:\n")
cat(sprintf("%4s | %-25s | %-25s\n", "rep", "default", "cnr/1e2"))
for (i in 1:3) {
  dat <- all_data[((i-1)*nn+1):(i*nn), ]
  e1 <- EBMRAlgorithmFast4$new("y", ps_spec_23, dat, W_sm)
  r1 <- e1$EBMR_IPW(h_nu_fn, true_ps = ps_true_fn(dat), se.fit = TRUE)
  e2 <- EBMRAlgorithmFast4$new("y", ps_spec_23, dat, W_sm)
  r2 <- e2$EBMR_IPW(h_nu_fn, true_ps = ps_true_fn(dat), se.fit = TRUE,
                    nu_optimizer = "constrained_nr", nu_cond_threshold = 1e2)
  cat(sprintf("%4d | mu=%.4f se=%.4f w1=%.3f | mu=%.4f se=%.4f w1=%.3f\n",
              i, r1$mu_ipw, r1$se_ipw, r1$w.hat[1],
                 r2$mu_ipw, r2$se_ipw, r2$w.hat[1]))
}
cat("\nDONE\n")
