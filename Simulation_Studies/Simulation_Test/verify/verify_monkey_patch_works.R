## Quick sanity check: does monkey-patching ebmr$.__enclos_env__$private$gmm
## actually re-route the nu fit through my wrapper?
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({ source("Basic_setup.r"); source("Data_Generation.r"); source("Simulation.r") })
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

W_sm <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
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
ps_model.true <- function(dat, alpha.true) {
  X <- cbind(rep(1, nrow(dat)), dat$y, dat$u1, dat$u2)
  1 / (1 + exp(X %*% alpha.true))
}
all_data <- readRDS("Simulation_Data/setting2.B1_n2000_replicate1000.RDS")
nn <- 2000

call_log <- new.env()
call_log$nu_optimizer <- character(0)
call_log$nu_cond      <- numeric(0)

# Wrapper that RECORDS what optimizer/cond it was called with, then forwards.
spy_gmm <- function(g, W, n, esteq_dim, param_dim, init, se.fit = TRUE,
                    dg = NULL, d2g = NULL, lower = -Inf, upper = Inf,
                    stall_limit = 20L, max_outer_override = NULL,
                    Gamma_direct = NULL, optimizer = "constrained_nr",
                    cond_threshold = 1e2, trust_radius = 2.0) {
  call_log$nu_optimizer <- c(call_log$nu_optimizer, optimizer)
  call_log$nu_cond      <- c(call_log$nu_cond, cond_threshold)
  EBMRalgorithmFast4:::gmm(g, W, n, esteq_dim, param_dim, init, se.fit, dg, d2g,
                           lower, upper, stall_limit, max_outer_override,
                           Gamma_direct, optimizer, cond_threshold, trust_radius)
}

cat("===== Rep 1: WITHOUT patch (baseline package behavior) =====\n")
dat <- all_data[1:nn, ]
ebmr <- EBMRAlgorithmFast4$new("y", ps_spec_23, dat, W_sm)
out_un <- ebmr$EBMR_IPW(h_nu_fn, true_ps = ps_model.true(dat, alpha.true), se.fit = TRUE, type = "HT")
cat(sprintf("  mu_ipw=%.4f  se_ipw=%.4f  nu1=%.3f  nu2=%.3f  w1=%.3f\n",
            out_un$mu_ipw, out_un$se_ipw, out_un$nu.hat[1], out_un$nu.hat[2], out_un$w.hat[1]))

cat("\n===== Rep 1: WITH patch (cnr/1e2 for nu fit) =====\n")
ebmr2 <- EBMRAlgorithmFast4$new("y", ps_spec_23, dat, W_sm)
priv_env <- ebmr2$.__enclos_env__$private
patched <- tryCatch({
  unlockBinding("gmm", priv_env)
  assign("gmm", spy_gmm, envir = priv_env)
  lockBinding("gmm", priv_env)
  TRUE
}, error = function(e) { cat("PATCH FAILED:", conditionMessage(e), "\n"); FALSE })
if (patched) {
  out_p <- ebmr2$EBMR_IPW(h_nu_fn, true_ps = ps_model.true(dat, alpha.true), se.fit = TRUE, type = "HT")
  cat(sprintf("  mu_ipw=%.4f  se_ipw=%.4f  nu1=%.3f  nu2=%.3f  w1=%.3f\n",
              out_p$mu_ipw, out_p$se_ipw, out_p$nu.hat[1], out_p$nu.hat[2], out_p$w.hat[1]))
  cat(sprintf("  spy log: optimizer(s)=[%s]  cond(s)=[%s]\n",
              paste(call_log$nu_optimizer, collapse=", "),
              paste(call_log$nu_cond, collapse=", ")))
  cat(sprintf("  (Expected: 1 call, optimizer='constrained_nr', cond=1e2)\n"))
}

cat("\n===== Reps 2-3: confirm reproducibility =====\n")
for (i in 2:3) {
  dat_i <- all_data[((i-1)*nn+1):(i*nn), ]
  ebmr_i <- EBMRAlgorithmFast4$new("y", ps_spec_23, dat_i, W_sm)
  pe_i <- ebmr_i$.__enclos_env__$private
  unlockBinding("gmm", pe_i); assign("gmm", spy_gmm, envir = pe_i); lockBinding("gmm", pe_i)
  out_i <- ebmr_i$EBMR_IPW(h_nu_fn, true_ps = ps_model.true(dat_i, alpha.true), se.fit = TRUE, type = "HT")
  cat(sprintf("  Rep %d: mu_ipw=%.4f  se_ipw=%.4f  w1=%.3f\n",
              i, out_i$mu_ipw, out_i$se_ipw, out_i$w.hat[1]))
}
cat(sprintf("\nTotal nu-gmm calls intercepted: %d (expected 3 across reps 1-3)\n",
            length(call_log$nu_optimizer)))
cat("DONE\n")
