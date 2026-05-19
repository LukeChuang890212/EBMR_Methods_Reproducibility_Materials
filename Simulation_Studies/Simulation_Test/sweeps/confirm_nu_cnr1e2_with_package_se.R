## Confirm cnr/1e2 nu fix for s2.B1 9-3 _23 using PACKAGE'S sandwich SE (R + S_mat).
## Approach: monkey-patch ebmr$.__enclos_env__$private$gmm AFTER $new() so the nu fit
## call inside EBMR_IPW (Methods.r:496) uses constrained_nr/1e2 by default, while
## alpha fits (done inside $new()) still use ps_spec_23's per-model override.
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({ source("Basic_setup.r"); source("Data_Generation.r"); source("Simulation.r") })
devtools::load_all("../EPS", quiet = TRUE)
library(parallel); library(foreach); library(doSNOW)

W_sm    <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_full  <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2)
h_nu_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1_u2=dat$u1*dat$u2)
inv_link_fn <- function(eta) 1/(1+exp(eta))

ps_spec_23 <- list(
  formula.list   = list(r ~ y + u1 + z2, r ~ y + u2 + z2),
  h_alpha.list   = list(h_full, h_full),
  inv_link       = inv_link_fn, outcome = "y",
  optimizer      = list("L-BFGS-B",       "constrained_nr"),
  cond_threshold = list(1e8,              1e4)
)

alpha.true <- misspecified_model_alpha.true.list$setting2$miss50[[1]]
ps_model.true <- function(dat, alpha.true) {
  X <- cbind(rep(1, nrow(dat)), dat$y, dat$u1, dat$u2)
  1 / (1 + exp(X %*% alpha.true))
}
all_data <- readRDS("Simulation_Data/setting2.B1_n2000_replicate1000.RDS")

n_cores <- min(detectCores() - 1, 10); n_reps <- 1000; nn <- 2000
mu_true <- mean(all_data$y)
cat(sprintf("\nConfig: s2.B1 _23, nu cnr/1e2 via monkey-patch, %d reps, mu_true=%.4f\n",
            n_reps, mu_true))

t0 <- proc.time()
cl <- makeCluster(n_cores); registerDoSNOW(cl)
clusterExport(cl, c("nn","ps_spec_23","W_sm","h_full","h_nu_fn","inv_link_fn",
                    "all_data","alpha.true","ps_model.true"),
              envir = environment())
clusterEvalQ(cl, {
  setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
  devtools::load_all("../EPS", quiet = TRUE)
  # Wrapper that defaults the nu-fit optimizer/cond to cnr/1e2.
  # Methods.r:496 calls private$gmm(...) WITHOUT optimizer/cond_threshold, so the
  # wrapper's defaults apply. For alpha fits this wrapper is NEVER called because
  # they were already finished inside $new() before the patch is installed.
  my_gmm_nu_cnr1e2 <- function(g, W, n, esteq_dim, param_dim, init, se.fit = TRUE,
                               dg = NULL, d2g = NULL, lower = -Inf, upper = Inf,
                               stall_limit = 20L, max_outer_override = NULL,
                               Gamma_direct = NULL, optimizer = "constrained_nr",
                               cond_threshold = 1e2, trust_radius = 2.0) {
    EPS:::gmm(g, W, n, esteq_dim, param_dim, init, se.fit, dg, d2g,
                             lower, upper, stall_limit, max_outer_override,
                             Gamma_direct, optimizer, cond_threshold, trust_radius)
  }
})
pb <- txtProgressBar(max = n_reps, style = 3)
opts <- list(progress = function(n) setTxtProgressBar(pb, n))
res <- foreach(i = 1:n_reps, .combine = 'cbind', .options.snow = opts,
               .packages = c("stringr","Matrix")) %dopar% {
  tryCatch({
    dat <- all_data[((i - 1) * nn + 1):(i * nn), ]
    ebmr <- EPS$new("y", ps_spec_23, dat, W_sm)
    # Monkey-patch the private gmm AFTER alpha fits complete (R6 locks bindings,
    # so go via unlockBinding). Only the nu fit inside EBMR_IPW (Methods.r:496)
    # calls private$gmm without optimizer/cond_threshold, so cnr/1e2 defaults apply.
    pe <- ebmr$.__enclos_env__$private
    unlockBinding("gmm", pe); assign("gmm", my_gmm_nu_cnr1e2, envir = pe); lockBinding("gmm", pe)
    out <- ebmr$EBMR_IPW(h_nu_fn, true_ps = ps_model.true(dat, alpha.true),
                         se.fit = TRUE, type = "HT")
    nu_hat <- as.vector(out$nu.hat); w_hat <- as.vector(out$w.hat)
    # Canonical row order expected by clean_sim_result: rows 1-4 = mu_ipw, mu_ipw.true, se_ipw, se_ipw.true
    c(mu_ipw      = unname(out$mu_ipw),
      mu_ipw.true = unname(out$mu_ipw.true),
      se_ipw      = unname(out$se_ipw),
      se_ipw.true = unname(out$se_ipw.true),
      nu1 = nu_hat[1], nu2 = nu_hat[2],
      w1  = w_hat[1],  w2  = w_hat[2])
  }, error = function(e) c(mu_ipw=NA_real_, mu_ipw.true=NA_real_,
                           se_ipw=NA_real_, se_ipw.true=NA_real_,
                           nu1=NA_real_, nu2=NA_real_, w1=NA_real_, w2=NA_real_))
}
close(pb); stopCluster(cl)
elapsed <- (proc.time() - t0)["elapsed"]

saveRDS(res, "Simulation_Test/confirm_nu_cnr1e2_with_package_se_res.RDS")

report <- function(label, sc) {
  mu_v <- sc["mu_ipw", ]; se_v <- sc["se_ipw", ]; w1_v <- sc["w1", ]
  bias <- mean(mu_v) - mu_true; esd <- sd(mu_v)
  ese_mean <- mean(se_v); ese_med <- median(se_v); ratio <- ese_mean / esd
  ci_lo <- mu_v - 1.96*se_v; ci_hi <- mu_v + 1.96*se_v
  cp <- mean(ci_lo <= mu_true & mu_true <= ci_hi)
  n_M3 <- sum(w1_v < 0.1, na.rm = TRUE)
  n_mid <- sum(w1_v >= 0.1 & w1_v < 0.9, na.rm = TRUE)
  n_M2 <- sum(w1_v >= 0.9, na.rm = TRUE)
  cat(sprintf("\n=== %s ===\n", label))
  cat(sprintf("n=%d  Bias=%+.4f  ESD=%.4f  ESE(mean)=%.4f  ESE(median)=%.4f  ratio=%.3f  CP=%.3f\n",
              ncol(sc), bias, esd, ese_mean, ese_med, ratio, cp))
  cat(sprintf("w_M2 dist: [0,0.1):%d  [0.1,0.9):%d  [0.9,1]:%d  mean=%.3f  sd=%.3f\n",
              n_M3, n_mid, n_M2, mean(w1_v, na.rm=T), sd(w1_v, na.rm=T)))
}

# Raw (NA-only filtered, no outlier removal)
raw_mask <- apply(res[1:4, , drop = FALSE], 2, function(x) !any(is.na(x) | !is.finite(x)))
report("RAW (NA/Inf only)", res[, raw_mask, drop = FALSE])

# clean_sim_result with default multiplier=2, max_pct=0.01
cleaned <- clean_sim_result(res, multiplier = 2, max_pct = 0.01, verbose = TRUE)
cat(sprintf("\nclean_sim_result: n_total=%d  n_na=%d  n_outliers=%d  n_successful=%d\n",
            cleaned$n_total, cleaned$n_na, cleaned$n_outliers, cleaned$n_successful))
report("CLEANED (mult=2, max_pct=0.01)", cleaned$result)

cat(sprintf("\nelapsed: %.0fs\nDONE\n", elapsed))
