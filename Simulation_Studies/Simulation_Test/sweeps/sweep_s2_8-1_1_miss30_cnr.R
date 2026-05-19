## Sweep cnr cond for M1 alpha fit, scenario 8-1 setting2 miss30 _1 n=2000.
## J=1 (single model = full = r ~ y + u1 + u2). Per cnr cond workflow.
##
## Motivation: max|alpha| > 10 in ~3% of reps drives CP from 0.928 -> 0.958 if dropped.
## Likely binary-outcome separation problem at n=2000 with miss30 (~600 observed r=1).
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({ source("Basic_setup.r"); source("Data_Generation.r"); source("Simulation.r") })
library(EBMRalgorithmFast4); library(parallel); library(foreach); library(doSNOW)

W_sm <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_full <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2)
inv_link_fn <- function(eta) 1/(1+exp(eta))

make_ps_spec_M1 <- function(opt, cond) {
  list(formula.list = list(r ~ y + u1 + u2),
       h_alpha.list = list(h_full),
       inv_link = inv_link_fn, outcome = "y",
       optimizer = list(opt), cond_threshold = list(cond))
}

set.seed(12345); d <- setting2.A1(1e7); mu_true <- mean(d$y)
cat(sprintf("mu_true(setting2.A1, n=1e7) = %.4f\n\n", mu_true))

all_data <- readRDS("Simulation_Data/setting2.A1_n2000_replicate1000.RDS")
alpha.true <- correct_model_alpha.true.list$setting2$miss30[[1]]
cat(sprintf("alpha.true (setting2 miss30 correct model): %s\n\n", paste(round(alpha.true, 4), collapse=", ")))
ps_model.true <- function(dat) {
  X <- cbind(rep(1, nrow(dat)), dat$y, dat$u1, dat$u2)
  1 / (1 + exp(X %*% alpha.true))
}

nn <- 2000; n_reps <- 500
n_cores <- min(detectCores() - 2, 10)

run_one_cond <- function(opt, cond, tag) {
  cl <- makeCluster(n_cores); registerDoSNOW(cl)
  clusterExport(cl, c("nn","W_sm","h_full","inv_link_fn","all_data",
                      "alpha.true","ps_model.true","make_ps_spec_M1","opt","cond"),
                envir = environment())
  clusterEvalQ(cl, { library(EBMRalgorithmFast4) })
  pb <- txtProgressBar(max = n_reps, style = 3)
  opts <- list(progress = function(n) setTxtProgressBar(pb, n))
  res <- foreach(i = 1:n_reps, .combine = 'cbind', .options.snow = opts,
                 .packages = c("stringr","Matrix")) %dopar% {
    tryCatch({
      dat <- all_data[((i-1)*nn+1):(i*nn), ]
      sp <- make_ps_spec_M1(opt, cond)
      ebmr <- EBMRAlgorithmFast4$new("y", sp, dat, W_sm)
      out <- ebmr$EBMR_IPW(h_nu = h_full, true_ps = ps_model.true(dat),
                           se.fit = TRUE, type = "HT")
      a <- unname(ebmr$ps_fit.list[[1]]$coefficients)
      c(mu_ipw   = unname(out$mu_ipw),
        se_ipw   = unname(out$se_ipw),
        max_a    = max(abs(a)),
        a1=a[1], a2=a[2], a3=a[3], a4=a[4])
    }, error = function(e) c(mu_ipw=NA_real_, se_ipw=NA_real_, max_a=NA_real_,
                             a1=NA_real_, a2=NA_real_, a3=NA_real_, a4=NA_real_))
  }
  close(pb); stopCluster(cl)
  res
}

configs <- list(
  list(opt="L-BFGS-B",       cond=1e8, tag="default L-BFGS-B/1e8 (current)"),
  list(opt="constrained_nr", cond=1e2, tag="cnr/1e2"),
  list(opt="constrained_nr", cond=2e2, tag="cnr/2e2"),
  list(opt="constrained_nr", cond=5e2, tag="cnr/5e2"),
  list(opt="constrained_nr", cond=1e3, tag="cnr/1e3"),
  list(opt="constrained_nr", cond=1e4, tag="cnr/1e4"),
  list(opt="constrained_nr", cond=1e5, tag="cnr/1e5")
)

cat(sprintf("\n========== SWEEP (%d reps per config) ==========\n", n_reps))
cat(sprintf("%-32s | %s\n", "config",
            "Bias        ESD     ESE     ratio  CP    | max|a| med/q90/q99/max | n_blow"))
cat(strrep("-", 118), "\n", sep="")
for (cfg in configs) {
  t0 <- proc.time()
  r <- run_one_cond(cfg$opt, cfg$cond, cfg$tag)
  el <- (proc.time() - t0)["elapsed"]
  mu_v <- r["mu_ipw", ]; se_v <- r["se_ipw", ]; ma <- r["max_a", ]
  valid <- !is.na(mu_v) & !is.na(se_v) & is.finite(mu_v) & is.finite(se_v)
  mu_v <- mu_v[valid]; se_v <- se_v[valid]; ma_v <- ma[valid]
  bias <- mean(mu_v) - mu_true; esd <- sd(mu_v); ese <- mean(se_v)
  ci_lo <- mu_v - 1.96*se_v; ci_hi <- mu_v + 1.96*se_v
  cp <- mean(ci_lo <= mu_true & mu_true <= ci_hi)
  n_blow <- sum(ma_v > 10, na.rm = TRUE)
  cat(sprintf("%-32s | %+.4f  %.4f  %.4f  %.3f  %.3f | %.2f/%.2f/%.2f/%.2f      | %3d (%.0fs)\n",
              cfg$tag, bias, esd, ese, ese/esd, cp,
              median(ma_v, na.rm=T), quantile(ma_v, 0.9, na.rm=T),
              quantile(ma_v, 0.99, na.rm=T), max(ma_v, na.rm=T), n_blow, el))
}
cat("\nDONE\n")
