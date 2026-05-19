## Sweep cnr cond for M2 alpha fit, scenario 9-1 setting1 miss50 _2 n=2000.
## M2 = u1_z2 + h_alpha=full (misspecified relative to truth r ~ y+u1+u2).
## Current default: 1 rep blows up to alpha=503.7, SE=106 -> ratio=1.91 after cleaning.
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({ source("Basic_setup.r"); source("Data_Generation.r"); source("Simulation.r") })
library(EBMRalgorithmFast4); library(parallel); library(foreach); library(doSNOW)

W_sm <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_full <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2)
inv_link_fn <- function(eta) 1/(1+exp(eta))

make_ps_spec_M2 <- function(opt, cond) {
  list(formula.list = list(r ~ y + u1 + z2),       # 9-1's M2 = u1_z2
       h_alpha.list = list(h_full),
       inv_link = inv_link_fn, outcome = "y",
       optimizer = list(opt), cond_threshold = list(cond))
}

set.seed(12345); mu_true <- mean(setting1.A1(1e7)$y)
cat(sprintf("mu_true(setting1.A1) = %.4f\n\n", mu_true))

all_data <- readRDS("Simulation_Data/setting1.A1_n2000_replicate1000.RDS")
alpha.true <- correct_model_alpha.true.list$setting1$miss50[[1]]
cat(sprintf("alpha.true: %s\n\n", paste(round(alpha.true, 4), collapse=", ")))
ps_model.true <- function(dat) {
  X <- cbind(rep(1, nrow(dat)), dat$y, dat$u1, dat$u2)
  1 / (1 + exp(X %*% alpha.true))
}

nn <- 2000; n_reps <- 1000
n_cores <- min(detectCores() - 2, 10)

run_one <- function(opt, cond) {
  t0 <- proc.time()
  cl <- makeCluster(n_cores); registerDoSNOW(cl)
  clusterExport(cl, c("nn","W_sm","h_full","inv_link_fn","all_data",
                      "alpha.true","ps_model.true","make_ps_spec_M2","opt","cond"),
                envir = environment())
  clusterEvalQ(cl, { library(EBMRalgorithmFast4) })
  pb <- txtProgressBar(max = n_reps, style = 3)
  opts <- list(progress = function(n) setTxtProgressBar(pb, n))
  res <- foreach(i = 1:n_reps, .combine = 'cbind', .options.snow = opts,
                 .packages = c("stringr","Matrix")) %dopar% {
    tryCatch({
      dat <- all_data[((i-1)*nn+1):(i*nn), ]
      sp <- make_ps_spec_M2(opt, cond)
      ebmr <- EBMRAlgorithmFast4$new("y", sp, dat, W_sm)
      out <- ebmr$EBMR_IPW(h_nu = h_full, true_ps = ps_model.true(dat),
                           se.fit = TRUE, type = "HT")
      a <- unname(ebmr$ps_fit.list[[1]]$coefficients)
      c(mu_ipw=unname(out$mu_ipw), se_ipw=unname(out$se_ipw),
        mu_ipw.true=unname(out$mu_ipw.true), se_ipw.true=unname(out$se_ipw.true),
        max_a=max(abs(a)))
    }, error = function(e) c(mu_ipw=NA_real_, se_ipw=NA_real_,
                             mu_ipw.true=NA_real_, se_ipw.true=NA_real_, max_a=NA_real_))
  }
  close(pb); stopCluster(cl)
  list(res = res, elapsed = (proc.time()-t0)["elapsed"])
}

report <- function(tag, out) {
  r <- out$res
  # Add the standard 4 rows for clean_sim_result
  m <- matrix(NA, 4, ncol(r), dimnames=list(c("mu_ipw","mu_ipw.true","se_ipw","se_ipw.true"), NULL))
  m[1,] <- r["mu_ipw",]; m[2,] <- r["mu_ipw.true",]; m[3,] <- r["se_ipw",]; m[4,] <- r["se_ipw.true",]
  cl <- clean_sim_result(m, multiplier = 2, max_pct = 0.01, verbose = FALSE)
  xc <- cl$result
  mu_v <- xc[1,]; se_v <- xc[3,]; ma_full <- r["max_a", ]
  # restrict ma to cleaned set: not directly possible, just report overall ma stats
  ma_v <- ma_full[!is.na(ma_full)]
  bias <- mean(mu_v)-mu_true; esd <- sd(mu_v); ese <- mean(se_v); ese_med <- median(se_v)
  ci_lo <- mu_v - 1.96*se_v; ci_hi <- mu_v + 1.96*se_v
  cp <- mean(ci_lo <= mu_true & mu_true <= ci_hi)
  n_blow <- sum(ma_v > 10, na.rm = TRUE)
  cat(sprintf("%-26s | Bias=%+.4f ESD=%.4f ESE(m)=%.4f ESE(med)=%.4f ratio=%.3f CP=%.3f | max|a| q99=%.2f max=%.2f | blow=%d | clean removed=%d (%.0fs)\n",
              tag, bias, esd, ese, ese_med, ese/esd, cp,
              quantile(ma_v, 0.99, na.rm=T), max(ma_v, na.rm=T), n_blow,
              cl$n_outliers + cl$n_na, out$elapsed))
}

configs <- list(
  list(opt="L-BFGS-B",       cond=1e8, tag="default L-BFGS-B/1e8"),
  list(opt="constrained_nr", cond=1e2, tag="cnr/1e2"),
  list(opt="constrained_nr", cond=2e2, tag="cnr/2e2"),
  list(opt="constrained_nr", cond=5e2, tag="cnr/5e2"),
  list(opt="constrained_nr", cond=1e3, tag="cnr/1e3"),
  list(opt="constrained_nr", cond=1e4, tag="cnr/1e4"),
  list(opt="constrained_nr", cond=1e5, tag="cnr/1e5"),
  list(opt="constrained_nr", cond=1e6, tag="cnr/1e6"),
  list(opt="constrained_nr", cond=1e8, tag="cnr/1e8")
)

cat("===== Sweep (cleaned via clean_sim_result mult=2, max_pct=0.01) =====\n")
for (cfg in configs) {
  o <- run_one(cfg$opt, cfg$cond)
  report(cfg$tag, o)
}
cat("\nDONE\n")
