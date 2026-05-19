## Head-to-head: does cnr/1e8 fix CP for the n=500 M1 fits?
## Tests both miss50 (A1 data) and miss30 (A2 data) for scenario 8-1 setting2.
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

# mu_true (same for A1/A2 since m(x) determines E[y])
set.seed(12345); mu_true <- mean(setting2.A1(1e7)$y)
cat(sprintf("mu_true(setting2) = %.4f\n\n", mu_true))

nn <- 500; n_reps <- 1000
n_cores <- min(detectCores() - 2, 10)

run_one <- function(data_file, alpha.true, opt, cond) {
  all_data <- readRDS(data_file)
  ps_model.true <- function(dat) {
    X <- cbind(rep(1, nrow(dat)), dat$y, dat$u1, dat$u2)
    1 / (1 + exp(X %*% alpha.true))
  }
  t0 <- proc.time()
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
      c(mu_ipw=unname(out$mu_ipw), se_ipw=unname(out$se_ipw), max_a=max(abs(a)))
    }, error = function(e) c(mu_ipw=NA_real_, se_ipw=NA_real_, max_a=NA_real_))
  }
  close(pb); stopCluster(cl)
  list(res = res, elapsed = (proc.time()-t0)["elapsed"])
}

report <- function(tag, out) {
  r <- out$res
  mu_v <- r["mu_ipw", ]; se_v <- r["se_ipw", ]; ma <- r["max_a", ]
  valid <- !is.na(mu_v) & !is.na(se_v) & is.finite(mu_v) & is.finite(se_v)
  mu_v <- mu_v[valid]; se_v <- se_v[valid]; ma_v <- ma[valid]
  bias <- mean(mu_v)-mu_true; esd <- sd(mu_v); ese <- mean(se_v)
  ci_lo <- mu_v - 1.96*se_v; ci_hi <- mu_v + 1.96*se_v
  cp <- mean(ci_lo <= mu_true & mu_true <= ci_hi)
  n_blow <- sum(ma_v > 10, na.rm = TRUE)
  zsc <- (mu_v - mu_true) / se_v
  cat(sprintf("%-22s | Bias=%+.4f ESD=%.4f ESE=%.4f ratio=%.3f CP=%.3f | max|a| q99=%.2f max=%.2f | blow=%d | z@0.975=%+.2f (%.0fs)\n",
              tag, bias, esd, ese, ese/esd, cp,
              quantile(ma_v, 0.99, na.rm=T), max(ma_v, na.rm=T), n_blow,
              quantile(zsc, 0.975), out$elapsed))
}

# ============= setting2 miss50 (A1 data) =============
cat("===== setting2 miss50 n=500 (A1 data) =====\n")
df50 <- "Simulation_Data/setting2.A1_n2000_replicate1000.RDS"
at50 <- correct_model_alpha.true.list$setting2$miss50[[1]]
cat(sprintf("alpha.true (miss50): %s\n", paste(round(at50, 4), collapse=", ")))
o50d <- run_one(df50, at50, "L-BFGS-B", 1e8);       report("default L-BFGS-B/1e8", o50d)
o50c <- run_one(df50, at50, "constrained_nr", 1e8); report("cnr/1e8",              o50c)

# ============= setting2 miss30 (A2 data) =============
cat("\n===== setting2 miss30 n=500 (A2 data) =====\n")
df30 <- "Simulation_Data/setting2.A2_n2000_replicate1000.RDS"
at30 <- correct_model_alpha.true.list$setting2$miss30[[1]]
cat(sprintf("alpha.true (miss30): %s\n", paste(round(at30, 4), collapse=", ")))
o30d <- run_one(df30, at30, "L-BFGS-B", 1e8);       report("default L-BFGS-B/1e8", o30d)
o30c <- run_one(df30, at30, "constrained_nr", 1e8); report("cnr/1e8",              o30c)

cat("\nDONE\n")
