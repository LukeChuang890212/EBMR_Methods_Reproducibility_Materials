## Sweep cond for s2.B1 n=2000, M3 only (matches test65 file
## EBMR_IPW_setting2-miss50-scenario9-3_3_n2000_replicate1000_test65.RDS).
## M3: formula r ~ y + u2 + z2, h_alpha = cbind(u1,u2,z1,z2).
## Compare default L-BFGS-B/1e8 vs constrained_nr/{1e2, 5e2, 1e3, 1e4, 1e5}.
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({ source("Basic_setup.r"); source("Data_Generation.r"); source("Simulation.r") })
devtools::load_all("../EIPS", quiet = TRUE)
library(parallel); library(foreach); library(doSNOW)

W_sm <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_full <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2)
h_nu_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1_u2=dat$u1*dat$u2)
inv_link_fn <- function(eta) 1/(1+exp(eta))

make_spec <- function(optimizer, cond) list(
  formula.list = list(r ~ y + u2 + z2),
  h_alpha.list = list(h_full),
  inv_link = inv_link_fn, outcome = "y",
  optimizer = optimizer, cond_threshold = cond
)

n_cores <- min(detectCores() - 1, 10)
n_reps  <- 500
cfg <- list(label = "s2.B1 n=2000 (9-3 M3)",
            dgp = "setting2.B1", n = 2000,
            mu_true = get_mu_true("setting2"))

configs_opt <- list(
  list(opt = "L-BFGS-B",       cond = 1e8, tag = "L-BFGS-B/1e8 (default)"),
  list(opt = "constrained_nr", cond = 1e2, tag = "cnr/1e2"),
  list(opt = "constrained_nr", cond = 5e2, tag = "cnr/500"),
  list(opt = "constrained_nr", cond = 1e3, tag = "cnr/1e3"),
  list(opt = "constrained_nr", cond = 1e4, tag = "cnr/1e4"),
  list(opt = "constrained_nr", cond = 1e5, tag = "cnr/1e5")
)

run_one <- function(cfg, optspec) {
  ps_spec <- make_spec(optspec$opt, optspec$cond)
  t0 <- proc.time()
  cl <- makeCluster(n_cores); registerDoSNOW(cl)
  clusterExport(cl, c("cfg","ps_spec","W_sm","h_full","h_nu_fn","inv_link_fn"),
                envir = environment())
  clusterEvalQ(cl, {
    setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
    devtools::load_all("../EIPS", quiet = TRUE)
    source("Data_Generation.r")
  })
  pb <- txtProgressBar(max = n_reps, style = 3)
  opts <- list(progress = function(n) setTxtProgressBar(pb, n))
  res <- foreach(i = 1:n_reps, .combine = 'rbind', .options.snow = opts,
                 .packages = c("stringr","Matrix")) %dopar% {
    out <- tryCatch({
      set.seed(i)
      dat <- get(cfg$dgp)(cfg$n)
      ebmr <- EIPS$new("y", ps_spec, dat, W_sm)
      ipw <- ebmr$EBMR_IPW(h_nu_fn, model_indices = 1, se.fit = TRUE)
      a <- unname(ebmr$ps_fit.list[[1]]$coefficients)
      c(mu_ipw = ipw$mu_ipw, se_ipw = ipw$se_ipw,
        a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4],
        max_abs_alpha = max(abs(a)))
    }, error = function(e) c(mu_ipw = NA, se_ipw = NA,
                             a1 = NA, a2 = NA, a3 = NA, a4 = NA, max_abs_alpha = NA))
    out
  }
  close(pb); stopCluster(cl)
  list(res = res, elapsed = (proc.time() - t0)["elapsed"])
}

cat(sprintf("\n========= SWEEP: %s (mu_true=%.4f) =========\n", cfg$label, cfg$mu_true))
out_all <- list()
for (optspec in configs_opt) {
  out <- run_one(cfg, optspec)
  res <- out$res
  # Use project's clean_sim_result()
  sim_mat <- rbind(
    mu_ipw      = res[, "mu_ipw"],
    mu_ipw.true = res[, "mu_ipw"],
    se_ipw      = res[, "se_ipw"],
    se_ipw.true = res[, "se_ipw"],
    alpha.hat1  = res[, "a1"],
    alpha.hat2  = res[, "a2"],
    alpha.hat3  = res[, "a3"],
    alpha.hat4  = res[, "a4"]
  )
  cleaned <- clean_sim_result(sim_mat, multiplier = 2, max_pct = 0.01, verbose = FALSE)
  sc <- cleaned$result
  mu_v <- sc["mu_ipw", ]; se_v <- sc["se_ipw", ]
  bias <- mean(mu_v) - cfg$mu_true; esd <- sd(mu_v); ese <- mean(se_v)
  ratio <- ese / esd
  ci_lo <- mu_v - 1.96*se_v; ci_hi <- mu_v + 1.96*se_v
  cp <- mean(ci_lo <= cfg$mu_true & cfg$mu_true <= ci_hi)
  alpha_mat <- res[, c("a1","a2","a3","a4")]
  n_zero <- sum(rowSums(abs(alpha_mat), na.rm = TRUE) == 0, na.rm = TRUE)
  n_sep  <- sum(res[, "max_abs_alpha"] > 50, na.rm = TRUE)
  sd_a   <- apply(alpha_mat, 2, sd, na.rm = TRUE)
  cat(sprintf("  %-22s | stuck@0=%4d/%d sep(>50)=%3d NA=%d Out=%d | Bias=%+.4f ESD=%.4f ESE=%.4f ESE/ESD=%.3f CP=%.3f | sd(a)=(%.2f,%.2f,%.2f,%.2f) max|a|=%.1f (%.0fs)\n",
              optspec$tag, n_zero, n_reps, n_sep,
              cleaned$n_na, cleaned$n_outliers,
              bias, esd, ese, ratio, cp,
              sd_a[1], sd_a[2], sd_a[3], sd_a[4],
              max(res[,"max_abs_alpha"], na.rm = TRUE), out$elapsed))
  out_all[[optspec$tag]] <- list(res = res)
}
saveRDS(out_all, "Simulation_Test/sweep_cnr_s2_miss50_n2000_M3_results.RDS")
cat("\n========= DONE =========\n")
