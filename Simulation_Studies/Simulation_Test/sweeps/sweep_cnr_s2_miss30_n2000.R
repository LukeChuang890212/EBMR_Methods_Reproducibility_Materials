## Sweep cond for s2.B2 n=2000, M3 only (matches test61 file
## EBMR_IPW_setting2-miss30-scenario9-3_3_n2000_replicate1000_test61.RDS).
## Tests cond ∈ {500, 1e3, 1e4, 1e5} with 500 reps to find the sweet spot.
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

clean_iqr <- function(x, mult = 2, max_pct = 0.01) {
  Q1 <- quantile(x, 0.25); Q3 <- quantile(x, 0.75); iqr <- Q3 - Q1
  is_out <- x < (Q1 - mult*iqr) | x > (Q3 + mult*iqr)
  max_rm <- floor(max_pct * length(x))
  if (sum(is_out) > max_rm && max_rm > 0) {
    d <- abs(x - median(x)); oi <- which(is_out)
    keep <- oi[order(d[oi], decreasing = TRUE)[(max_rm+1):length(oi)]]
    is_out[keep] <- FALSE
  }
  !is_out
}

n_cores <- min(detectCores() - 1, 10)
n_reps  <- 500
cfg <- list(label = "s2.B2 n=2000 (miss30, 9-3 M3)",
            dgp = "setting2.B2", n = 2000,
            mu_true = get_mu_true("setting2"))

# Include default for baseline; then sweep constrained_nr conds
configs_opt <- list(
  list(opt = "constrained_nr", cond = 150, tag = "cnr/150"),
  list(opt = "constrained_nr", cond = 200, tag = "cnr/200"),
  list(opt = "constrained_nr", cond = 300, tag = "cnr/300"),
  list(opt = "constrained_nr", cond = 400, tag = "cnr/400")
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
  mu_v <- res[,"mu_ipw"]; se_v <- res[,"se_ipw"]
  valid <- !is.na(mu_v) & !is.na(se_v) & is.finite(mu_v) & is.finite(se_v)
  mu_v <- mu_v[valid]; se_v <- se_v[valid]
  keep <- clean_iqr(mu_v); mu_c <- mu_v[keep]; se_c <- se_v[keep]
  ci_lo <- mu_c - 1.96*se_c; ci_hi <- mu_c + 1.96*se_c
  cp <- mean(ci_lo <= cfg$mu_true & cfg$mu_true <= ci_hi)
  alpha_mat <- res[, c("a1","a2","a3","a4")]
  n_zero <- sum(rowSums(abs(alpha_mat), na.rm = TRUE) == 0, na.rm = TRUE)
  n_sep  <- sum(res[, "max_abs_alpha"] > 50, na.rm = TRUE)
  sd_a   <- apply(alpha_mat, 2, sd, na.rm = TRUE)
  cat(sprintf("  %-22s | stuck@0=%4d/%d sep(>50)=%3d | Bias=%+.4f ESD=%.4f ESE=%.4f ESE/ESD=%.3f CP=%.3f | sd(a)=(%.2f,%.2f,%.2f,%.2f) max|a|=%.1f (%.0fs)\n",
              optspec$tag, n_zero, n_reps, n_sep,
              mean(mu_c) - cfg$mu_true, sd(mu_c), mean(se_c),
              mean(se_c)/sd(mu_c), cp,
              sd_a[1], sd_a[2], sd_a[3], sd_a[4],
              max(res[,"max_abs_alpha"], na.rm = TRUE), out$elapsed))
  out_all[[optspec$tag]] <- list(res = res)
}
saveRDS(out_all, "Simulation_Test/sweep_cnr_s2_miss30_n2000_results.RDS")
cat("\n========= DONE =========\n")
