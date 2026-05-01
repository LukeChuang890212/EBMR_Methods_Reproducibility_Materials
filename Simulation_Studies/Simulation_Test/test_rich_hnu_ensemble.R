## Test S4 M2+M3 ensemble with enriched h_alpha AND h_nu
## h_alpha = (1, u1, u2, z1, z2, u1u2, u1z1, u1z2, u2z1, u2z2, z1z2)
## h_nu   = (u1, u2, z1, z2, u1u2, u1z1, u1z2, u2z1, u2z2, z1z2, u1u2z1)
## Uses the package directly
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)
library(parallel)
library(foreach)
library(doSNOW)

n_val <- 2000; n_reps <- 1000
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

## Get base ps_spec and override h_alpha with enriched version
ps_spec_base <- get_ps_spec("9-alt1")

h_alpha_rich <- function(dat) cbind(
  u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2,
  u1u2 = dat$u1*dat$u2, u1z1 = dat$u1*dat$z1, u1z2 = dat$u1*dat$z2,
  u2z1 = dat$u2*dat$z1, u2z2 = dat$u2*dat$z2, z1z2 = dat$z1*dat$z2
)

ps_spec <- list(
  formula.list = ps_spec_base$formula.list,
  h_alpha.list = list(h_alpha_rich, h_alpha_rich, h_alpha_rich),
  inv_link = ps_spec_base$inv_link,
  outcome = ps_spec_base$outcome
)

# Note: ensemble() auto-adds intercept, so h_nu should NOT include 1
h_nu_fn <- function(dat) cbind(
  u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2,
  u1u2 = dat$u1*dat$u2, u1z1 = dat$u1*dat$z1, u1z2 = dat$u1*dat$z2,
  u2z1 = dat$u2*dat$z1, u2z2 = dat$u2*dat$z2, z1z2 = dat$z1*dat$z2,
  u2sq = dat$u2^2, z2sq = dat$z2^2,
  u1u2z1 = dat$u1*dat$u2*dat$z1, u1u2z2 = dat$u1*dat$u2*dat$z2,
  u1z1z2 = dat$u1*dat$z1*dat$z2, u2z1z2 = dat$u2*dat$z1*dat$z2,
  u2sq_u1 = dat$u2^2*dat$u1, u2sq_z1 = dat$u2^2*dat$z1, u2sq_z2 = dat$u2^2*dat$z2,
  z2sq_u1 = dat$z2^2*dat$u1, z2sq_u2 = dat$z2^2*dat$u2, z2sq_z1 = dat$z2^2*dat$z1,
  u2cu = dat$u2^3, z2cu = dat$z2^3
)

## Report function using clean_sim_result
report_summary <- function(sim_result, mu_true, pe_idx, se_idx, label) {
  n_total <- ncol(sim_result)

  # --- Raw (no cleaning) ---
  na_mask <- apply(sim_result[c(pe_idx, se_idx), , drop = FALSE], 2, function(x) !any(is.na(x)))
  raw <- sim_result[, na_mask, drop = FALSE]
  n_valid <- ncol(raw)
  mu_v <- raw[pe_idx, ]; se_v <- raw[se_idx, ]

  cat(sprintf("  --- %s (raw, n=%d/%d) ---\n", label, n_valid, n_total))
  cat(sprintf("  Bias     = %.4f\n", mean(mu_v) - mu_true))
  cat(sprintf("  ESD      = %.4f\n", sd(mu_v)))
  cat(sprintf("  ESE      = %.4f\n", mean(se_v)))
  cat(sprintf("  ESE/ESD  = %.3f\n", mean(se_v) / sd(mu_v)))
  ci_lo <- mu_v - 1.96 * se_v; ci_hi <- mu_v + 1.96 * se_v
  cat(sprintf("  CP       = %.3f\n", mean(mu_true >= ci_lo & mu_true <= ci_hi)))

  # --- After clean_sim_result ---
  cleaned <- clean_sim_result(sim_result, multiplier = 3, verbose = FALSE)
  sim_c <- cleaned$result
  mu_c <- sim_c[pe_idx, ]; se_c <- sim_c[se_idx, ]

  cat(sprintf("\n  --- %s (cleaned: NA=%d, outliers=%d, n=%d) ---\n",
      label, cleaned$n_na, cleaned$n_outliers, cleaned$n_successful))
  cat(sprintf("  Bias     = %.4f\n", mean(mu_c) - mu_true))
  cat(sprintf("  ESD      = %.4f\n", sd(mu_c)))
  cat(sprintf("  ESE      = %.4f\n", mean(se_c)))
  cat(sprintf("  ESE/ESD  = %.3f\n", mean(se_c) / sd(mu_c)))
  ci_lo <- mu_c - 1.96 * se_c; ci_hi <- mu_c + 1.96 * se_c
  cat(sprintf("  CP       = %.3f\n", mean(mu_true >= ci_lo & mu_true <= ci_hi)))
  cat("\n")
}

cat("=== S4 M2+M3 ensemble: enriched h_alpha (11 instr) + enriched h_nu (12 instr) ===\n\n")

setting <- "setting4"
data_file <- misspecified_model_all_data_file.list[[setting]][["miss50"]][[1]]
all_data <- readRDS(data_file)
mu_true <- get_mu_true(setting)
alpha_true <- misspecified_model_alpha.true.list[[setting]][["miss50"]][[1]]
ps_model_true <- function(dat, alpha.true) {
  X <- cbind(rep(1, nrow(dat)), dat$y, dat$u1, dat$u2)
  1 / (1 + exp(X %*% alpha.true))
}
cat(sprintf("True mean: %.6f\n", mu_true))

n_cores <- min(detectCores() - 1, 10)
cl <- makeCluster(n_cores)
registerDoSNOW(cl)
clusterExport(cl, c("ps_spec", "W_fn", "h_nu_fn", "all_data", "n_val",
                     "ps_model_true", "alpha_true", "h_alpha_rich"),
              envir = environment())

pb <- txtProgressBar(max = n_reps, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

sim_result <- foreach(
  i = 1:n_reps,
  .combine = 'cbind',
  .options.snow = opts,
  .packages = c("EBMRalgorithmFast4", "stringr", "Matrix", "numDeriv")
) %dopar% {
  tryCatch({
    dat <- all_data[((i-1)*n_val + 1):(i*n_val), ]
    ebmr <- EBMRAlgorithmFast4$new("y", ps_spec, dat, W_fn)
    result <- ebmr$EBMR_IPW(
      h_nu = h_nu_fn,
      model_indices = c(2, 3),
      true_ps = ps_model_true(dat, alpha_true),
      type = "HT"
    )
    c(unlist(result[1:4]), nu1 = result$nu.hat[1], nu2 = result$nu.hat[2],
      w1 = result$w.hat[1], w2 = result$w.hat[2])
  }, error = function(e) {
    rep(NA, 8)
  })
}

close(pb)
stopCluster(cl)

# sim_result rows: mu_ipw(1), mu_ipw.true(2), se_ipw(3), se_ipw.true(4)
cat("\n")
report_summary(sim_result, mu_true, pe_idx = 1, se_idx = 3, label = "M2+M3 IPW (estimated PS)")
report_summary(sim_result, mu_true, pe_idx = 2, se_idx = 4, label = "M2+M3 IPW (true PS)")

# Nu and weight diagnostics
nu1 <- sim_result[5, ]; nu2 <- sim_result[6, ]
w1 <- sim_result[7, ]; w2 <- sim_result[8, ]
valid_w <- !is.na(w1)
nu1v <- nu1[valid_w]; nu2v <- nu2[valid_w]
w1v <- w1[valid_w]; w2v <- w2[valid_w]

cat("=== Nu diagnostics ===\n")
cat(sprintf("  nu1 (M2): mean=%.4f, sd=%.4f, min=%.4f, Q1=%.4f, Q5=%.4f, Q25=%.4f, med=%.4f, Q75=%.4f, Q95=%.4f, Q99=%.4f, max=%.4f\n",
    mean(nu1v), sd(nu1v), min(nu1v), quantile(nu1v,0.01), quantile(nu1v,0.05),
    quantile(nu1v,0.25), median(nu1v), quantile(nu1v,0.75), quantile(nu1v,0.95), quantile(nu1v,0.99), max(nu1v)))
cat(sprintf("  nu2 (M3): mean=%.4f, sd=%.4f, min=%.4f, Q1=%.4f, Q5=%.4f, Q25=%.4f, med=%.4f, Q75=%.4f, Q95=%.4f, Q99=%.4f, max=%.4f\n",
    mean(nu2v), sd(nu2v), min(nu2v), quantile(nu2v,0.01), quantile(nu2v,0.05),
    quantile(nu2v,0.25), median(nu2v), quantile(nu2v,0.75), quantile(nu2v,0.95), quantile(nu2v,0.99), max(nu2v)))

cat("\n=== Weight diagnostics ===\n")
cat(sprintf("  w1 (M2): mean=%.4f, sd=%.4f, min=%.4f, Q25=%.4f, med=%.4f, Q75=%.4f, max=%.4f\n",
    mean(w1v), sd(w1v), min(w1v), quantile(w1v,0.25), median(w1v), quantile(w1v,0.75), max(w1v)))
cat(sprintf("  w2 (M3): mean=%.4f, sd=%.4f, min=%.4f, Q25=%.4f, med=%.4f, Q75=%.4f, max=%.4f\n",
    mean(w2v), sd(w2v), min(w2v), quantile(w2v,0.25), median(w2v), quantile(w2v,0.75), max(w2v)))
cat(sprintf("  Fraction w1 > 0.5: %.3f\n", mean(w1v > 0.5)))
cat(sprintf("  Fraction w1 > 0.9: %.3f\n", mean(w1v > 0.9)))
cat(sprintf("  Fraction w2 > 0.9: %.3f\n", mean(w2v > 0.9)))
mu_v <- sim_result[1, valid_w]
cat(sprintf("  cor(w1, mu): %.3f\n", cor(w1v, mu_v, use="complete.obs")))
cat(sprintf("  cor(w2, mu): %.3f\n", cor(w2v, mu_v, use="complete.obs")))
