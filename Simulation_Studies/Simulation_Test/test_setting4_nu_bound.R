setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
source("Simulation.r")
library(EBMRalgorithmFast4)
library(parallel)
library(foreach)
library(doSNOW)

# Setting 4, scenario 9-3, model combination "13" (models 1 and 3)
mu_true <- get_mu_true("setting4")

ps_spec <- get_ps_spec("9-alt1")
ps_spec_13 <- list(
  formula.list = ps_spec[["formula.list"]][c(1, 3)],
  h_alpha.list = ps_spec[["h_alpha.list"]][c(1, 3)],
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]]
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

n_val <- 2000
n_reps <- 1000

cat(sprintf("\n=== Setting4, models 1&3, n=%d, %d reps (nu >= 0 bound) ===\n", n_val, n_reps))
cat(sprintf("mu_true = %.6f\n", mu_true))

cores <- max(1, detectCores() - 2)
cl <- makeCluster(cores)
registerDoSNOW(cl)
pb <- txtProgressBar(max = n_reps, style = 3)
opts <- list(progress = function(i) setTxtProgressBar(pb, i))

clusterExport(cl, c("setting4.B1", "ps_spec_13", "W_func", "n_val"), envir = environment())

results_raw <- foreach(rep_i = 1:n_reps,
                       .combine = cbind,
                       .options.snow = opts,
                       .packages = c("EBMRalgorithmFast4", "numDeriv"),
                       .errorhandling = "pass") %dopar% {
  set.seed(12345 + rep_i)
  dat <- setting4.B1(n_val)

  tryCatch({
    ebmr <- EBMRAlgorithmFast4[["new"]]("y", ps_spec_13, dat, W_func)
    res <- ebmr[["EBMR_IPW"]](
      h_nu = function(dat) cbind(u1 = dat[["u1"]], u2 = dat[["u2"]], z1 = dat[["z1"]], z2 = dat[["z2"]], u1_u2 = dat[["u1"]]*dat[["u2"]]),
      se.fit = TRUE, type = "HT"
    )
    nu_hat <- res[["nu.hat"]]
    w_hat <- res[["w.hat"]]
    c(mu_ipw = res[["mu_ipw"]], mu_ipw.true = 0, se_ipw = res[["se_ipw"]], se_ipw.true = 0,
      nu1 = nu_hat[1], nu2 = nu_hat[2], w1 = w_hat[1], w2 = w_hat[2])
  }, error = function(e) {
    rep(NA, 8)
  })
}

close(pb)
stopCluster(cl)

# Handle errors returned as list elements
if (is.list(results_raw)) {
  results_raw <- do.call(cbind, lapply(results_raw, function(x) if (is.numeric(x) && length(x) == 8) x else rep(NA, 8)))
}
rownames(results_raw) <- c("mu_ipw", "mu_ipw.true", "se_ipw", "se_ipw.true", "nu1", "nu2", "w1", "w2")

n_err <- sum(is.na(results_raw[1, ]))
cat(sprintf("\nErrors: %d / %d\n", n_err, n_reps))

# --- BEFORE cleaning ---
valid <- !is.na(results_raw[1, ])
mu_all <- results_raw[1, valid]
se_all <- results_raw[3, valid]
esd_all <- sd(mu_all)
cat(sprintf("\n--- BEFORE cleaning ---\n"))
cat(sprintf("N: %d, Bias: %.4f, ESD: %.4f\n", sum(valid), mean(mu_all) - mu_true, esd_all))
cat(sprintf("Mean ESE: %.4f, ESE/ESD: %.4f\n", mean(se_all), mean(se_all)/esd_all))

# Nu/w distribution
nu1 <- results_raw[5, valid]; nu2 <- results_raw[6, valid]
w1 <- results_raw[7, valid]; w2 <- results_raw[8, valid]
cat(sprintf("\nnu1: mean=%.4f, median=%.4f, min=%.4f, max=%.4f\n", mean(nu1), median(nu1), min(nu1), max(nu1)))
cat(sprintf("nu2: mean=%.4f, median=%.4f, min=%.4f, max=%.4f\n", mean(nu2), median(nu2), min(nu2), max(nu2)))
cat(sprintf("w1:  mean=%.4f, median=%.4f\n", mean(w1), median(w1)))
cat(sprintf("w2:  mean=%.4f, median=%.4f\n", mean(w2), median(w2)))
cat(sprintf("nu1==0: %d, nu2==0: %d\n", sum(nu1 == 0), sum(nu2 == 0)))

# --- AFTER cleaning ---
cleaned <- clean_sim_result(results_raw[1:4, ], multiplier = 3, verbose = TRUE)
rc_idx <- !is.na(cleaned$result[1, ])
# Map back to full results
na_mask <- apply(results_raw[1:4, , drop = FALSE], 2, function(x) !any(is.na(x)))
results_nona <- results_raw[, na_mask, drop = FALSE]
mu_ipw <- results_nona[1, ]
Q1 <- quantile(mu_ipw, 0.25); Q3 <- quantile(mu_ipw, 0.75)
IQR_value <- Q3 - Q1
is_outlier <- mu_ipw < (Q1 - 3 * IQR_value) | mu_ipw > (Q3 + 3 * IQR_value)
max_remove <- floor(0.01 * length(mu_ipw))
if (sum(is_outlier) > max_remove && max_remove > 0) {
  dist_from_median <- abs(mu_ipw - median(mu_ipw))
  outlier_idx <- which(is_outlier)
  keep_idx <- outlier_idx[order(dist_from_median[outlier_idx], decreasing = TRUE)[(max_remove + 1):length(outlier_idx)]]
  is_outlier[keep_idx] <- FALSE
}
rc <- results_nona[, !is_outlier, drop = FALSE]

mu_c <- rc[1, ]
se_c <- rc[3, ]
esd_c <- sd(mu_c)
cat(sprintf("\n--- AFTER cleaning ---\n"))
cat(sprintf("N: %d (removed %d NA, %d outlier)\n", ncol(rc), sum(!na_mask), sum(is_outlier)))
cat(sprintf("Bias: %.4f, ESD: %.4f\n", mean(mu_c) - mu_true, esd_c))
cat(sprintf("Mean ESE: %.4f, ESE/ESD: %.4f\n", mean(se_c), mean(se_c)/esd_c))

ci_l <- mu_c - 1.96 * se_c
ci_u <- mu_c + 1.96 * se_c
cp <- mean(mu_true >= ci_l & mu_true <= ci_u)
cat(sprintf("CP: %.4f\n", cp))
