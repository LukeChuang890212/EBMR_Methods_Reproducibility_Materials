setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
source("Simulation.r")
library(EBMRalgorithmFast4)
library(parallel)
library(foreach)
library(doSNOW)

# Setting 4, scenario 9-3, model 3 only
mu_true <- get_mu_true("setting4")

ps_spec <- get_ps_spec("9-alt1")
# Model 3 only: u2_z2 formula, full h_alpha
single_ps_spec <- list(
  formula.list = ps_spec[["formula.list"]][3],
  h_alpha.list = ps_spec[["h_alpha.list"]][3],
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]]
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

n_val <- 2000
n_reps <- 1000

cat(sprintf("\n=== Setting4, scenario 9-3, model 3 only, n=%d, %d reps ===\n", n_val, n_reps))
cat(sprintf("mu_true = %.6f\n", mu_true))

cores <- max(1, detectCores() - 2)
cl <- makeCluster(cores)
registerDoSNOW(cl)
pb <- txtProgressBar(max = n_reps, style = 3)
opts <- list(progress = function(i) setTxtProgressBar(pb, i))

clusterExport(cl, c("setting4.B1", "single_ps_spec", "W_func", "n_val"), envir = environment())

results_raw <- foreach(rep_i = 1:n_reps,
                       .combine = cbind,
                       .options.snow = opts,
                       .packages = c("EBMRalgorithmFast4", "numDeriv"),
                       .errorhandling = "pass") %dopar% {
  set.seed(12345 + rep_i)
  dat <- setting4.B1(n_val)

  tryCatch({
    ebmr <- EBMRAlgorithmFast4[["new"]]("y", single_ps_spec, dat, W_func)
    res <- ebmr[["EBMR_IPW"]](
      h_nu = function(dat) cbind(u1 = dat[["u1"]], u2 = dat[["u2"]], z1 = dat[["z1"]], z2 = dat[["z2"]]),
      se.fit = TRUE, type = "HT"
    )
    c(mu_ipw = res[["mu_ipw"]], mu_ipw.true = 0, se_ipw = res[["se_ipw"]], se_ipw.true = 0)
  }, error = function(e) {
    rep(NA, 4)
  })
}

close(pb)
stopCluster(cl)

# Handle errors returned as list elements
if (is.list(results_raw)) {
  results_raw <- do.call(cbind, lapply(results_raw, function(x) if (is.numeric(x) && length(x) == 4) x else rep(NA, 4)))
}
rownames(results_raw) <- c("mu_ipw", "mu_ipw.true", "se_ipw", "se_ipw.true")

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
cat(sprintf("mu quantiles:\n"))
print(quantile(mu_all, c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1)))
cat(sprintf("se quantiles:\n"))
print(quantile(se_all, c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1)))

# --- AFTER cleaning ---
cleaned <- clean_sim_result(results_raw, multiplier = 3, verbose = TRUE)
rc <- cleaned$result
mu_c <- rc[1, ]
se_c <- rc[3, ]
esd_c <- sd(mu_c)
cat(sprintf("\n--- AFTER cleaning ---\n"))
cat(sprintf("N: %d (removed %d NA, %d outlier)\n", ncol(rc), cleaned$n_na, cleaned$n_outliers))
cat(sprintf("Bias: %.4f, ESD: %.4f\n", mean(mu_c) - mu_true, esd_c))
cat(sprintf("Mean ESE: %.4f, ESE/ESD: %.4f\n", mean(se_c), mean(se_c)/esd_c))

ci_l <- mu_c - 1.96 * se_c
ci_u <- mu_c + 1.96 * se_c
cp <- mean(mu_true >= ci_l & mu_true <= ci_u)
cat(sprintf("CP: %.4f\n", cp))
