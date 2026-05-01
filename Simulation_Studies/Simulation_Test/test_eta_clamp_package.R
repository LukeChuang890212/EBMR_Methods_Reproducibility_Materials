## Test: Verify eta clamping in package on S3, scenario 9-2, model 2, n=500, miss30
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)
library(parallel)
library(foreach)
library(doSNOW)

n_val <- 500; n_reps <- 1000
ps_spec <- get_ps_spec("9")
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
mu_true <- get_mu_true("setting3")

h_nu <- function(dat) cbind(
  u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2,
  u1u2 = dat$u1*dat$u2
)

all_data <- readRDS("Simulation_Data/Setting3.B2_n500_replicate1000.RDS")

cat(sprintf("=== Package test: S3, 9-2, M2, n=500, miss30 (with eta clamping) ===\n"))
cat(sprintf("True mean: %.4f\n\n", mu_true))

# Use model 2 only
m_idx <- 2
single_ps <- list(
  formula.list = list(ps_spec[["formula.list"]][[m_idx]]),
  h_alpha.list = list(ps_spec[["h_alpha.list"]][[m_idx]]),
  inv_link = ps_spec[["inv_link"]], outcome = ps_spec[["outcome"]]
)

n_cores <- min(detectCores() - 1, 10)
cl <- makeCluster(n_cores)
registerDoSNOW(cl)
clusterExport(cl, c("all_data", "n_val", "single_ps", "W_fn", "h_nu"),
              envir = environment())
clusterEvalQ(cl, devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE))

pb <- txtProgressBar(max = n_reps, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

sim_result <- foreach(
  i = 1:n_reps,
  .combine = 'cbind',
  .options.snow = opts,
  .packages = c("stringr", "Matrix", "numDeriv")
) %dopar% {
  tryCatch({
    dat <- all_data[((i-1)*n_val + 1):(i*n_val), ]
    ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
    result <- ebmr$EBMR_IPW(h_nu = h_nu, true_ps = NULL, type = "HT")
    unlist(result[1:4])
  }, error = function(e) {
    rep(NA, 4)
  })
}

close(pb)
stopCluster(cl)

# Report — only check rows 1 and 3 (mu_ipw and se_ipw) since true_ps=NULL makes rows 2,4 NA
cat("\n\n=== Results (package with eta clamping) ===\n\n")

mu_all <- sim_result[1, ]; se_all <- sim_result[3, ]
valid <- !is.na(mu_all) & !is.na(se_all)
n_na <- sum(!valid)
mu_v0 <- mu_all[valid]; se_v0 <- se_all[valid]

# Outlier removal (IQR x 3, cap 1%)
Q1 <- quantile(mu_v0, 0.25); Q3 <- quantile(mu_v0, 0.75); IQR_val <- Q3 - Q1
is_out <- mu_v0 < (Q1 - 3*IQR_val) | mu_v0 > (Q3 + 3*IQR_val)
max_rm <- floor(0.01 * length(mu_v0))
if (sum(is_out) > max_rm && max_rm > 0) {
  dm <- abs(mu_v0 - median(mu_v0))
  oi <- which(is_out)
  keep <- oi[order(dm[oi], decreasing=TRUE)[(max_rm+1):length(oi)]]
  is_out[keep] <- FALSE
}
n_out <- sum(is_out)
mu_v <- mu_v0[!is_out]; se_v <- se_v0[!is_out]

cat(sprintf("  Total: %d, NA: %d, Outliers: %d, Used: %d\n",
    n_reps, n_na, n_out, length(mu_v)))
cat(sprintf("  Bias     = %.4f\n", mean(mu_v) - mu_true))
cat(sprintf("  ESD      = %.4f\n", sd(mu_v)))
cat(sprintf("  ESE      = %.4f\n", mean(se_v)))
cat(sprintf("  ESE/ESD  = %.3f\n", mean(se_v) / sd(mu_v)))
ci_lo <- mu_v - 1.96 * se_v; ci_hi <- mu_v + 1.96 * se_v
cat(sprintf("  CP       = %.3f\n", mean(mu_true >= ci_lo & mu_true <= ci_hi)))

# Compare with original
cat("\n=== Original (no clamping, from saved file) ===\n\n")
sim_orig <- readRDS("Simulation_Results/EBMR_IPW_setting3-miss30-scenario9-2_2_n500_replicate1000_test59.RDS")
cleaned_orig <- clean_sim_result(sim_orig, multiplier = 3, verbose = TRUE)
sim_co <- cleaned_orig$result

cat(sprintf("  Total: %d, NA: %d, Outliers: %d, Used: %d\n",
    cleaned_orig$n_total, cleaned_orig$n_na, cleaned_orig$n_outliers, cleaned_orig$n_successful))

mu_o <- sim_co[1, ]; se_o <- sim_co[3, ]
cat(sprintf("  Bias     = %.4f\n", mean(mu_o) - mu_true))
cat(sprintf("  ESD      = %.4f\n", sd(mu_o)))
cat(sprintf("  ESE      = %.4f\n", mean(se_o)))
cat(sprintf("  ESE/ESD  = %.3f\n", mean(se_o) / sd(mu_o)))
ci_lo <- mu_o - 1.96 * se_o; ci_hi <- mu_o + 1.96 * se_o
cat(sprintf("  CP       = %.3f\n", mean(mu_true >= ci_lo & mu_true <= ci_hi)))
