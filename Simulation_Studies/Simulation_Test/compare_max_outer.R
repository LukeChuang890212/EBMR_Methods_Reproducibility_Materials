setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
source("Data_Generation.r")
source("config/scenarios.R")

old_wd <- getwd()
setwd("../EBMRalgorithmFast4")
source("R/EBMRAlgorithm.r")
setwd(old_wd)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
inv_link_fn <- function(eta) 1 / (1 + exp(eta))

ps_spec_m2 <- list(
  formula.list = list(FORMULAS$u1_z2),
  h_alpha.list = list(H_ALPHA$full),
  inv_link = inv_link_fn,
  outcome = "y",
  alpha_init.list = list(NULL)
)

h_nu_func <- function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1,
                                  z2 = dat$z2, u1_u2 = dat$u1 * dat$u2)

set.seed(999)
big <- setting3.A1(1e6)
mu_true <- mean(big$y)
rm(big); gc()
cat(sprintf("mu_true = %.6f\n\n", mu_true))

set.seed(42)
n_reps <- 200
n_val <- 2000
datasets <- vector("list", n_reps)
for (i in 1:n_reps) datasets[[i]] <- setting3.B1(n_val)

gmm_source <- readLines("../EBMRalgorithmFast4/R/GMM_functions.R")

max_outer_values <- c(10, 20, 50, 100, 200, 500, 1000)

for (mo in max_outer_values) {
  cat(sprintf("\n========== max_outer = %d ==========\n", mo))

  patched <- gsub("max_outer = 200L", sprintf("max_outer = %dL", mo), gmm_source)
  tmp_file <- tempfile(fileext = ".R")
  writeLines(patched, tmp_file)
  source(tmp_file)
  source("../EBMRalgorithmFast4/R/Methods.r")
  old_wd2 <- getwd()
  setwd("../EBMRalgorithmFast4")
  source("R/EBMRAlgorithm.r")
  setwd(old_wd2)

  mu_hat <- numeric(n_reps); se_hat <- numeric(n_reps)
  converged <- logical(n_reps); iterations <- integer(n_reps)
  gnorm_vec <- numeric(n_reps)

  t_start <- proc.time()
  for (i in 1:n_reps) {
    dat <- datasets[[i]]
    tryCatch({
      ebmr <- EBMRAlgorithmFast4$new("y", ps_spec_m2, dat, W_func, bounds = Inf)
      ipw_result <- ebmr$EBMR_IPW(h_nu = h_nu_func)
      fit <- ebmr$ps_fit.list[[1]]
      gmm_fit <- fit$gmm_fit

      mu_hat[i] <- ipw_result$mu_ipw
      se_hat[i] <- ipw_result$se_ipw
      converged[i] <- gmm_fit$opt$converged
      iterations[i] <- gmm_fit$opt$iterations
      gnorm_vec[i] <- ifelse(is.null(gmm_fit$opt$final_grad_norm), NA, gmm_fit$opt$final_grad_norm)
    }, error = function(e) {
      mu_hat[i] <<- NA; se_hat[i] <<- NA; converged[i] <<- NA
      cat(sprintf("Rep %d: ERROR - %s\n", i, conditionMessage(e)))
    })
  }
  elapsed <- (proc.time() - t_start)[3]

  valid <- !is.na(mu_hat) & !is.na(se_hat)
  mu_v <- mu_hat[valid]; se_v <- se_hat[valid]
  q1 <- quantile(mu_v, 0.25); q3 <- quantile(mu_v, 0.75); iqr <- q3 - q1
  is_outlier <- mu_v < (q1 - 3*iqr) | mu_v > (q3 + 3*iqr)
  mu_v <- mu_v[!is_outlier]; se_v <- se_v[!is_outlier]; n_kept <- length(mu_v)
  empirical_sd <- sd(mu_v)
  valid_idx <- which(valid); kept_idx <- valid_idx[!is_outlier]
  ci_lower <- mu_v - 1.96*se_v; ci_upper <- mu_v + 1.96*se_v
  covers <- (ci_lower <= mu_true) & (mu_true <= ci_upper)
  fgn <- gnorm_vec[kept_idx]; fgn <- fgn[!is.na(fgn)]

  cat(sprintf("  Time: %.1fs | Outliers: %d | Kept: %d\n", elapsed, sum(is_outlier), n_kept))
  cat(sprintf("  Bias: %.4f | RMSE: %.4f | SD: %.4f\n",
              mean(mu_v) - mu_true, sqrt(mean((mu_v - mu_true)^2)), empirical_sd))
  cat(sprintf("  mean(SE): %.4f | mean(SE)/SD: %.4f\n", mean(se_v), mean(se_v)/empirical_sd))
  cat(sprintf("  Coverage: %.1f%% | Conv: %.0f%%\n", 100*mean(covers), 100*mean(converged[kept_idx])))
  cat(sprintf("  Iters: mean=%.1f, median=%.0f, range=[%d, %d]\n",
              mean(iterations[kept_idx]), median(iterations[kept_idx]),
              min(iterations[kept_idx]), max(iterations[kept_idx])))
  if (length(fgn) > 0) cat(sprintf("  Grad norm: median=%.2e, max=%.2e\n", median(fgn), max(fgn)))

  file.remove(tmp_file)
}

cat("\nDone!\n")
