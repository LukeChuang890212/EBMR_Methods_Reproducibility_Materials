setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
source("Simulation.r")
library(EBMRalgorithmFast4)

data_file <- misspecified_model_all_data_file.list$setting3$miss50[[1]]
all_data <- readRDS(data_file)
n_val <- 2000

ps_spec <- get_ps_spec("9-alt1")
subset_ps_spec <- list(
  formula.list = ps_spec$formula.list[2],
  h_alpha.list = ps_spec$h_alpha.list[2],
  inv_link = ps_spec$inv_link,
  outcome = ps_spec$outcome
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
mu.true <- get_mu_true("setting3")

n_reps <- 1000

# Storage for both methods
results_twostep <- matrix(NA, 4, n_reps)
results_iterative <- matrix(NA, 4, n_reps)
rownames(results_twostep) <- rownames(results_iterative) <-
  c("mu_ipw", "mu_ipw.true", "se_ipw", "se_ipw.true")

cat("Running 1000 reps for both two-step and iterative GMM...\n")
n_err_ts <- 0
n_err_it <- 0
for (rep_i in 1:n_reps) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]

  # Two-step
  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", subset_ps_spec, dat, W_func, gmm_type = "twostep")
    res <- ebmr$EBMR_IPW(
      h_nu = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2),
      se.fit = TRUE, type = "HT"
    )
    results_twostep[1, rep_i] <- res$mu_ipw
    results_twostep[3, rep_i] <- res$se_ipw
  }, error = function(e) { n_err_ts <<- n_err_ts + 1 })

  # Iterative
  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", subset_ps_spec, dat, W_func, gmm_type = "iterative")
    res <- ebmr$EBMR_IPW(
      h_nu = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2),
      se.fit = TRUE, type = "HT"
    )
    results_iterative[1, rep_i] <- res$mu_ipw
    results_iterative[3, rep_i] <- res$se_ipw
  }, error = function(e) { n_err_it <<- n_err_it + 1 })

  if (rep_i %% 200 == 0) cat(sprintf("  %d / %d done (err: ts=%d, it=%d)\n",
                                      rep_i, n_reps, n_err_ts, n_err_it))
}

results_twostep[2, ] <- ifelse(is.na(results_twostep[1, ]), NA, 0)
results_twostep[4, ] <- ifelse(is.na(results_twostep[3, ]), NA, 0)
results_iterative[2, ] <- ifelse(is.na(results_iterative[1, ]), NA, 0)
results_iterative[4, ] <- ifelse(is.na(results_iterative[3, ]), NA, 0)

cat(sprintf("\nErrors: twostep=%d, iterative=%d\n", n_err_ts, n_err_it))

# Report: before outlier removal (raw valid reps)
cat("\n===== BEFORE OUTLIER REMOVAL =====\n")
for (label in c("Two-step", "Iterative")) {
  res <- if (label == "Two-step") results_twostep else results_iterative
  valid <- !is.na(res[1, ])
  mu_v <- res[1, valid]
  se_v <- res[3, valid]

  bias <- mean(mu_v) - mu.true
  esd <- sd(mu_v)
  ese <- mean(se_v)
  ci_lower <- mu_v - 1.96 * se_v
  ci_upper <- mu_v + 1.96 * se_v
  cp <- mean(mu.true >= ci_lower & mu.true <= ci_upper)

  cat(sprintf("\n%s (n=%d):\n", label, length(mu_v)))
  cat(sprintf("  Bias = %.4f, ESD = %.4f, ESE = %.4f, ESE/ESD = %.2f, CP = %.4f\n",
              bias, esd, ese, ese/esd, cp))
}

# Report: after outlier removal with different multipliers
for (mult in c(1.5, 3)) {
  cat(sprintf("\n===== AFTER OUTLIER REMOVAL (mult=%.1f) =====\n", mult))
  for (label in c("Two-step", "Iterative")) {
    res <- if (label == "Two-step") results_twostep else results_iterative
    cleaned <- clean_sim_result(res, multiplier = mult, verbose = FALSE)
    sim <- cleaned$result
    mu_v <- sim[1, ]
    se_v <- sim[3, ]

    bias <- mean(mu_v) - mu.true
    esd <- sd(mu_v)
    ese <- mean(se_v)
    ci_lower <- mu_v - 1.96 * se_v
    ci_upper <- mu_v + 1.96 * se_v
    cp <- mean(mu.true >= ci_lower & mu.true <= ci_upper)

    n_removed <- sum(!is.na(res[1, ])) - length(mu_v)
    cat(sprintf("\n%s (n=%d, removed=%d):\n", label, length(mu_v), n_removed))
    cat(sprintf("  Bias = %.4f, ESD = %.4f, ESE = %.4f, ESE/ESD = %.2f, CP = %.4f\n",
                bias, esd, ese, ese/esd, cp))
  }
}

cat("\nDone!\n")
