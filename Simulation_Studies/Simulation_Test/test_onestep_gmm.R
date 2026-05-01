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

# Use W = I (identity) — one-step GMM, no efficient weighting
W_identity <- function(g.matrix) diag(ncol(g.matrix))
# Use W = optimal (two-step)
W_optimal <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

mu.true <- get_mu_true("setting3")

n_reps <- 1000
results_identity <- matrix(NA, 4, n_reps)
rownames(results_identity) <- c("mu_ipw", "mu_ipw.true", "se_ipw", "se_ipw.true")

results_optimal <- matrix(NA, 4, n_reps)
rownames(results_optimal) <- c("mu_ipw", "mu_ipw.true", "se_ipw", "se_ipw.true")

cat("Running 1000 reps with W=I (one-step) and W=optimal (two-step)...\n")
n_err_I <- 0
n_err_O <- 0
for (rep_i in 1:n_reps) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]

  # W = I
  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", subset_ps_spec, dat, W_identity)
    res <- ebmr$EBMR_IPW(
      h_nu = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2),
      se.fit = TRUE, type = "HT"
    )
    results_identity[1, rep_i] <- res$mu_ipw
    results_identity[3, rep_i] <- res$se_ipw
  }, error = function(e) {
    n_err_I <<- n_err_I + 1
  })

  # W = optimal (two-step)
  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", subset_ps_spec, dat, W_optimal)
    res <- ebmr$EBMR_IPW(
      h_nu = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2),
      se.fit = TRUE, type = "HT"
    )
    results_optimal[1, rep_i] <- res$mu_ipw
    results_optimal[3, rep_i] <- res$se_ipw
  }, error = function(e) {
    n_err_O <<- n_err_O + 1
  })

  if (rep_i %% 200 == 0) cat(sprintf("  %d / %d done (err: I=%d, O=%d)\n", rep_i, n_reps, n_err_I, n_err_O))
}

results_identity[2, ] <- ifelse(is.na(results_identity[1, ]), NA, 0)
results_identity[4, ] <- ifelse(is.na(results_identity[3, ]), NA, 0)
results_optimal[2, ] <- ifelse(is.na(results_optimal[1, ]), NA, 0)
results_optimal[4, ] <- ifelse(is.na(results_optimal[3, ]), NA, 0)

# Report for both
for (label in c("W=I (one-step)", "W=optimal (two-step)")) {
  if (label == "W=I (one-step)") {
    res <- results_identity
  } else {
    res <- results_optimal
  }

  for (mult in c(1.5, 3)) {
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

    cat(sprintf("\n%s, mult=%.1f (n=%d, removed=%d):\n", label, mult, length(mu_v), n_reps - length(mu_v) - sum(is.na(res[1,]))))
    cat(sprintf("  Bias = %.4f, ESD = %.4f, ESE = %.4f, ESE/ESD = %.2f, CP = %.4f\n",
                bias, esd, ese, ese/esd, cp))
  }
}

cat("\nDone!\n")
