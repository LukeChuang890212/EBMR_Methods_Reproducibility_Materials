setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Data_Generation.r")
source("config/scenarios.R")
library(EPS)
source("Basic_setup.r")

data_file <- misspecified_model_all_data_file.list$setting3$miss50[[2]]  # n=500
all_data <- readRDS(data_file)
n_val <- 500

ps_spec <- get_ps_spec("9-alt1")
subset_ps_spec <- list(
  formula.list = ps_spec$formula.list[2],
  h_alpha.list = ps_spec$h_alpha.list[2],
  inv_link = ps_spec$inv_link,
  outcome = ps_spec$outcome
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

# Type 1: all-NA rep (total failure)
cat("=== Type 1: All-NA rep (rep 93) ===\n")
dat <- all_data[((93 - 1) * n_val + 1):(93 * n_val), ]
tryCatch({
  ebmr <- EPS$new("y", subset_ps_spec, dat, W_func)
  cat("Constructor succeeded\n")
  res <- ebmr$EBMR_IPW(
    h_nu = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2),
    se.fit = TRUE, type = "HT"
  )
  cat(sprintf("mu_ipw=%.4f, se_ipw=%.4f\n", res$mu_ipw, res$se_ipw))
}, error = function(e) {
  cat(sprintf("ERROR: %s\n", conditionMessage(e)))
  cat(sprintf("Call: %s\n", deparse(conditionCall(e))))
})

# Try another all-NA rep
cat("\n=== Type 1: All-NA rep (rep 140) ===\n")
dat <- all_data[((140 - 1) * n_val + 1):(140 * n_val), ]
tryCatch({
  ebmr <- EPS$new("y", subset_ps_spec, dat, W_func)
  cat("Constructor succeeded\n")
  res <- ebmr$EBMR_IPW(
    h_nu = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2),
    se.fit = TRUE, type = "HT"
  )
  cat(sprintf("mu_ipw=%.4f, se_ipw=%.4f\n", res$mu_ipw, res$se_ipw))
}, error = function(e) {
  cat(sprintf("ERROR: %s\n", conditionMessage(e)))
  cat(sprintf("Call: %s\n", deparse(conditionCall(e))))
})

# Type 2: se-only NA rep
cat("\n=== Type 2: se_ipw-only NA rep (rep 327) ===\n")
dat <- all_data[((327 - 1) * n_val + 1):(327 * n_val), ]
tryCatch({
  ebmr <- EPS$new("y", subset_ps_spec, dat, W_func)
  cat("Constructor succeeded\n")

  # First without SE
  res_no_se <- ebmr$EBMR_IPW(
    h_nu = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2),
    se.fit = FALSE, type = "HT"
  )
  cat(sprintf("Without SE: mu_ipw=%.4f\n", res_no_se$mu_ipw))

  # Then with SE
  res <- ebmr$EBMR_IPW(
    h_nu = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2),
    se.fit = TRUE, type = "HT"
  )
  cat(sprintf("With SE: mu_ipw=%.4f, se_ipw=%.4f\n", res$mu_ipw, res$se_ipw))
}, error = function(e) {
  cat(sprintf("ERROR: %s\n", conditionMessage(e)))
  cat(sprintf("Call: %s\n", deparse(conditionCall(e))))
})

cat("\nDone!\n")
