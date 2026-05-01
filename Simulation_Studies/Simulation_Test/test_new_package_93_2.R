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

n_reps <- 1000
results <- matrix(NA, 4, n_reps)
rownames(results) <- c("mu_ipw", "mu_ipw.true", "se_ipw", "se_ipw.true")

mu.true <- get_mu_true("setting3")

# True PS for this setting
alpha_true <- correct_model_alpha.true.list$setting3$miss50[[1]]

cat("Running 1000 reps with new package (eta_max=20)...\n")
n_err <- 0
for (rep_i in 1:n_reps) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]

  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", subset_ps_spec, dat, W_func)
    res <- ebmr$EBMR_IPW(
      h_nu = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2),
      se.fit = TRUE, type = "HT"
    )
    results[1, rep_i] <- res$mu_ipw
    results[3, rep_i] <- res$se_ipw
  }, error = function(e) {
    n_err <<- n_err + 1
  })

  if (rep_i %% 200 == 0) cat(sprintf("  %d / %d done (errors: %d)\n", rep_i, n_reps, n_err))
}

cat(sprintf("\nErrors: %d / %d\n", n_err, n_reps))

# Fill rows 2,4 with 0 so clean_sim_result doesn't drop them as NA
results[2, ] <- ifelse(is.na(results[1, ]), NA, 0)
results[4, ] <- ifelse(is.na(results[3, ]), NA, 0)

# Apply clean_sim_result
cleaned <- clean_sim_result(results, verbose = TRUE)
sim <- cleaned$result

mu_v <- sim[1, ]
se_v <- sim[3, ]

bias <- mean(mu_v) - mu.true
esd <- sd(mu_v)
ese <- mean(se_v)
ci_lower <- mu_v - 1.96 * se_v
ci_upper <- mu_v + 1.96 * se_v
cp <- mean(mu.true >= ci_lower & mu.true <= ci_upper)

cat(sprintf("\nmu.true = %.4f\n", mu.true))
cat(sprintf("Bias    = %.4f\n", bias))
cat(sprintf("ESD     = %.4f\n", esd))
cat(sprintf("ESE     = %.4f\n", ese))
cat(sprintf("ESE/ESD = %.2f\n", ese / esd))
cat(sprintf("CP      = %.4f\n", cp))
