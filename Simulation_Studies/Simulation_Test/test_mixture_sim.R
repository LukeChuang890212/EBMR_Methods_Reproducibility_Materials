setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
source("Simulation.r")
library(EBMRalgorithmFast4)

# Setting 3, model 2 alone
data_file <- misspecified_model_all_data_file.list[["setting3"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
n_val <- 2000
mu_true <- get_mu_true("setting3")

ps_spec <- get_ps_spec("9-alt1")
single_ps_spec <- list(
  formula.list = ps_spec[["formula.list"]][2],
  h_alpha.list = ps_spec[["h_alpha.list"]][2],
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]]
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

n_reps <- 1000
results <- matrix(NA, 4, n_reps)
rownames(results) <- c("mu_ipw", "mu_ipw.true", "se_ipw", "se_ipw.true")

cat(sprintf("mu_true = %.4f\n", mu_true))
cat("=== RUNNING MIXTURE MODEL (1000 reps) ===\n")
n_err <- 0

for (rep_i in 1:n_reps) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]
  tryCatch({
    ebmr <- EBMRAlgorithmFast4[["new"]]("y", single_ps_spec, dat, W_func)
    res <- ebmr[["EBMR_IPW"]](
      h_nu = function(dat) cbind(u1 = dat[["u1"]], u2 = dat[["u2"]], z1 = dat[["z1"]], z2 = dat[["z2"]]),
      se.fit = TRUE, type = "HT"
    )
    results[1, rep_i] <- res[["mu_ipw"]]
    results[3, rep_i] <- res[["se_ipw"]]
  }, error = function(e) {
    n_err <<- n_err + 1
    if (n_err <= 5) cat(sprintf("  Error rep %d: %s\n", rep_i, conditionMessage(e)))
  })
  if (rep_i %% 200 == 0) cat(sprintf("  %d / %d done (errors: %d)\n", rep_i, n_reps, n_err))
}

cat(sprintf("\nErrors: %d / %d\n", n_err, n_reps))

valid <- !is.na(results[1, ])
mu <- results[1, valid]
se <- results[3, valid]

bias <- mean(mu) - mu_true
esd <- sd(mu)
ese_mean <- mean(se)
ese_median <- median(se)

# Trim extreme SEs for robust mean
se_trimmed <- se[se < quantile(se, 0.99)]
ese_trimmed <- mean(se_trimmed)

ci_l <- mu - 1.96 * se
ci_u <- mu + 1.96 * se
cp <- mean(mu_true >= ci_l & mu_true <= ci_u)

cat(sprintf("\n=== RESULTS (model 2 alone, setting 3) ===\n"))
cat(sprintf("Valid reps: %d / %d\n", sum(valid), n_reps))
cat(sprintf("Bias:       %.4f\n", bias))
cat(sprintf("ESD:        %.4f\n", esd))
cat(sprintf("ESE mean:   %.4f\n", ese_mean))
cat(sprintf("ESE median: %.4f\n", ese_median))
cat(sprintf("ESE trim99: %.4f\n", ese_trimmed))
cat(sprintf("ESE/ESD (mean):   %.4f\n", ese_mean / esd))
cat(sprintf("ESE/ESD (median): %.4f\n", ese_median / esd))
cat(sprintf("ESE/ESD (trim):   %.4f\n", ese_trimmed / esd))
cat(sprintf("CP:         %.4f\n", cp))

# Check distribution of SE
cat(sprintf("\nSE quantiles: %.4f %.4f %.4f %.4f %.4f\n",
    quantile(se, 0.05), quantile(se, 0.25), quantile(se, 0.50),
    quantile(se, 0.75), quantile(se, 0.95)))
cat(sprintf("SE > 1: %d reps\n", sum(se > 1)))
cat(sprintf("SE > 10: %d reps\n", sum(se > 10)))

# Also check alpha norms
cat("\n=== Alpha norm distribution ===\n")
# Re-run a few to check
alpha_norms <- numeric(min(200, n_reps))
for (rep_i in 1:length(alpha_norms)) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]
  tryCatch({
    ebmr <- EBMRAlgorithmFast4[["new"]]("y", single_ps_spec, dat, W_func)
    alpha_norms[rep_i] <- sqrt(sum(ebmr[["ps_fit.list"]][[1]][["coefficients"]]^2))
  }, error = function(e) {})
}
cat(sprintf("||alpha|| quantiles: %.2f %.2f %.2f %.2f %.2f\n",
    quantile(alpha_norms, 0.05), quantile(alpha_norms, 0.25),
    quantile(alpha_norms, 0.50), quantile(alpha_norms, 0.75),
    quantile(alpha_norms, 0.95)))
