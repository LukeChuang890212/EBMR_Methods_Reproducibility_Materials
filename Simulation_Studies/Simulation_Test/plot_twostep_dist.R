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

cat("Running 1000 reps...\n")
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
  if (rep_i %% 200 == 0) cat(sprintf("  %d / %d done\n", rep_i, n_reps))
}

results[2, ] <- ifelse(is.na(results[1, ]), NA, 0)
results[4, ] <- ifelse(is.na(results[3, ]), NA, 0)

# Raw (no NA, no outlier removal)
valid <- !is.na(results[1, ])
mu_raw <- results[1, valid]
se_raw <- results[3, valid]

# After clean_sim_result with multiplier=1.5
cleaned <- clean_sim_result(results, multiplier = 1.5, verbose = TRUE)
sim <- cleaned$result
mu_clean <- sim[1, ]
se_clean <- sim[3, ]

# Plot
png("Simulation_Results/twostep_distribution.png", width = 1200, height = 800)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# mu_ipw before
hist(mu_raw, breaks = 50, main = "mu_ipw (before outlier removal)",
     xlab = "mu_ipw", col = "lightblue", border = "white")
abline(v = mu.true, col = "red", lwd = 2)
abline(v = mean(mu_raw), col = "blue", lwd = 2, lty = 2)
legend("topright", c(sprintf("mu.true = %.3f", mu.true),
                     sprintf("mean = %.3f", mean(mu_raw)),
                     sprintf("sd = %.3f", sd(mu_raw)),
                     sprintf("n = %d", length(mu_raw))),
       col = c("red", "blue", NA, NA), lwd = c(2, 2, NA, NA), lty = c(1, 2, NA, NA), cex = 0.9)

# mu_ipw after
hist(mu_clean, breaks = 50, main = "mu_ipw (after outlier removal, mult=1.5)",
     xlab = "mu_ipw", col = "lightgreen", border = "white")
abline(v = mu.true, col = "red", lwd = 2)
abline(v = mean(mu_clean), col = "blue", lwd = 2, lty = 2)
legend("topright", c(sprintf("mu.true = %.3f", mu.true),
                     sprintf("mean = %.3f", mean(mu_clean)),
                     sprintf("sd = %.3f", sd(mu_clean)),
                     sprintf("n = %d", length(mu_clean))),
       col = c("red", "blue", NA, NA), lwd = c(2, 2, NA, NA), lty = c(1, 2, NA, NA), cex = 0.9)

# se_ipw before
hist(se_raw, breaks = 50, main = "se_ipw (before outlier removal)",
     xlab = "se_ipw", col = "lightyellow", border = "white")
abline(v = mean(se_raw), col = "blue", lwd = 2, lty = 2)
legend("topright", c(sprintf("mean(se) = %.3f", mean(se_raw)),
                     sprintf("sd(mu) = %.3f", sd(mu_raw)),
                     sprintf("n = %d", length(se_raw))),
       col = c("blue", NA, NA), lwd = c(2, NA, NA), lty = c(2, NA, NA), cex = 0.9)

# se_ipw after
hist(se_clean, breaks = 50, main = "se_ipw (after outlier removal, mult=1.5)",
     xlab = "se_ipw", col = "mistyrose", border = "white")
abline(v = mean(se_clean), col = "blue", lwd = 2, lty = 2)
legend("topright", c(sprintf("ESE = %.3f", mean(se_clean)),
                     sprintf("ESD = %.3f", sd(mu_clean)),
                     sprintf("ESE/ESD = %.2f", mean(se_clean) / sd(mu_clean)),
                     sprintf("n = %d", length(se_clean))),
       col = c("blue", NA, NA, NA), lwd = c(2, NA, NA, NA), lty = c(2, NA, NA, NA), cex = 0.9)

dev.off()
cat("Plot saved to Simulation_Results/twostep_distribution.png\n")
