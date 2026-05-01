setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
source("Simulation.r")
library(EBMRalgorithmFast4)

# Scenario 9-3, model 2 only => ps_spec "9-alt1", model index 2
ps_spec <- get_ps_spec("9-alt1")
# "_2" = single model 2
ps_spec_2 <- list(
  formula.list = ps_spec[["formula.list"]][2],
  h_alpha.list = ps_spec[["h_alpha.list"]][2],
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]]
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
n <- 2000

# setting3, miss50 => setting3.B1
mu_true <- get_mu_true("setting3")

# Load pre-generated data
data_file <- "Simulation_Data/Setting3.B1_n2000_replicate1000.RDS"
if (file.exists(data_file)) {
  all_data <- readRDS(data_file)
  cat("Loaded pre-generated data\n")
} else {
  cat("Data file not found, generating on the fly\n")
  set.seed(999)
  all_data <- setting3.B1(n * 100)  # generate enough for 100 reps
}

# Run 100 reps to check
n_reps <- min(100, nrow(all_data) %/% n)
cat(sprintf("Running %d reps, mu_true=%.4f\n\n", n_reps, mu_true))

mu_vals <- numeric(n_reps)
se_vals <- numeric(n_reps)
alpha_norms <- numeric(n_reps)
hessian_pd <- logical(n_reps)
pi_extreme <- numeric(n_reps)

alpha.true <- misspecified_model_alpha.true.list$setting3$miss50$setting3_1_mild_1000
ps_model.true <- function(dat, alpha.true) {
  1 / (1 + exp(cbind(rep(1, nrow(dat)), dat$y, dat$u1, dat$u2) %*% alpha.true))
}

t0 <- Sys.time()
for (i in 1:n_reps) {
  dat <- all_data[((i-1)*n + 1):(i*n), ]

  ebmr <- EBMRAlgorithmFast4$new("y", ps_spec_2, dat, W_func)
  result <- ebmr$EBMR_IPW(
    h_nu = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2, u1_u2 = dat$u1*dat$u2),
    true_ps = ps_model.true(dat, alpha.true),
    type = "HT",
    se.fit = TRUE
  )

  mu_vals[i] <- result$mu_ipw
  se_vals[i] <- result$se_ipw
  ps_fit <- ebmr$ps_fit.list[[1]]
  alpha_norms[i] <- sqrt(sum(ps_fit$coefficients^2))
  hessian_pd[i] <- ps_fit$gmm_fit$opt$hessian_pd
  pi_hat <- ps_fit$fitted.values
  pi_extreme[i] <- sum(pi_hat < 0.05 | pi_hat > 0.95)

  if (i %% 20 == 0) cat(sprintf("Rep %d done\n", i))
}

elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("\nTime: %.1f sec (%.2f sec/rep)\n", elapsed, elapsed/n_reps))
cat("\n=== Results (new package) ===\n")
cat(sprintf("mu_true: %.4f\n", mu_true))
cat(sprintf("Mean mu: %.4f\n", mean(mu_vals)))
cat(sprintf("Bias: %.4f\n", mean(mu_vals) - mu_true))
cat(sprintf("ESD: %.4f\n", sd(mu_vals)))
cat(sprintf("ASE: %.4f\n", mean(se_vals)))
cat(sprintf("ASE/ESD: %.4f\n", mean(se_vals) / sd(mu_vals)))

# Outliers
q <- quantile(mu_vals, c(0.25, 0.75))
iqr <- q[2] - q[1]
outliers <- which(mu_vals < q[1] - 3*iqr | mu_vals > q[2] + 3*iqr)
cat(sprintf("Outliers (3*IQR): %d out of %d\n", length(outliers), n_reps))

if (length(outliers) > 0) {
  cat("Outlier mu values:", sort(round(mu_vals[outliers], 4)), "\n")
  cat("Outlier ||alpha||:", round(alpha_norms[outliers], 2), "\n")
}

cat(sprintf("\nWithout outliers: ESD=%.4f, ASE=%.4f, ratio=%.4f\n",
    sd(mu_vals[-outliers]), mean(se_vals[-outliers]),
    mean(se_vals[-outliers]) / sd(mu_vals[-outliers])))

cat(sprintf("\nhessian_pd: %d TRUE, %d FALSE\n", sum(hessian_pd), sum(!hessian_pd)))
cat(sprintf("Mean pi_extreme: %.1f\n", mean(pi_extreme)))
