## Root cause analysis: ESE/ESD discrepancy for mu_010, scenario 7-1, setting 3
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)
library(parallel); library(foreach); library(doSNOW)

mu_true <- get_mu_true("setting3")
ps_spec_base <- get_ps_spec("7")
h_alpha_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2)
W_sm <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_nu_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1u2=dat$u1*dat$u2)

data_file <- correct_model_all_data_file.list[["setting3"]][["miss50"]]
all_data <- NULL
for (fi in seq_along(data_file)) {
  if (file.exists(data_file[[fi]])) {
    test <- readRDS(data_file[[fi]])
    if (nrow(test) / 1000 == 2000) { all_data <- test; break }
  }
}

ps_spec_m2 <- list(
  formula.list = list(ps_spec_base$formula.list[[2]]),
  h_alpha.list = list(h_alpha_fn),
  inv_link = ps_spec_base$inv_link,
  outcome = ps_spec_base$outcome
)

nn <- 2000; n_reps <- 1000; n_cores <- min(detectCores() - 1, 10)

cl <- makeCluster(n_cores)
registerDoSNOW(cl)
clusterExport(cl, c("all_data", "nn", "ps_spec_m2", "W_sm", "h_nu_fn", "n_reps"),
              envir = environment())
clusterEvalQ(cl, devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE))

pb <- txtProgressBar(max = n_reps, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

sim_result <- foreach(
  i = 1:n_reps, .combine = 'rbind', .options.snow = opts,
  .packages = c("stringr", "Matrix", "numDeriv")
) %dopar% {
  tryCatch({
    dat <- all_data[((i-1)*nn + 1):(i*nn), ]
    ebmr <- EBMRAlgorithmFast4$new("y", ps_spec_m2, dat, W_sm)
    fit <- ebmr$ps_fit.list[[1]]
    ipw <- ebmr$EBMR_IPW(h_nu_fn, model_indices = 1, se.fit = TRUE)

    conv <- fit$gmm_fit$opt$converged
    grad <- fit$gmm_fit$opt$final_grad_norm
    ps_min <- min(fit$fitted.values)

    c(mu = ipw$mu_ipw, se = ipw$se_ipw, conv = as.numeric(conv),
      grad = grad, ps_min = ps_min)
  }, error = function(e) c(mu = NA, se = NA, conv = NA, grad = NA, ps_min = NA))
}

close(pb); stopCluster(cl)

mu_vals <- sim_result[,1]; se_vals <- sim_result[,2]
conv_vals <- sim_result[,3]; grad_vals <- sim_result[,4]
ps_min_vals <- sim_result[,5]

valid <- !is.na(mu_vals)

cat(sprintf("\n\nTrue mean: %.4f\n", mu_true))
cat(sprintf("Total: %d, NA: %d\n", n_reps, sum(!valid)))
cat(sprintf("Converged: %d/%d (%.1f%%)\n", sum(conv_vals[valid]==1), sum(valid),
    100*mean(conv_vals[valid]==1)))

# Split by convergence
conv <- conv_vals[valid] == 1
cat(sprintf("\n--- Converged (%d) ---\n", sum(conv)))
cat(sprintf("  Bias=%.4f, ESD=%.4f, ESE=%.4f, ESE/ESD=%.3f\n",
    mean(mu_vals[valid][conv]) - mu_true, sd(mu_vals[valid][conv]),
    mean(se_vals[valid][conv]), mean(se_vals[valid][conv])/sd(mu_vals[valid][conv])))

cat(sprintf("\n--- Non-converged (%d) ---\n", sum(!conv)))
cat(sprintf("  Bias=%.4f, ESD=%.4f, ESE=%.4f, ESE/ESD=%.3f\n",
    mean(mu_vals[valid][!conv]) - mu_true, sd(mu_vals[valid][!conv]),
    mean(se_vals[valid][!conv]), mean(se_vals[valid][!conv])/sd(mu_vals[valid][!conv])))

# Grad distribution for non-converged
cat(sprintf("\n  Non-converged grad: min=%.2e, med=%.2e, max=%.2e\n",
    min(grad_vals[valid][!conv]), median(grad_vals[valid][!conv]), max(grad_vals[valid][!conv])))

# PS min distribution
cat(sprintf("\nPS min distribution (all):\n"))
print(quantile(ps_min_vals[valid], c(0, 0.01, 0.05, 0.1, 0.25, 0.5)))

# Correlation between ps_min and extreme mu
cat(sprintf("\ncor(ps_min, |mu - mu_true|): %.4f\n",
    cor(ps_min_vals[valid], abs(mu_vals[valid] - mu_true))))

# Among converged, check if SE is correct
cat(sprintf("\n--- Converged, after IQR*1.5 cleaning ---\n"))
mu_c <- mu_vals[valid][conv]; se_c <- se_vals[valid][conv]
Q1 <- quantile(mu_c, 0.25); Q3 <- quantile(mu_c, 0.75)
iqr <- Q3 - Q1
keep <- mu_c >= (Q1 - 1.5*iqr) & mu_c <= (Q3 + 1.5*iqr)
cat(sprintf("  Outliers: %d/%d\n", sum(!keep), length(mu_c)))
cat(sprintf("  Bias=%.4f, ESD=%.4f, ESE=%.4f, ESE/ESD=%.3f\n",
    mean(mu_c[keep]) - mu_true, sd(mu_c[keep]),
    mean(se_c[keep]), mean(se_c[keep])/sd(mu_c[keep])))
