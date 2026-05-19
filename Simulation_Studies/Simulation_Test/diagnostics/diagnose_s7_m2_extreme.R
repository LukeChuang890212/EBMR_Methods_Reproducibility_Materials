## Check extreme reps from the saved result
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

ps_spec_base <- get_ps_spec("7")
h_alpha_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2)
W_sm <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

# Load setting 4 data
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]]
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

nn <- 2000

# Check the extreme reps: 158, 207, 225
cat("=== Checking extreme reps ===\n\n")
for (i in c(158, 207, 225, 242, 267)) {
  if (i * nn > nrow(all_data)) { cat(sprintf("Rep %d: data not available\n\n", i)); next }
  dat <- all_data[((i-1)*nn + 1):(i*nn), ]
  ebmr <- EBMRAlgorithmFast4$new("y", ps_spec_m2, dat, W_sm)
  fit <- ebmr$ps_fit.list[[1]]

  alpha <- fit$coefficients
  ps <- fit$fitted.values
  conv <- fit$gmm_fit$opt$converged
  grad <- fit$gmm_fit$opt$final_grad_norm
  iter <- fit$gmm_fit$opt$iterations

  mu_ipw <- mean(dat$r * dat$y / ps)

  cat(sprintf("Rep %3d: conv=%5s grad=%.2e iter=%d\n", i, conv, grad, iter))
  cat(sprintf("         alpha=(%s)\n", paste(round(alpha, 4), collapse=", ")))
  cat(sprintf("         PS: min=%.6f Q1=%.4f med=%.4f max=%.4f | mu=%.4f\n\n",
      min(ps), quantile(ps, 0.25), median(ps), quantile(ps, 0.75), mu_ipw))
}

# Also check: how many reps have data available?
cat("Total data rows:", nrow(all_data), "\n")
cat("Max rep with n=2000:", nrow(all_data) / nn, "\n")
