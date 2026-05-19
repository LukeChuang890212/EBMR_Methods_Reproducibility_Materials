## Diagnose M2 alpha convergence for scenario 7-1, setting 4
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EPS", quiet = TRUE)

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
  formula.list = list(ps_spec_base$formula.list[[2]]),  # M2: r ~ y + u1 + z1
  h_alpha.list = list(h_alpha_fn),
  inv_link = ps_spec_base$inv_link,
  outcome = ps_spec_base$outcome
)

cat("M2 formula:", deparse(ps_spec_base$formula.list[[2]]), "\n\n")

nn <- 2000
cat("=== M2 alpha convergence check, 10 reps ===\n\n")
for (i in 1:10) {
  dat <- all_data[((i-1)*nn + 1):(i*nn), ]
  ebmr <- EPS$new("y", ps_spec_m2, dat, W_sm)
  fit <- ebmr$ps_fit.list[[1]]

  alpha <- fit$coefficients
  ps <- fit$fitted.values
  conv <- fit$gmm_fit$opt$converged
  grad <- fit$gmm_fit$opt$final_grad_norm

  mu_ipw <- mean(dat$r * dat$y / ps)

  cat(sprintf("Rep %2d: conv=%5s grad=%.2e | alpha=(%s)\n",
      i, conv, grad, paste(round(alpha, 4), collapse=", ")))
  cat(sprintf("        PS: min=%.4f Q1=%.4f med=%.4f Q3=%.4f max=%.4f | mu=%.4f\n\n",
      min(ps), quantile(ps, 0.25), median(ps), quantile(ps, 0.75), max(ps), mu_ipw))
}
