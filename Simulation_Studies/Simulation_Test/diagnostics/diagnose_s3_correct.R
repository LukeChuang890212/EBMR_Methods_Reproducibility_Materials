## Diagnose mu_010 with correct setting 3 data
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

ps_spec_base <- get_ps_spec("7")
h_alpha_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2)
W_sm <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

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

nn <- 2000
cat("M2:", deparse(ps_spec_base$formula.list[[2]]), "\n\n")

# Check 10 reps including some that might be extreme
cat("=== Checking reps ===\n\n")
for (i in 1:20) {
  dat <- all_data[((i-1)*nn + 1):(i*nn), ]
  ebmr <- EBMRAlgorithmFast4$new("y", ps_spec_m2, dat, W_sm)
  fit <- ebmr$ps_fit.list[[1]]

  alpha <- fit$coefficients
  ps <- fit$fitted.values
  conv <- fit$gmm_fit$opt$converged
  grad <- fit$gmm_fit$opt$final_grad_norm

  mu_ipw <- mean(dat$r * dat$y / ps)
  se_ipw <- tryCatch({
    ipw <- ebmr$EBMR_IPW(function(dat) cbind(u1=dat$u1,u2=dat$u2,z1=dat$z1,z2=dat$z2,u1u2=dat$u1*dat$u2),
                          model_indices = 1, se.fit = TRUE)
    ipw$se_ipw
  }, error = function(e) NA)

  extreme <- if (abs(mu_ipw) > 3) " *** EXTREME ***" else ""
  cat(sprintf("Rep %2d: conv=%s grad=%.2e | mu=%.4f se=%.4f | PS min=%.4f max=%.4f | alpha=(%s)%s\n",
      i, conv, grad, mu_ipw, se_ipw,
      min(ps), max(ps), paste(round(alpha, 3), collapse=","), extreme))
}
