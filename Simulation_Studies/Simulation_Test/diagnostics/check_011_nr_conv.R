## Check if mu_011's ESE/ESD=0.439 with constrained_nr is due to non-convergence
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EPS", quiet = TRUE)

ps_spec_base <- get_ps_spec("9-alt1")
h_alpha_fn <- function(dat) cbind(
  u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2,
  u1u2 = dat$u1*dat$u2, z1z2 = dat$z1*dat$z2,
  u1z1 = dat$u1*dat$z1, z1u2 = dat$z1*dat$u2,
  u1z2 = dat$u1*dat$z2, u2z2 = dat$u2*dat$z2)
ps_spec <- list(
  formula.list = ps_spec_base$formula.list,
  h_alpha.list = list(h_alpha_fn, h_alpha_fn, h_alpha_fn),
  inv_link = ps_spec_base$inv_link,
  outcome = ps_spec_base$outcome,
  optimizer = "constrained_nr"
)
W_fn <- function(g.matrix) solve(var(g.matrix))
h_nu_fn <- h_alpha_fn

data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]]
all_data <- NULL
for (fi in seq_along(data_file)) {
  if (file.exists(data_file[[fi]])) {
    test <- readRDS(data_file[[fi]])
    if (nrow(test) / 1000 == 2000) { all_data <- test; break }
  }
}

cat("=== Checking convergence of M2, M3, and nu for mu_011 with constrained_nr ===\n\n")

for (i in 1:10) {
  dat <- all_data[((i-1)*2000 + 1):(i*2000), ]
  ebmr <- EPS$new("y", ps_spec, dat, W_fn)
  ipw <- ebmr$EBMR_IPW(h_nu_fn, model_indices = 2:3, se.fit = TRUE)

  # Check convergence of M2 and M3
  conv2 <- ebmr$ps_fit.list[[2]]$gmm_fit$opt$converged
  conv3 <- ebmr$ps_fit.list[[3]]$gmm_fit$opt$converged
  grad2 <- ebmr$ps_fit.list[[2]]$gmm_fit$opt$final_grad_norm
  grad3 <- ebmr$ps_fit.list[[3]]$gmm_fit$opt$final_grad_norm

  # Nu convergence
  nu_conv <- ipw$ensemble_fit$gmm_fit$opt$converged
  nu_grad <- ipw$ensemble_fit$gmm_fit$opt$final_grad_norm

  cat(sprintf("Rep %2d: mu=%.4f, se=%.4f | M2: conv=%s grad=%.2e | M3: conv=%s grad=%.2e | nu: conv=%s grad=%.2e | w=(%.3f,%.3f)\n",
      i, ipw$mu_ipw, ipw$se_ipw,
      conv2, grad2, conv3, grad3,
      nu_conv, nu_grad,
      ipw$w.hat[1], ipw$w.hat[2]))
}
