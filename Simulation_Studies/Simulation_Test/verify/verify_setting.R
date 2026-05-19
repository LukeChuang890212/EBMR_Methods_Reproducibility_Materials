## Verify: compare my NR(1e8) vs package NR(1e8) for 20 reps
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

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

nn <- 2000

ps_spec_pkg <- list(
  formula.list = list(ps_spec_base$formula.list[[2]]),
  h_alpha.list = list(h_alpha_fn),
  inv_link = ps_spec_base$inv_link,
  outcome = ps_spec_base$outcome,
  optimizer = "constrained_nr"
)

cat("Verify: package NR vs my test script, 20 reps\n\n")
cat(sprintf("%-5s | %-20s | %-20s\n", "Rep", "Package NR(1e8)", "Package L-BFGS-B+fb"))
cat(strrep("-", 50), "\n")

ps_spec_lb <- list(
  formula.list = list(ps_spec_base$formula.list[[2]]),
  h_alpha.list = list(h_alpha_fn),
  inv_link = ps_spec_base$inv_link,
  outcome = ps_spec_base$outcome
)

for (i in 1:20) {
  dat <- all_data[((i-1)*nn + 1):(i*nn), ]

  ebmr_nr <- EBMRAlgorithmFast4$new("y", ps_spec_pkg, dat, W_sm)
  ipw_nr <- ebmr_nr$EBMR_IPW(h_nu_fn, model_indices = 1, se.fit = TRUE)

  ebmr_lb <- EBMRAlgorithmFast4$new("y", ps_spec_lb, dat, W_sm)
  ipw_lb <- ebmr_lb$EBMR_IPW(h_nu_fn, model_indices = 1, se.fit = TRUE)

  cat(sprintf("%5d | mu=%.2f se=%.4f | mu=%.2f se=%.4f\n",
      i, ipw_nr$mu_ipw, ipw_nr$se_ipw, ipw_lb$mu_ipw, ipw_lb$se_ipw))
}
