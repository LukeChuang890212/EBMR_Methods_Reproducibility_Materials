## Check if PS collapses to 0.5 with cond=1e2 for scenario 7-1, setting 3
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
    if (nrow(test) / 1000 >= 2) { all_data <- test; break }
  }
}

nn <- 2000
ps_spec <- list(
  formula.list = list(ps_spec_base$formula.list[[2]]),
  h_alpha.list = list(h_alpha_fn),
  inv_link = ps_spec_base$inv_link,
  outcome = ps_spec_base$outcome,
  optimizer = "constrained_nr",
  cond_threshold = 1e2
)

cat("Setting 3, M2 (r ~ y + u1 + z1), NR cond=1e2, 10 reps:\n\n")
for (i in 1:10) {
  dat <- all_data[((i-1)*nn + 1):(i*nn), ]
  ebmr <- EBMRAlgorithmFast4$new("y", ps_spec, dat, W_sm)
  ps <- ebmr$ps_fit.list[[1]]$fitted.values
  alpha <- ebmr$ps_fit.list[[1]]$coefficients
  cat(sprintf("  Rep %2d: PS min=%.4f med=%.4f max=%.4f sd=%.6f | alpha=(%s)\n",
      i, min(ps), median(ps), max(ps), sd(ps),
      paste(round(alpha, 4), collapse=",")))
}
