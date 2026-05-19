## Check PS distribution for scenario 9-2, setting 4, miss50, M3, NR cond=1e2
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EPS", quiet = TRUE)

ps_spec_base <- get_ps_spec("9")
h_alpha_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2)
W_sm <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]]
all_data <- NULL
for (fi in seq_along(data_file)) {
  if (file.exists(data_file[[fi]])) {
    test <- readRDS(data_file[[fi]])
    if (nrow(test) / 1000 >= 2) { all_data <- test; break }
  }
}

nn <- 2000

ps_spec_m3 <- list(
  formula.list = list(ps_spec_base$formula.list[[3]]),
  h_alpha.list = list(h_alpha_fn),
  inv_link = ps_spec_base$inv_link,
  outcome = ps_spec_base$outcome,
  optimizer = "constrained_nr",
  cond_threshold = 1e2
)

cat("M3:", deparse(ps_spec_base$formula.list[[3]]), "\n\n")
for (cond_val in c(1e4, 1e6)) {
  ps_spec_m3 <- list(
    formula.list = list(ps_spec_base$formula.list[[3]]),
    h_alpha.list = list(h_alpha_fn),
    inv_link = ps_spec_base$inv_link,
    outcome = ps_spec_base$outcome,
    optimizer = "constrained_nr",
    cond_threshold = cond_val
  )

  cat(sprintf("=== NR cond=%.0e, 5 reps ===\n", cond_val))
  for (i in 1:5) {
    dat <- all_data[((i-1)*nn + 1):(i*nn), ]
    ebmr <- EPS$new("y", ps_spec_m3, dat, W_sm)
    ps <- ebmr$ps_fit.list[[1]]$fitted.values
    alpha <- ebmr$ps_fit.list[[1]]$coefficients
    mu <- mean(dat$r * dat$y / ps)
    cat(sprintf("  Rep %d: PS min=%.4f med=%.4f max=%.4f sd=%.4f mu=%.4f | alpha=(%s)\n",
        i, min(ps), median(ps), max(ps), sd(ps), mu,
        paste(round(alpha, 4), collapse=",")))
  }
  cat("\n")
}
