setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
  devtools::load_all("../EBMRalgorithmFast4")
})

W_func  <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_nu_fn <- function(dat) cbind(u1=dat[["u1"]], u2=dat[["u2"]], z1=dat[["z1"]], z2=dat[["z2"]], u1_u2=dat[["u1"]]*dat[["u2"]])
n_val   <- 2000
n_reps  <- 200

ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
mu_true <- get_mu_true("setting4")

cat(sprintf("mu.true = %.4f\n\n", mu_true))

for (model_idx in 2:3) {
  ps_sub <- list(
    formula.list = ps_spec[["formula.list"]][model_idx],
    h_alpha.list = ps_spec[["h_alpha.list"]][model_idx],
    inv_link     = ps_spec[["inv_link"]],
    outcome      = ps_spec[["outcome"]]
  )

  cat(sprintf("=== Model %d only, %d reps ===\n", model_idx, n_reps))

  mu_vals  <- rep(NA_real_, n_reps)
  se1_vals <- rep(NA_real_, n_reps)

  n_err <- 0
  for (i in 1:n_reps) {
    dat <- all_data[((i-1)*n_val + 1):(i*n_val), ]
    tryCatch({
      ebmr <- EBMRAlgorithmFast4[["new"]]("y", ps_sub, dat, W_func)
      res  <- ebmr[["EBMR_IPW"]](h_nu=h_nu_fn, type="HT", se.fit=TRUE)
      mu_vals[i]  <- res[["mu_ipw"]]
      se1_vals[i] <- res[["se_ipw"]]
    }, error = function(e) {
      n_err <<- n_err + 1
      if (n_err <= 3) cat(sprintf("  Error rep %d: %s\n", i, conditionMessage(e)))
    })
    if (i %% 50 == 0) cat(sprintf("  %d / %d done (errors: %d)\n", i, n_reps, n_err))
  }

  valid <- !is.na(mu_vals)
  mu_v  <- mu_vals[valid]
  se1_v <- se1_vals[valid]
  esd   <- sd(mu_v)
  ese1  <- mean(se1_v)

  ci_l <- mu_v - 1.96 * se1_v; ci_u <- mu_v + 1.96 * se1_v
  cp   <- mean(mu_true >= ci_l & mu_true <= ci_u)

  cat(sprintf("\n  Errors: %d / %d, Valid: %d\n", n_err, n_reps, sum(valid)))
  cat(sprintf("  Bias:  %.4f\n", mean(mu_v) - mu_true))
  cat(sprintf("  ESD:   %.4f\n", esd))
  cat(sprintf("  ESE1:  %.4f  (ESE1/ESD = %.3f)\n", ese1, ese1/esd))
  cat(sprintf("  CP:    %.3f\n\n", cp))
}
