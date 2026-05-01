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
ps_sub <- list(
  formula.list = ps_spec[["formula.list"]][2:3],
  h_alpha.list = ps_spec[["h_alpha.list"]][2:3],
  inv_link     = ps_spec[["inv_link"]],
  outcome      = ps_spec[["outcome"]]
)

data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
mu_true <- get_mu_true("setting4")

cat(sprintf("mu.true = %.4f\n", mu_true))
cat(sprintf("=== Newton-Raphson ensemble (no bounds), %d reps ===\n\n", n_reps))

mu_vals  <- rep(NA_real_, n_reps)
se1_vals <- rep(NA_real_, n_reps)
w1_vals  <- rep(NA_real_, n_reps)
nu1_vals <- rep(NA_real_, n_reps)
nu2_vals <- rep(NA_real_, n_reps)

n_err <- 0
for (i in 1:n_reps) {
  dat <- all_data[((i-1)*n_val + 1):(i*n_val), ]
  tryCatch({
    ebmr <- EBMRAlgorithmFast4[["new"]]("y", ps_sub, dat, W_func)
    res  <- ebmr[["EBMR_IPW"]](h_nu=h_nu_fn, type="HT", se.fit=TRUE)
    mu_vals[i]  <- res[["mu_ipw"]]
    se1_vals[i] <- res[["se_ipw"]]
    w_hat <- res[["w.hat"]]
    w1_vals[i]  <- w_hat[1]
    nu_hat <- res[["nu.hat"]]
    nu1_vals[i] <- nu_hat[1]
    nu2_vals[i] <- nu_hat[2]
  }, error = function(e) {
    n_err <<- n_err + 1
    if (n_err <= 5) cat(sprintf("  Error rep %d: %s\n", i, conditionMessage(e)))
  })
  if (i %% 50 == 0) cat(sprintf("  %d / %d done (errors: %d)\n", i, n_reps, n_err))
}

cat(sprintf("\nErrors: %d / %d\n", n_err, n_reps))

valid <- !is.na(mu_vals)
mu_v  <- mu_vals[valid]
se1_v <- se1_vals[valid]
w1_v  <- w1_vals[valid]

esd  <- sd(mu_v)
ese1 <- mean(se1_v)

cat(sprintf("\n=== RESULTS (Newton-Raphson ensemble, models 2+3, scenario 9-3) ===\n"))
cat(sprintf("Valid reps: %d\n", sum(valid)))
cat(sprintf("Bias: %.4f\n", mean(mu_v) - mu_true))
cat(sprintf("ESD:  %.4f\n", esd))
cat(sprintf("ESE1: %.4f  (ESE1/ESD = %.3f)\n", ese1, ese1/esd))

ci1_l <- mu_v - 1.96 * se1_v; ci1_u <- mu_v + 1.96 * se1_v
cp1   <- mean(mu_true >= ci1_l & mu_true <= ci1_u)
cat(sprintf("CP:   %.3f\n", cp1))

cat(sprintf("\nw1 (model 1 weight): mean=%.3f, sd=%.3f, range=[%.3f, %.3f]\n",
    mean(w1_v), sd(w1_v), min(w1_v), max(w1_v)))
cat(sprintf("nu1: mean=%.3f, sd=%.3f, range=[%.3f, %.3f]\n",
    mean(nu1_vals[valid]), sd(nu1_vals[valid]), min(nu1_vals[valid]), max(nu1_vals[valid])))
cat(sprintf("nu2: mean=%.3f, sd=%.3f, range=[%.3f, %.3f]\n",
    mean(nu2_vals[valid]), sd(nu2_vals[valid]), min(nu2_vals[valid]), max(nu2_vals[valid])))

# Distribution of w1
cat(sprintf("\nw1 distribution:\n"))
cat(sprintf("  w1 < 0.1: %d reps\n", sum(w1_v < 0.1)))
cat(sprintf("  w1 in [0.1, 0.9]: %d reps\n", sum(w1_v >= 0.1 & w1_v <= 0.9)))
cat(sprintf("  w1 > 0.9: %d reps\n", sum(w1_v > 0.9)))
