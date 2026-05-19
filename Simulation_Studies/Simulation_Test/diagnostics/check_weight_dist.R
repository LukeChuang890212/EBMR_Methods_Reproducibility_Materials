setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
  library(EPS)
})

W_func  <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_nu_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1_u2=dat$u1*dat$u2)

ps_spec <- get_ps_spec("9-alt1")
ps_sub  <- list(formula.list = ps_spec[["formula.list"]][2],
                h_alpha.list = ps_spec[["h_alpha.list"]][2],
                inv_link     = ps_spec[["inv_link"]],
                outcome      = ps_spec[["outcome"]])
mu_true <- get_mu_true("setting3")
all_data <- readRDS("Simulation_Data/Setting3.B1_n2000_replicate1000.RDS")
n_val <- 2000
n_reps <- 500

mu_vals  <- numeric(n_reps)
w_mean   <- w_max <- w_sd <- w_p99 <- numeric(n_reps)

for (i in 1:n_reps) {
  set.seed(12345 + i)
  dat  <- all_data[((i-1)*n_val + 1):(i*n_val), ]
  ebmr <- EPS$new("y", ps_sub, dat, W_func)
  res  <- ebmr$EBMR_IPW(h_nu=h_nu_fn, type="HT", se.fit=FALSE)
  mu_vals[i] <- res$mu_ipw

  pi_hat <- ebmr$ps_fit.list[[1]]$fitted.values
  r      <- dat$r
  wt     <- r / pi_hat   # IPW weights (only nonzero for R=1)
  wt_obs <- wt[r == 1]

  w_mean[i] <- mean(wt_obs)
  w_sd[i]   <- sd(wt_obs)
  w_max[i]  <- max(wt_obs)
  w_p99[i]  <- quantile(wt_obs, 0.99)
}

resid  <- mu_vals - mu_true
sd_mu  <- sd(mu_vals)
is_out <- abs(resid) > 3 * sd_mu

cat(sprintf("Outliers: %d/%d\n\n", sum(is_out), n_reps))

cat(sprintf("%-10s  %8s  %8s  %8s  %8s\n", "Group", "mean(w)", "sd(w)", "max(w)", "p99(w)"))
cat(sprintf("%-10s  %8.3f  %8.3f  %8.1f  %8.3f\n", "Normal",
    mean(w_mean[!is_out]), mean(w_sd[!is_out]),
    mean(w_max[!is_out]),  mean(w_p99[!is_out])))
cat(sprintf("%-10s  %8.3f  %8.3f  %8.1f  %8.3f\n", "Outlier",
    mean(w_mean[is_out]),  mean(w_sd[is_out]),
    mean(w_max[is_out]),   mean(w_p99[is_out])))

cat("\n=== Per outlier rep ===\n")
cat(sprintf("%-5s  %8s  %8s  %8s  %8s  %8s\n", "Rep", "mu", "mean(w)", "sd(w)", "max(w)", "p99(w)"))
for (i in which(is_out)) {
  cat(sprintf("%-5d  %8.3f  %8.3f  %8.3f  %8.1f  %8.3f\n",
      i, mu_vals[i], w_mean[i], w_sd[i], w_max[i], w_p99[i]))
}
