setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
  library(EBMRalgorithmFast4)
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
n_reps <- 200

mu_vals <- se_vals <- numeric(n_reps)
hess_pd <- logical(n_reps)
for (i in 1:n_reps) {
  set.seed(12345 + i)
  dat  <- all_data[((i-1)*n_val + 1):(i*n_val), ]
  ebmr <- EBMRAlgorithmFast4$new("y", ps_sub, dat, W_func)
  res  <- ebmr$EBMR_IPW(h_nu=h_nu_fn, type="HT", se.fit=TRUE)
  mu_vals[i] <- res$mu_ipw
  se_vals[i] <- res$se_ipw
  hess_pd[i] <- isTRUE(ebmr$ps_fit.list[[1]]$gmm_fit$opt$hessian_pd)
}

cat(sprintf("setting3-miss50-9-3_2 (M2 only), n=2000, %d reps\n", n_reps))
cat(sprintf("mu_true = %.4f\n", mu_true))
cat(sprintf("Bias:    %.4f\n", mean(mu_vals) - mu_true))
cat(sprintf("ESD:     %.4f\n", sd(mu_vals)))
cat(sprintf("ESE:     %.4f\n", mean(se_vals)))
cat(sprintf("ESE/ESD: %.3f\n", mean(se_vals)/sd(mu_vals)))
cat(sprintf("hess_pd: %d/%d\n", sum(hess_pd), n_reps))
