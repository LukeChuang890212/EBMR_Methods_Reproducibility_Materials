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
all_data <- readRDS("Simulation_Data/Setting3.B1_n2000_replicate1000.RDS")
n_val <- 2000

outlier_reps <- c(125, 154, 180, 216, 221, 232, 255, 263, 285, 308, 314, 353, 420, 430, 438)

cat(sprintf("%-5s %-8s %-10s %-12s %-8s %-8s %s\n",
    "Rep", "hpd", "converged", "grad_norm", "iters", "mu", "alpha"))

for (i in outlier_reps) {
  set.seed(12345 + i)
  dat  <- all_data[((i-1)*n_val + 1):(i*n_val), ]
  ebmr <- EPS$new("y", ps_sub, dat, W_func)
  res  <- ebmr$EBMR_IPW(h_nu=h_nu_fn, type="HT", se.fit=FALSE)
  fit  <- ebmr$ps_fit.list[[1]]$gmm_fit
  opt  <- fit$opt

  cat(sprintf("%-5d %-8s %-10s %-12.2e %-8d %8.3f [%s]\n",
      i,
      isTRUE(opt$hessian_pd),
      isTRUE(opt$converged),
      ifelse(is.null(opt$final_grad_norm) || is.na(opt$final_grad_norm), NA, opt$final_grad_norm),
      ifelse(is.null(opt$iterations), NA_integer_, opt$iterations),
      res$mu_ipw,
      paste(round(fit$estimates, 2), collapse=", ")))
}
