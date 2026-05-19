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

mu_vals <- se_vals <- obj_vals <- numeric(n_reps)
hess_pd <- logical(n_reps)
alpha_list <- vector("list", n_reps)
for (i in 1:n_reps) {
  set.seed(12345 + i)
  dat  <- all_data[((i-1)*n_val + 1):(i*n_val), ]
  ebmr <- EBMRAlgorithmFast4$new("y", ps_sub, dat, W_func)
  res  <- ebmr$EBMR_IPW(h_nu=h_nu_fn, type="HT", se.fit=TRUE)
  mu_vals[i] <- res$mu_ipw
  se_vals[i] <- res$se_ipw
  fit <- ebmr$ps_fit.list[[1]]$gmm_fit
  hess_pd[i] <- isTRUE(fit$opt$hessian_pd)
  obj_vals[i] <- fit$opt$objective
  alpha_list[[i]] <- fit$estimates
}

cat("=== All reps ===\n")
cat(sprintf("ESE/ESD: %.3f  (ESD=%.4f  ESE=%.4f)\n",
    mean(se_vals)/sd(mu_vals), sd(mu_vals), mean(se_vals)))

cat("\n=== hess_pd=TRUE (%d/%d) ===\n", sum(hess_pd), n_reps)
pd  <- hess_pd
cat(sprintf("ESE/ESD: %.3f  (ESD=%.4f  ESE=%.4f  Bias=%.4f)\n",
    mean(se_vals[pd])/sd(mu_vals[pd]),
    sd(mu_vals[pd]), mean(se_vals[pd]), mean(mu_vals[pd]) - mu_true))

cat(sprintf("\n=== hess_pd=FALSE (%d/%d) ===\n", sum(!pd), n_reps))
cat(sprintf("ESE/ESD: %.3f  (ESD=%.4f  ESE=%.4f  Bias=%.4f)\n",
    mean(se_vals[!pd])/sd(mu_vals[!pd]),
    sd(mu_vals[!pd]), mean(se_vals[!pd]), mean(mu_vals[!pd]) - mu_true))

cat("\n=== Objective values (CUE) ===\n")
cat(sprintf("PD:    mean=%.4f  min=%.4f  max=%.4f\n",
    mean(obj_vals[pd]), min(obj_vals[pd]), max(obj_vals[pd])))
cat(sprintf("Non-PD: mean=%.4f  min=%.4f  max=%.4f\n",
    mean(obj_vals[!pd]), min(obj_vals[!pd]), max(obj_vals[!pd])))

cat("\n=== Alpha estimates (a few non-PD cases) ===\n")
non_pd_idx <- which(!pd)[1:min(5, sum(!pd))]
for (i in non_pd_idx) {
  cat(sprintf("rep %d: alpha=[%s]  mu=%.3f  se=%.4f  obj=%.6f\n",
      i, paste(round(alpha_list[[i]], 3), collapse=", "),
      mu_vals[i], se_vals[i], obj_vals[i]))
}

cat("\n=== Alpha estimates (a few PD cases) ===\n")
pd_idx <- which(pd)[1:min(5, sum(pd))]
for (i in pd_idx) {
  cat(sprintf("rep %d: alpha=[%s]  mu=%.3f  se=%.4f  obj=%.6f\n",
      i, paste(round(alpha_list[[i]], 3), collapse=", "),
      mu_vals[i], se_vals[i], obj_vals[i]))
}
