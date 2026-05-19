setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
  library(EPS)
})

W_func  <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_nu_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1_u2=dat$u1*dat$u2)

# Models 2 and 3 from scenario 9-3 (ps_spec = 9-alt1)
ps_spec <- get_ps_spec("9-alt1")
ps_sub  <- list(formula.list = ps_spec[["formula.list"]][2:3],
                h_alpha.list = ps_spec[["h_alpha.list"]][2:3],
                inv_link     = ps_spec[["inv_link"]],
                outcome      = ps_spec[["outcome"]])
mu_true <- 0.75
all_data <- readRDS("Simulation_Data/setting4.B1_n2000_replicate1000.RDS")
n_val <- 2000
n_reps <- 500

mu_vals <- se_vals <- nu1_vals <- numeric(n_reps)
hpd_nu  <- logical(n_reps)

for (i in 1:n_reps) {
  set.seed(12345 + i)
  dat  <- all_data[((i-1)*n_val + 1):(i*n_val), ]
  ebmr <- EPS$new("y", ps_sub, dat, W_func)
  res  <- ebmr$EBMR_IPW(h_nu=h_nu_fn, type="HT", se.fit=TRUE)
  mu_vals[i]  <- res$mu_ipw
  se_vals[i]  <- res$se_ipw
  nu1_vals[i] <- res$nu.hat[1]
  # Check nu hessian_pd if available
  nu_fit <- tryCatch(ebmr$nu_fit, error=function(e) NULL)
  if (!is.null(nu_fit) && !is.null(nu_fit$gmm_fit)) {
    hpd_nu[i] <- isTRUE(nu_fit$gmm_fit$opt$hessian_pd)
  } else {
    hpd_nu[i] <- NA
  }
}

# Detect outliers using same IQR logic as clean_sim_result
sim_result <- rbind(mu_ipw      = mu_vals,
                    mu_ipw.true = rep(mu_true, n_reps),
                    se_ipw      = se_vals,
                    se_ipw.true = se_vals)
colnames(sim_result) <- as.character(seq_len(n_reps))
cleaned  <- clean_sim_result(sim_result, multiplier = 3, max_pct = 0.01, verbose = FALSE)
kept_idx <- as.integer(colnames(cleaned$result))
is_out   <- !(seq_len(n_reps) %in% kept_idx)
norm_idx <- !is_out

cat(sprintf("n_reps=%d  mu_true=%.4f  ESD=%.4f  ESE=%.4f  ESE/ESD=%.3f\n",
    n_reps, mu_true, sd(mu_vals), mean(se_vals), mean(se_vals)/sd(mu_vals)))
cat(sprintf("Outliers (IQR method, cap 1%%): %d/%d\n\n", sum(is_out), n_reps))
cat(sprintf("Excluding outliers: ESD=%.4f  ESE=%.4f  ESE/ESD=%.3f  (n=%d)\n\n",
    sd(mu_vals[norm_idx]), mean(se_vals[norm_idx]),
    mean(se_vals[norm_idx])/sd(mu_vals[norm_idx]), sum(norm_idx)))

# Nu switching analysis
cat("=== Nu (ensemble weight for model 2) distribution ===\n")
cat(sprintf("nu1: mean=%.3f  sd=%.3f  min=%.3f  max=%.3f\n",
    mean(nu1_vals), sd(nu1_vals), min(nu1_vals), max(nu1_vals)))

m2_dom <- nu1_vals > 0.9
m3_dom <- nu1_vals < 0.1
mixed  <- nu1_vals >= 0.1 & nu1_vals <= 0.9

cat(sprintf("nu1 > 0.9 (model2 dom): %d/%d (%.1f%%)  mean_mu=%.4f\n",
    sum(m2_dom), n_reps, 100*mean(m2_dom), mean(mu_vals[m2_dom])))
cat(sprintf("nu1 < 0.1 (model3 dom): %d/%d (%.1f%%)  mean_mu=%.4f\n",
    sum(m3_dom), n_reps, 100*mean(m3_dom), mean(mu_vals[m3_dom])))
cat(sprintf("0.1<=nu1<=0.9 (mixed):  %d/%d (%.1f%%)  mean_mu=%.4f\n",
    sum(mixed),  n_reps, 100*mean(mixed),  mean(mu_vals[mixed])))
