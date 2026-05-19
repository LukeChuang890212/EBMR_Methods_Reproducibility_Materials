## Diagnose setting2 mu_110 miss30: CP=0.897 with Bias~0, ESE/ESD~0.99
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast5", quiet = TRUE)
library(parallel); library(foreach); library(doSNOW)

mu_true <- get_mu_true("setting2")
ps_spec_base <- get_ps_spec("9")
h_alpha_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2)
W_sm <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_nu_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1_u2=dat$u1*dat$u2)

nn <- 2000; n_reps <- 1000; n_cores <- min(detectCores() - 1, 10)
ps_spec_sub <- list(
  formula.list = ps_spec_base$formula.list[c(1,2)],
  h_alpha.list = list(h_alpha_fn, h_alpha_fn),
  inv_link = ps_spec_base$inv_link, outcome = ps_spec_base$outcome
)

cat(sprintf("mu_true (setting2) = %.6f\n", mu_true))

cl <- makeCluster(n_cores); registerDoSNOW(cl)
clusterExport(cl, c("nn","ps_spec_sub","W_sm","h_nu_fn","n_reps"), envir = environment())
clusterEvalQ(cl, {
  setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
  devtools::load_all("../EBMRalgorithmFast5", quiet = TRUE)
  source("Data_Generation.r")
})
pb <- txtProgressBar(max = n_reps, style = 3)
opts <- list(progress = function(n) setTxtProgressBar(pb, n))

sim_result <- foreach(
  i = 1:n_reps, .combine = 'rbind', .options.snow = opts,
  .packages = c("stringr","Matrix")
) %dopar% {
  tryCatch({
    set.seed(i)
    dat <- setting2.A2(nn)
    ebmr <- EBMRAlgorithmFast5$new("y", ps_spec_sub, dat, W_sm)
    ipw <- ebmr$EBMR_IPW(h_nu_fn, model_indices = 1:2, se.fit = TRUE)
    a1 <- ebmr$ps_fit.list[[1]]$coefficients
    ps1 <- ebmr$ps_fit.list[[1]]$fitted.values
    c(mu = ipw$mu_ipw, se = ipw$se_ipw, w1 = ipw$w.hat[1],
      ps1_min = min(ps1), ps1_max = max(ps1),
      a1_int = a1[1], a1_y = a1[2])
  }, error = function(e) rep(NA, 7))
}
close(pb); stopCluster(cl)

cn <- c("mu","se","w1","ps1_min","ps1_max","a1_int","a1_y")
colnames(sim_result) <- cn
valid <- !is.na(sim_result[,"mu"]) & !is.na(sim_result[,"se"])
sr <- sim_result[valid,]

mu_v <- sr[,"mu"]; se_v <- sr[,"se"]; w1_v <- sr[,"w1"]
cat(sprintf("\n\nn valid: %d\n", nrow(sr)))
cat(sprintf("mu:  mean=%.4f sd=%.4f  (mu_true=%.4f, bias=%.4f)\n",
            mean(mu_v), sd(mu_v), mu_true, mean(mu_v)-mu_true))
cat(sprintf("se:  mean=%.4f sd=%.4f min=%.4f max=%.4f\n",
            mean(se_v), sd(se_v), min(se_v), max(se_v)))
cat(sprintf("w1:  mean=%.4f sd=%.4f min=%.4f max=%.4f\n",
            mean(w1_v), sd(w1_v), min(w1_v), max(w1_v)))

# Coverage (raw)
covered <- (mu_v - 1.96*se_v <= mu_true) & (mu_true <= mu_v + 1.96*se_v)
cat(sprintf("\nRaw CP = %.3f\n", mean(covered)))

# Standardized residual
z <- (mu_v - mu_true) / se_v
cat(sprintf("z=(mu-mu_true)/se:  mean=%.4f sd=%.4f  (ideal: 0, 1)\n", mean(z), sd(z)))
cat(sprintf("  |z|>1.96: %.3f (ideal 0.05)\n", mean(abs(z) > 1.96)))
cat(sprintf("  z quantiles: %s\n",
            paste(sprintf("%.2f", quantile(z, c(.01,.05,.25,.5,.75,.95,.99))), collapse=" ")))

# Is se mismatched? Compare mean(se) vs sd(mu)
cat(sprintf("\nmean(se)=%.4f vs sd(mu)=%.4f  ratio=%.3f\n",
            mean(se_v), sd(mu_v), mean(se_v)/sd(mu_v)))

# Non-covered reps
nc <- !covered
cat(sprintf("\nNon-covered (%d reps):\n", sum(nc)))
cat(sprintf("  mu-mu_true: mean=%.4f range=[%.4f,%.4f]\n",
            mean(mu_v[nc]-mu_true), min(mu_v[nc]-mu_true), max(mu_v[nc]-mu_true)))
cat(sprintf("  se: mean=%.4f (covered: %.4f)\n", mean(se_v[nc]), mean(se_v[!nc])))
cat(sprintf("  w1: mean=%.4f (covered: %.4f)\n", mean(w1_v[nc]), mean(w1_v[!nc])))
cat(sprintf("  above mu_true: %d, below: %d\n",
            sum(mu_v[nc] > mu_true), sum(mu_v[nc] < mu_true)))

# se vs |error|
cat(sprintf("\ncor(se, |mu-mu_true|) = %.3f\n", cor(se_v, abs(mu_v-mu_true))))
cat(sprintf("PS1 fitted range: min med=%.4f, max med=%.4f\n",
            median(sr[,"ps1_min"]), median(sr[,"ps1_max"])))

# Skewness of mu
m3 <- mean((mu_v - mean(mu_v))^3) / sd(mu_v)^3
cat(sprintf("skewness(mu) = %.3f\n", m3))
