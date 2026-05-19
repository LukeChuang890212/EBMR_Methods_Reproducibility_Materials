## Diagnose mu_110 under miss30: CP=0.897 anomaly
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EIPS", quiet = TRUE)
library(parallel); library(foreach); library(doSNOW)

mu_true <- get_mu_true("setting4")
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

cl <- makeCluster(n_cores)
registerDoSNOW(cl)
clusterExport(cl, c("nn", "ps_spec_sub", "W_sm", "h_nu_fn", "n_reps"), envir = environment())
clusterEvalQ(cl, {
  setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
  devtools::load_all("../EIPS", quiet = TRUE)
  source("Data_Generation.r")
})
pb <- txtProgressBar(max = n_reps, style = 3)
opts <- list(progress = function(n) setTxtProgressBar(pb, n))

sim_result <- foreach(
  i = 1:n_reps, .combine = 'rbind', .options.snow = opts,
  .packages = c("stringr", "Matrix")
) %dopar% {
  tryCatch({
    set.seed(i)
    dat <- setting4.A2(nn)
    ebmr <- EIPS$new("y", ps_spec_sub, dat, W_sm)
    ipw <- ebmr$EBMR_IPW(h_nu_fn, model_indices = 1:2, se.fit = TRUE)
    c(mu = ipw$mu_ipw, se = ipw$se_ipw, w1 = ipw$w.hat[1], w2 = ipw$w.hat[2])
  }, error = function(e) c(mu = NA, se = NA, w1 = NA, w2 = NA))
}
close(pb); stopCluster(cl)

mu_v <- sim_result[,1]; se_v <- sim_result[,2]
w1_v <- sim_result[,3]; w2_v <- sim_result[,4]
valid <- !is.na(mu_v) & !is.na(se_v)
mu_v <- mu_v[valid]; se_v <- se_v[valid]; w1_v <- w1_v[valid]; w2_v <- w2_v[valid]

cat(sprintf("\n\nmu_true = %.4f\n", mu_true))
cat(sprintf("mu:  mean=%.4f sd=%.4f min=%.4f med=%.4f max=%.4f\n",
            mean(mu_v), sd(mu_v), min(mu_v), median(mu_v), max(mu_v)))
cat(sprintf("se:  mean=%.4f sd=%.4f min=%.4f med=%.4f max=%.4f\n",
            mean(se_v), sd(se_v), min(se_v), median(se_v), max(se_v)))
cat(sprintf("w1 (M1): mean=%.4f sd=%.4f min=%.4f med=%.4f max=%.4f\n",
            mean(w1_v), sd(w1_v), min(w1_v), median(w1_v), max(w1_v)))
cat(sprintf("w2 (M2): mean=%.4f sd=%.4f min=%.4f med=%.4f max=%.4f\n",
            mean(w2_v), sd(w2_v), min(w2_v), median(w2_v), max(w2_v)))

# Coverage
ci_lo <- mu_v - 1.96*se_v; ci_hi <- mu_v + 1.96*se_v
covered <- ci_lo <= mu_true & mu_true <= ci_hi
cat(sprintf("\nRaw CP = %.3f (%d/%d)\n", mean(covered), sum(covered), length(covered)))

# Where do non-covered reps fall?
nc <- !covered
cat(sprintf("\nNon-covered reps (%d):\n", sum(nc)))
cat(sprintf("  mu:  mean=%.4f sd=%.4f range=[%.4f, %.4f]\n",
            mean(mu_v[nc]), sd(mu_v[nc]), min(mu_v[nc]), max(mu_v[nc])))
cat(sprintf("  se:  mean=%.4f (vs covered mean se=%.4f)\n",
            mean(se_v[nc]), mean(se_v[!nc])))
cat(sprintf("  w1:  mean=%.4f (vs covered mean w1=%.4f)\n",
            mean(w1_v[nc]), mean(w1_v[!nc])))
cat(sprintf("  |mu-mu_true|: mean=%.4f (vs covered=%.4f)\n",
            mean(abs(mu_v[nc]-mu_true)), mean(abs(mu_v[!nc]-mu_true))))

# Is there bimodality in w1?
cat(sprintf("\nw1 distribution: <0.5: %d, 0.5-0.9: %d, >0.9: %d\n",
            sum(w1_v < 0.5), sum(w1_v >= 0.5 & w1_v < 0.9), sum(w1_v >= 0.9)))

# Among reps where M2 gets notable weight (w2 > 0.1)
hi_w2 <- w2_v > 0.1
if (sum(hi_w2) > 0) {
  cat(sprintf("\nReps with w2 > 0.1 (%d): mu mean=%.4f, CP=%.3f\n",
              sum(hi_w2), mean(mu_v[hi_w2]), mean(covered[hi_w2])))
}
cat(sprintf("Reps with w2 <= 0.1 (%d): mu mean=%.4f, CP=%.3f\n",
            sum(!hi_w2), mean(mu_v[!hi_w2]), mean(covered[!hi_w2])))

# se vs |error| correlation
cat(sprintf("\ncor(se, |mu-mu_true|) = %.3f\n", cor(se_v, abs(mu_v - mu_true))))
