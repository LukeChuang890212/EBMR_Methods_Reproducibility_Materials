setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
source("Simulation.r")
library(EBMRalgorithmFast4)
library(parallel)
library(foreach)
library(doSNOW)

mu_true <- get_mu_true("setting4")

ps_spec <- get_ps_spec("9-alt1")
ps_spec_13 <- list(
  formula.list = ps_spec[["formula.list"]][c(1, 3)],
  h_alpha.list = ps_spec[["h_alpha.list"]][c(1, 3)],
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]]
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
n_val <- 2000
n_reps <- 500

cores <- max(1, detectCores() - 2)
cl <- makeCluster(cores)
registerDoSNOW(cl)
pb <- txtProgressBar(max = n_reps, style = 3)
opts <- list(progress = function(i) setTxtProgressBar(pb, i))
clusterExport(cl, c("setting4.B1", "ps_spec_13", "W_func", "n_val"), envir = environment())

results_raw <- foreach(rep_i = 1:n_reps,
                       .combine = cbind,
                       .options.snow = opts,
                       .packages = c("EBMRalgorithmFast4", "numDeriv"),
                       .errorhandling = "pass") %dopar% {
  set.seed(12345 + rep_i)
  dat <- setting4.B1(n_val)
  tryCatch({
    ebmr <- EBMRAlgorithmFast4[["new"]]("y", ps_spec_13, dat, W_func)
    res <- ebmr[["EBMR_IPW"]](
      h_nu = function(dat) cbind(u1 = dat[["u1"]], u2 = dat[["u2"]], z1 = dat[["z1"]], z2 = dat[["z2"]], u1_u2 = dat[["u1"]]*dat[["u2"]]),
      se.fit = TRUE, type = "HT"
    )
    c(res[["mu_ipw"]], res[["se_ipw"]], res[["nu.hat"]], res[["w.hat"]])
  }, error = function(e) rep(NA, 6))
}

close(pb)
stopCluster(cl)

if (is.list(results_raw)) {
  results_raw <- do.call(cbind, lapply(results_raw, function(x) if (is.numeric(x) && length(x) == 6) x else rep(NA, 6)))
}
rownames(results_raw) <- c("mu_ipw", "se_ipw", "nu1", "nu2", "w1", "w2")

valid <- !is.na(results_raw[1, ])
r <- results_raw[, valid]
cat(sprintf("Valid: %d / %d\n\n", sum(valid), n_reps))

cat("=== w1 (model 1 weight) ===\n")
print(quantile(r["w1", ], c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1)))
cat(sprintf("w1 > 0.9: %d\n", sum(r["w1", ] > 0.9)))
cat(sprintf("w1 > 0.99: %d\n", sum(r["w1", ] > 0.99)))

cat("\n=== w2 (model 2 weight) ===\n")
print(quantile(r["w2", ], c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1)))
cat(sprintf("w2 > 0.5: %d\n", sum(r["w2", ] > 0.5)))
cat(sprintf("w2 > 0.9: %d\n", sum(r["w2", ] > 0.9)))
cat(sprintf("w2 > 0.99: %d\n", sum(r["w2", ] > 0.99)))

cat("\n=== mu_ipw by w2 group ===\n")
grp1 <- r["w2", ] < 0.1
grp2 <- r["w2", ] >= 0.1 & r["w2", ] < 0.9
grp3 <- r["w2", ] >= 0.9
cat(sprintf("w2<0.1:  n=%d, mean(mu)=%.4f, sd(mu)=%.4f\n", sum(grp1), mean(r["mu_ipw", grp1]), sd(r["mu_ipw", grp1])))
if (sum(grp2) > 1) cat(sprintf("0.1<=w2<0.9: n=%d, mean(mu)=%.4f, sd(mu)=%.4f\n", sum(grp2), mean(r["mu_ipw", grp2]), sd(r["mu_ipw", grp2])))
if (sum(grp3) > 1) cat(sprintf("w2>=0.9: n=%d, mean(mu)=%.4f, sd(mu)=%.4f\n", sum(grp3), mean(r["mu_ipw", grp3]), sd(r["mu_ipw", grp3])))
