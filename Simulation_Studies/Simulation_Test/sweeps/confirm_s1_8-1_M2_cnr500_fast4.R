## Confirm cnr/500 for s1 8-1 M2 (both miss50 and miss30), Fast4, 1000 reps each
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({ source("Basic_setup.r"); source("Data_Generation.r"); source("Simulation.r") })
devtools::load_all("../EPS", quiet = TRUE)
library(parallel); library(foreach); library(doSNOW)

W_sm <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_full <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2)
h_nu_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1_u2=dat$u1*dat$u2)
inv_link_fn <- function(eta) 1/(1+exp(eta))

ps_spec <- list(
  formula.list = list(r ~ y + u1 + z1),
  h_alpha.list = list(h_full),
  inv_link = inv_link_fn, outcome = "y",
  optimizer = "constrained_nr", cond_threshold = 500
)

ps_model.true <- function(dat, alpha.true) {
  X <- cbind(rep(1, nrow(dat)), dat$y, dat$u1, dat$u2)
  1 / (1 + exp(X %*% alpha.true))
}

n_cores <- min(detectCores() - 1, 10); n_reps <- 1000; nn <- 2000

configs_data <- list(
  list(miss = "miss50", dgp = "s1.A1",
       data_file = "Simulation_Data/setting1.A1_n2000_replicate1000.RDS",
       alpha.true = correct_model_alpha.true.list$setting1$miss50[[1]]),
  list(miss = "miss30", dgp = "s1.A2",
       data_file = "Simulation_Data/setting1.A2_n2000_replicate1000.RDS",
       alpha.true = correct_model_alpha.true.list$setting1$miss30[[1]])
)

for (cfg in configs_data) {
  all_data <- readRDS(cfg$data_file)
  alpha.true <- cfg$alpha.true
  mu_true <- mean(all_data$y)
  cat(sprintf("\n=== %s %s 8-1 M2, cnr/500 Fast4, %d reps (mu_true=%.4f) ===\n",
              cfg$dgp, cfg$miss, n_reps, mu_true))
  t0 <- proc.time()
  cl <- makeCluster(n_cores); registerDoSNOW(cl)
  clusterExport(cl, c("nn","ps_spec","W_sm","h_full","h_nu_fn","inv_link_fn",
                      "all_data","alpha.true","ps_model.true"), envir = environment())
  clusterEvalQ(cl, {
    setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
    devtools::load_all("../EPS", quiet = TRUE)
  })
  pb <- txtProgressBar(max = n_reps, style = 3)
  opts <- list(progress = function(n) setTxtProgressBar(pb, n))
  res <- foreach(i = 1:n_reps, .combine = 'cbind', .options.snow = opts,
                 .packages = c("stringr","Matrix")) %dopar% {
    tryCatch({
      dat <- all_data[((i - 1) * nn + 1):(i * nn), ]
      ebmr <- EPS$new("y", ps_spec, dat, W_sm)
      ipw <- ebmr$EBMR_IPW(h_nu_fn, model_indices = 1, se.fit = TRUE,
                            true_ps = ps_model.true(dat, alpha.true))
      a <- unname(ebmr$ps_fit.list[[1]]$coefficients)
      c(mu_ipw = ipw$mu_ipw, mu_ipw.true = ipw$mu_ipw.true,
        se_ipw = ipw$se_ipw, se_ipw.true = ipw$se_ipw.true,
        alpha.hat1 = a[1], alpha.hat2 = a[2], alpha.hat3 = a[3], alpha.hat4 = a[4])
    }, error = function(e) rep(NA_real_, 8))
  }
  close(pb); stopCluster(cl)
  elapsed <- (proc.time() - t0)["elapsed"]
  cleaned <- clean_sim_result(res, multiplier = 2, max_pct = 0.01, verbose = FALSE)
  sc <- cleaned$result
  mu_v <- sc["mu_ipw", ]; se_v <- sc["se_ipw", ]
  bias <- mean(mu_v) - mu_true; esd <- sd(mu_v); ese <- mean(se_v); ese_med <- median(se_v)
  ratio <- ese / esd
  ci_lo <- mu_v - 1.96*se_v; ci_hi <- mu_v + 1.96*se_v
  cp <- mean(ci_lo <= mu_true & mu_true <= ci_hi)
  alpha_mat <- res[c("alpha.hat1","alpha.hat2","alpha.hat3","alpha.hat4"), ]
  n_zero <- sum(colSums(abs(alpha_mat), na.rm = TRUE) == 0)
  max_abs <- apply(alpha_mat, 2, function(x) max(abs(x), na.rm = TRUE))
  n_sep <- sum(max_abs > 50, na.rm = TRUE)
  sd_a <- apply(alpha_mat, 1, sd, na.rm = TRUE)
  cat(sprintf("stuck@0=%d sep(>50)=%d NA=%d Out=%d  Bias=%+.4f ESD=%.4f ESE(mean)=%.4f ESE(median)=%.4f ratio=%.3f CP=%.3f\n",
              n_zero, n_sep, cleaned$n_na, cleaned$n_outliers,
              bias, esd, ese, ese_med, ratio, cp))
  cat(sprintf("sd(alpha)=(%.3f,%.3f,%.3f,%.3f)  max|alpha|=%.1f  elapsed=%.0fs\n",
              sd_a[1], sd_a[2], sd_a[3], sd_a[4], max(max_abs, na.rm = TRUE), elapsed))
}
cat("\nDONE\n")
