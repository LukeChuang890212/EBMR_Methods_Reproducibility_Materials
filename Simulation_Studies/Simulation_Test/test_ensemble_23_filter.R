setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

n_val <- 2000; n_reps <- 200
ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
mu_true <- get_mu_true("setting4")
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_nu_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1_u2=dat$u1*dat$u2)
model_set <- c(2, 3)

mu_vals <- se_vals <- rep(NA_real_, n_reps)
cond_m2 <- cond_m3 <- rep(NA_real_, n_reps)

for (rep_i in 1:n_reps) {
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
  ps_sub <- list(
    formula.list = ps_spec[["formula.list"]][model_set],
    h_alpha.list = ps_spec[["h_alpha.list"]][model_set],
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL, NULL),
    optimizer = "constrained_nr"
  )
  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", ps_sub, dat, W_fn)
    res <- ebmr$EBMR_IPW(h_nu = h_nu_fn, type = "HT", se.fit = TRUE)
    mu_vals[rep_i] <- res$mu_ipw
    se_vals[rep_i] <- res$se_ipw
    cond_m2[rep_i] <- ebmr$ps_fit.list[[1]]$gmm_fit$opt$solution_cond
    cond_m3[rep_i] <- ebmr$ps_fit.list[[2]]$gmm_fit$opt$solution_cond
  }, error = function(e) NULL)
}

v <- !is.na(mu_vals) & !is.na(se_vals)

cat(sprintf("All valid: n=%d, Bias=%.4f, ESD=%.4f, ESE=%.4f, ESE/ESD=%.3f\n",
    sum(v), mean(mu_vals[v]) - mu_true, sd(mu_vals[v]),
    mean(se_vals[v]), mean(se_vals[v]) / sd(mu_vals[v])))

no_m3d <- v & cond_m3 < 1e8
cat(sprintf("Remove M3 degen: n=%d, Bias=%.4f, ESD=%.4f, ESE=%.4f, ESE/ESD=%.3f\n",
    sum(no_m3d), mean(mu_vals[no_m3d]) - mu_true, sd(mu_vals[no_m3d]),
    mean(se_vals[no_m3d]), mean(se_vals[no_m3d]) / sd(mu_vals[no_m3d])))

no_m2d <- v & cond_m2 < 1e8
cat(sprintf("Remove M2 degen: n=%d\n", sum(no_m2d)))

no_either <- v & cond_m2 < 1e8 & cond_m3 < 1e8
cat(sprintf("Remove either degen: n=%d\n", sum(no_either)))

cat(sprintf("\nM2 degen count: %d/%d\n", sum(v & cond_m2 >= 1e8), sum(v)))
cat(sprintf("M3 degen count: %d/%d\n", sum(v & cond_m3 >= 1e8), sum(v)))
