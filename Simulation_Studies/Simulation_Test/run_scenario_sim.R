source("Data_Generation.r")
source("Basic_setup.r")
source("config/scenarios.R")

old_wd <- getwd()
setwd("../EBMRalgorithmFast2")
source("R/EBMRAlgorithm.r")
setwd(old_wd)

set.seed(12345)
big13 <- setting13.A1(1e7); mu_true_13 <- mean(big13$y); rm(big13)
big15 <- setting15.A1(1e7); mu_true_15 <- mean(big15$y); rm(big15)
gc()
cat(sprintf("mu_true: setting13=%.6f, setting15=%.6f\n\n", mu_true_13, mu_true_15))

inv_link_fn <- function(eta) 1 / (1 + exp(eta))

# 3 scenarios: differ only in model 1 h_alpha
scenarios <- list(
  "8-2" = list(
    h1 = c("u1","u2","z1","z2"),
    h2 = c("u1","u2","z1","z2"),
    h3 = c("u1","u2","z1","z2")
  ),
  "8-3" = list(
    h1 = function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u2_sq=dat$u2^2),
    h2 = c("u1","u2","z1","z2"),
    h3 = c("u1","u2","z1","z2")
  ),
  "8-4" = list(
    h1 = function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, z2_sq=dat$z2^2),
    h2 = c("u1","u2","z1","z2"),
    h3 = c("u1","u2","z1","z2")
  )
)

# Formulas: same across all 3 scenarios
# Model 1: r ~ y + u1 + u2 (full)
# Model 2: r ~ y + u1 + z1
# Model 3: r ~ y + u2 + z2

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
n_sim <- 1000
n_val <- 2000

settings <- list(
  "setting13.B1" = list(func="setting13.B1", mu=NULL, label="Setting13 B1 (50% miss)"),
  "setting13.B2" = list(func="setting13.B2", mu=NULL, label="Setting13 B2 (30% miss)"),
  "setting15.B1" = list(func="setting15.B1", mu=NULL, label="Setting15 B1 (50% miss)"),
  "setting15.B2" = list(func="setting15.B2", mu=NULL, label="Setting15 B2 (30% miss)")
)
settings[["setting13.B1"]]$mu <- mu_true_13
settings[["setting13.B2"]]$mu <- mu_true_13
settings[["setting15.B1"]]$mu <- mu_true_15
settings[["setting15.B2"]]$mu <- mu_true_15

for (sname in names(settings)) {
  sinfo <- settings[[sname]]
  mu_true <- sinfo$mu

  cat(sprintf("\n%s\n%s\n", strrep("=", 90), sinfo$label))
  cat(sprintf("%s\n\n", strrep("=", 90)))

  for (scen_name in names(scenarios)) {
    scen <- scenarios[[scen_name]]

    ps_spec <- list(
      formula.list = list(r ~ y + u1 + u2, r ~ y + u1 + z1, r ~ y + u2 + z2),
      h_alpha.list = list(scen$h1, scen$h2, scen$h3),
      inv_link = inv_link_fn,
      outcome = "y",
      alpha_init.list = list(NULL, NULL, NULL)
    )

    mu_m1 <- mu_m2 <- mu_m3 <- numeric(n_sim)
    # Store max absolute moment condition for each model
    max_mc_m1 <- max_mc_m2 <- max_mc_m3 <- numeric(n_sim)
    # Store all moment condition means
    mc_m1_all <- mc_m2_all <- mc_m3_all <- list()
    n_err <- 0

    set.seed(42)
    for (i in 1:n_sim) {
      dat <- get(sinfo$func)(n_val)
      tryCatch({
        ebmr <- EBMRAlgorithmFast2$new("y", ps_spec, dat, W_func)

        for (j in 1:3) {
          ps_j <- ebmr$ps_fit.list[[j]]$fitted.values
          h_x_j <- ebmr$ps_fit.list[[j]]$h_x
          mu_j <- mean(dat$r / ps_j * dat$y)
          rw <- dat$r / ps_j
          mc <- colMeans((rw - 1) * h_x_j)

          if (j == 1) { mu_m1[i] <- mu_j; max_mc_m1[i] <- max(abs(mc)); mc_m1_all[[i]] <- mc }
          if (j == 2) { mu_m2[i] <- mu_j; max_mc_m2[i] <- max(abs(mc)); mc_m2_all[[i]] <- mc }
          if (j == 3) { mu_m3[i] <- mu_j; max_mc_m3[i] <- max(abs(mc)); mc_m3_all[[i]] <- mc }
        }
      }, error = function(e) {
        mu_m1[i] <<- NA; mu_m2[i] <<- NA; mu_m3[i] <<- NA
        max_mc_m1[i] <<- NA; max_mc_m2[i] <<- NA; max_mc_m3[i] <<- NA
        n_err <<- n_err + 1
      })
      if (i %% 500 == 0) cat(sprintf("  Scenario %s: %d/%d (errors: %d)\n", scen_name, i, n_sim, n_err))
    }

    valid <- !is.na(mu_m1)
    cat(sprintf("\n--- Scenario %s (%d/%d valid) ---\n", scen_name, sum(valid), n_sim))

    # IPW bias
    cat(sprintf("  Model 1 (y+u1+u2):  bias = %8.4f, RMSE = %.4f\n",
                mean(mu_m1[valid]) - mu_true, sqrt(mean((mu_m1[valid] - mu_true)^2))))
    cat(sprintf("  Model 2 (y+u1+z1):  bias = %8.4f, RMSE = %.4f\n",
                mean(mu_m2[valid]) - mu_true, sqrt(mean((mu_m2[valid] - mu_true)^2))))
    cat(sprintf("  Model 3 (y+u2+z2):  bias = %8.4f, RMSE = %.4f\n",
                mean(mu_m3[valid]) - mu_true, sqrt(mean((mu_m3[valid] - mu_true)^2))))

    # Mean max|moment condition|
    cat(sprintf("\n  Mean max|MC|:  M1=%.6f  M2=%.6f  M3=%.6f\n",
                mean(max_mc_m1[valid]), mean(max_mc_m2[valid]), mean(max_mc_m3[valid])))

    # Average moment conditions across replicates
    mc1_avg <- rowMeans(do.call(cbind, mc_m1_all[valid]))
    mc2_avg <- rowMeans(do.call(cbind, mc_m2_all[valid]))
    mc3_avg <- rowMeans(do.call(cbind, mc_m3_all[valid]))

    h1_names <- names(mc_m1_all[[which(valid)[1]]])
    h2_names <- names(mc_m2_all[[which(valid)[1]]])
    h3_names <- names(mc_m3_all[[which(valid)[1]]])

    cat("\n  Avg moment conditions across 1000 reps:\n")
    cat("  Model 1 (y+u1+u2):", sprintf(" %s=%.6f", h1_names, mc1_avg), "\n")
    cat("  Model 2 (y+u1+z1):", sprintf(" %s=%.6f", h2_names, mc2_avg), "\n")
    cat("  Model 3 (y+u2+z2):", sprintf(" %s=%.6f", h3_names, mc3_avg), "\n")
    cat("\n")
  }
}
