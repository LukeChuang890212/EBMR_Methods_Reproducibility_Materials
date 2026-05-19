## Verify Fast5 SE formula against bootstrap
## Check on multiple reps: analytical SE vs bootstrap SE
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast5", quiet = TRUE)

mu_true <- get_mu_true("setting4")
ps_spec_base <- get_ps_spec("9")
h_alpha_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2)
W_sm <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_nu_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1u2=dat$u1*dat$u2)

nn <- 2000; B <- 500  # bootstrap reps

# Test on mu_011 (the problematic one) and mu_101
for (est_label in c("mu_011", "mu_101", "mu_001")) {
  if (est_label == "mu_011") {
    model_idx <- c(2, 3)
  } else if (est_label == "mu_101") {
    model_idx <- c(1, 3)
  } else {
    model_idx <- 3
  }
  J <- length(model_idx)

  ps_spec_sub <- list(
    formula.list = ps_spec_base$formula.list[model_idx],
    h_alpha.list = lapply(seq_len(J), function(j) h_alpha_fn),
    inv_link = ps_spec_base$inv_link,
    outcome = ps_spec_base$outcome
  )

  cat(sprintf("\n=== %s: analytical vs bootstrap SE (B=%d) ===\n", est_label, B))
  cat(sprintf("%-4s | %-10s | %-10s | %-10s | %s\n",
      "Rep", "SE_analyt", "SE_boot", "ratio", "w1_med"))
  cat(strrep("-", 60), "\n")

  for (rep in 1:10) {
    set.seed(rep)
    dat <- setting4.B1(nn)

    # Analytical SE
    ebmr <- EBMRAlgorithmFast5$new("y", ps_spec_sub, dat, W_sm)
    if (J == 1) {
      ipw <- ebmr$EBMR_IPW(h_nu_fn, se.fit = TRUE)
    } else {
      ipw <- ebmr$EBMR_IPW(h_nu_fn, model_indices = 1:J, se.fit = TRUE)
    }
    se_analyt <- ipw$se_ipw
    mu_hat <- ipw$mu_ipw

    # Bootstrap SE
    boot_mu <- numeric(B)
    for (b in 1:B) {
      idx <- sample(nn, nn, replace = TRUE)
      dat_b <- dat[idx, ]
      ebmr_b <- tryCatch({
        eb <- EBMRAlgorithmFast5$new("y", ps_spec_sub, dat_b, W_sm)
        if (J == 1) {
          ipw_b <- eb$EBMR_IPW(h_nu_fn, se.fit = FALSE)
        } else {
          ipw_b <- eb$EBMR_IPW(h_nu_fn, model_indices = 1:J, se.fit = FALSE)
        }
        ipw_b$mu_ipw
      }, error = function(e) NA)
      boot_mu[b] <- ebmr_b
    }
    se_boot <- sd(boot_mu, na.rm = TRUE)
    na_boot <- sum(is.na(boot_mu))

    w1_str <- if (J > 1) sprintf("%.3f", median(ipw$w.hat)) else "—"
    cat(sprintf("%4d | %.6f  | %.6f  | %.3f     | %s  (boot NA=%d)\n",
        rep, se_analyt, se_boot, se_analyt / se_boot, w1_str, na_boot))
  }
}
