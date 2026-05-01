setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
source("Simulation.r")
library(EBMRalgorithmFast4)

ps_spec <- get_ps_spec("9-alt1")
ps_spec_13 <- list(
  formula.list = ps_spec[["formula.list"]][c(1, 3)],
  h_alpha.list = ps_spec[["h_alpha.list"]][c(1, 3)],
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]]
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
n_val <- 2000

# Compare a "good" rep (1) and "bad" rep (2)
for (rep_i in c(1, 2, 10)) {
  set.seed(12345 + rep_i)
  dat <- setting4.B1(n_val)

  ebmr <- EBMRAlgorithmFast4[["new"]]("y", ps_spec_13, dat, W_func)

  cat(sprintf("\n============================================================\n"))
  cat(sprintf("=== Rep %d ===\n", rep_i))
  cat(sprintf("============================================================\n"))

  r_vec <- as.vector(dat[["r"]])
  y_vec <- dat[["y"]]
  n <- n_val
  cat(sprintf("Response rate: %.3f\n", mean(r_vec)))

  for (j in 1:2) {
    ps_fit <- ebmr[["ps_fit.list"]][[j]]
    alpha <- ps_fit[["coefficients"]]
    pi_hat <- ps_fit[["fitted.values"]]
    gmm_fit <- ps_fit[["gmm_fit"]]

    cat(sprintf("\n--- Model %d ---\n", j))
    cat(sprintf("  Formula: %s\n", deparse(ps_spec_13[["formula.list"]][[j]])))
    cat(sprintf("  h_alpha: (function, %d cols)\n", ncol(ps_spec_13[["h_alpha.list"]][[j]](dat))))
    cat(sprintf("  alpha: %s\n", paste(round(alpha, 4), collapse=", ")))
    cat(sprintf("  ||alpha||: %.4f\n", sqrt(sum(alpha^2))))
    cat(sprintf("  GMM obj: %.6f, converged: %s, iter: %d\n",
        gmm_fit[["opt"]][["objective"]], gmm_fit[["opt"]][["converged"]], gmm_fit[["opt"]][["iterations"]]))

    # Pi distribution
    cat(sprintf("  pi quantiles: %.4f %.4f %.4f %.4f %.4f\n",
        quantile(pi_hat, 0.01), quantile(pi_hat, 0.10),
        quantile(pi_hat, 0.50), quantile(pi_hat, 0.90), quantile(pi_hat, 0.99)))
    cat(sprintf("  pi at boundaries (<0.06 or >0.94): %d / %d\n",
        sum(pi_hat < 0.06 | pi_hat > 0.94), n))

    # How well does pi predict r?
    cat(sprintf("  Correlation(pi, r): %.4f\n", cor(pi_hat, r_vec)))

    # Check moment conditions: E[(r/pi - 1) * h_x]
    h_alpha_func <- ps_spec_13[["h_alpha.list"]][[j]]
    h_x_raw <- h_alpha_func(dat)
    h_x <- cbind(1, h_x_raw)
    g_mat <- (r_vec / pi_hat - 1) * h_x
    G_vec <- colMeans(g_mat)
    cat(sprintf("  ||G|| (moment conditions): %.6f\n", sqrt(sum(G_vec^2))))
    cat(sprintf("  G: %s\n", paste(round(G_vec, 6), collapse=", ")))

    # IPW weights: r/pi
    ipw <- r_vec / pi_hat
    ipw_r1 <- ipw[r_vec == 1]
    cat(sprintf("  IPW weights (r=1): mean=%.2f, max=%.2f, sd=%.2f\n",
        mean(ipw_r1), max(ipw_r1), sd(ipw_r1)))
    cat(sprintf("  sum(r/pi)/n = %.4f (should be ~1)\n", mean(ipw)))

    # mu_ipw for this model alone
    mu_j <- mean(r_vec / pi_hat * y_vec)
    cat(sprintf("  mu_ipw (this model): %.4f\n", mu_j))
  }

  # Ensemble result
  res <- ebmr[["EBMR_IPW"]](
    h_nu = function(dat) cbind(u1 = dat[["u1"]], u2 = dat[["u2"]], z1 = dat[["z1"]], z2 = dat[["z2"]], u1_u2 = dat[["u1"]]*dat[["u2"]]),
    se.fit = FALSE, type = "HT"
  )
  cat(sprintf("\n--- Ensemble ---\n"))
  cat(sprintf("  nu: %s\n", paste(round(res[["nu.hat"]], 4), collapse=", ")))
  cat(sprintf("  w:  %s\n", paste(round(res[["w.hat"]], 4), collapse=", ")))
  cat(sprintf("  mu_ipw: %.4f (true: 0.7536)\n", res[["mu_ipw"]]))

  # Check: what is mu_ipw at p=1 (model 1 only)?
  ps1 <- ebmr[["ps_fit.list"]][[1]][["fitted.values"]]
  ps2 <- ebmr[["ps_fit.list"]][[2]][["fitted.values"]]
  cat(sprintf("  mu_ipw(model1 only): %.4f\n", mean(r_vec / ps1 * y_vec)))
  cat(sprintf("  mu_ipw(model2 only): %.4f\n", mean(r_vec / ps2 * y_vec)))
}
