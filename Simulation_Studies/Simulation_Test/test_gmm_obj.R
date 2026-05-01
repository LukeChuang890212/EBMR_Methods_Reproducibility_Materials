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

for (rep_i in c(1, 2, 3, 10)) {
  set.seed(12345 + rep_i)
  dat <- setting4.B1(n_val)

  ebmr <- EBMRAlgorithmFast4[["new"]]("y", ps_spec_13, dat, W_func)

  ps.matrix <- do.call(cbind, lapply(ebmr[["ps_fit.list"]], function(x) x[["fitted.values"]]))
  r_vec <- as.vector(dat[["r"]])
  n <- n_val

  h_x2 <- cbind(u1 = dat[["u1"]], u2 = dat[["u2"]], z1 = dat[["z1"]], z2 = dat[["z2"]], u1_u2 = dat[["u1"]]*dat[["u2"]])
  h_x <- cbind(1, h_x2)
  h_dim <- ncol(h_x)

  gmm_obj <- function(p) {
    ps_ens <- p * ps.matrix[, 1] + (1 - p) * ps.matrix[, 2]
    g_mat <- (r_vec / ps_ens - 1) * h_x
    G_vec <- colMeans(g_mat)
    W_mat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(h_dim))
    as.numeric(t(G_vec) %*% W_mat %*% G_vec)
  }

  cat(sprintf("\n=== Rep %d ===\n", rep_i))
  cat(sprintf("  PS1 alpha obj: %.6f\n", ebmr[["ps_fit.list"]][[1]][["gmm_fit"]][["opt"]][["objective"]]))
  cat(sprintf("  PS2 alpha obj: %.6f\n", ebmr[["ps_fit.list"]][[2]][["gmm_fit"]][["opt"]][["objective"]]))

  # Evaluate ensemble GMM objective at various p values
  p_vals <- c(0, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1)
  for (p in p_vals) {
    cat(sprintf("  p=%.2f (w1=%.2f, w2=%.2f): obj=%.6f\n", p, p, 1-p, gmm_obj(p)))
  }

  opt <- optimize(gmm_obj, interval = c(0, 1))
  cat(sprintf("  optimal p=%.6f, obj=%.6f\n", opt$minimum, opt$objective))

  # Also check the unconstrained nu from the package
  res <- ebmr[["EBMR_IPW"]](
    h_nu = function(dat) cbind(u1 = dat[["u1"]], u2 = dat[["u2"]], z1 = dat[["z1"]], z2 = dat[["z2"]], u1_u2 = dat[["u1"]]*dat[["u2"]]),
    se.fit = FALSE, type = "HT"
  )
  cat(sprintf("  package nu=(%s), w=(%s)\n",
      paste(round(res[["nu.hat"]], 4), collapse=", "),
      paste(round(res[["w.hat"]], 4), collapse=", ")))
}
