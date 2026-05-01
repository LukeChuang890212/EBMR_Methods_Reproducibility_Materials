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
mu_true <- get_mu_true("setting4")

for (rep_i in c(1, 2, 3, 5, 10, 15, 20, 50)) {
  set.seed(12345 + rep_i)
  dat <- setting4.B1(n_val)

  # Standard (unconstrained) run
  ebmr <- EBMRAlgorithmFast4[["new"]]("y", ps_spec_13, dat, W_func)
  res_std <- ebmr[["EBMR_IPW"]](
    h_nu = function(dat) cbind(u1 = dat[["u1"]], u2 = dat[["u2"]], z1 = dat[["z1"]], z2 = dat[["z2"]], u1_u2 = dat[["u1"]]*dat[["u2"]]),
    se.fit = FALSE, type = "HT"
  )

  # Constrained: reparameterize with p in [0,1]
  # ensemble_ps = p * ps1 + (1-p) * ps2
  ps.matrix <- do.call(cbind, lapply(ebmr[["ps_fit.list"]], function(x) x[["fitted.values"]]))
  r_vec <- as.vector(dat[["r"]])
  y_vec <- dat[["y"]]
  n <- n_val

  h_x2 <- cbind(u1 = dat[["u1"]], u2 = dat[["u2"]], z1 = dat[["z1"]], z2 = dat[["z2"]], u1_u2 = dat[["u1"]]*dat[["u2"]])
  h_x <- cbind(1, h_x2)
  h_dim <- ncol(h_x)

  # GMM objective: Q(p) = G(p)' W(p) G(p)
  # where g_i(p) = (r_i / (p*ps1_i + (1-p)*ps2_i) - 1) * h_x_i
  gmm_obj <- function(p) {
    ps_ens <- p * ps.matrix[, 1] + (1 - p) * ps.matrix[, 2]
    g_mat <- (r_vec / ps_ens - 1) * h_x
    G_vec <- colMeans(g_mat)
    W_mat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(h_dim))
    as.numeric(t(G_vec) %*% W_mat %*% G_vec)
  }

  # Optimize p in [0, 1]
  opt <- optimize(gmm_obj, interval = c(0, 1))
  p_hat <- opt$minimum

  # Compute mu_ipw with constrained weight
  ps_ens_constrained <- p_hat * ps.matrix[, 1] + (1 - p_hat) * ps.matrix[, 2]
  mu_constrained <- mean(r_vec / ps_ens_constrained * y_vec)

  cat(sprintf("Rep %2d: std w=(%5.3f,%5.3f) mu=%.4f | constrained p=%.4f mu=%.4f | diff=%.4f\n",
      rep_i,
      res_std[["w.hat"]][1], res_std[["w.hat"]][2], res_std[["mu_ipw"]],
      p_hat, mu_constrained,
      mu_constrained - res_std[["mu_ipw"]]))
}
