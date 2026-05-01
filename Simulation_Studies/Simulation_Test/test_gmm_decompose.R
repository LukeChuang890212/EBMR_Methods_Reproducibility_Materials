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

for (rep_i in c(3, 1, 2, 10)) {
  set.seed(12345 + rep_i)
  dat <- setting4.B1(n_val)

  ebmr <- EBMRAlgorithmFast4[["new"]]("y", ps_spec_13, dat, W_func)
  ps_fit2 <- ebmr[["ps_fit.list"]][[2]]
  alpha2 <- ps_fit2[["coefficients"]]
  pi2 <- ps_fit2[["fitted.values"]]
  r_vec <- as.vector(dat$r)
  n <- n_val

  # Get h_x from ps_fit
  h_x <- ps_fit2[["h_x"]]
  h_dim <- ncol(h_x)

  # Compute g_mat = (r/pi - 1) * h_x
  rpi <- r_vec / pi2
  g_mat <- (rpi - 1) * h_x

  # Moment conditions G = colMeans(g_mat)
  G_vec <- colMeans(g_mat)

  # W = inv(g'g/n)
  W_mat <- solve(crossprod(g_mat) / n)

  # GMM objective = G' W G
  obj <- as.numeric(t(G_vec) %*% W_mat %*% G_vec)

  cat(sprintf("\n============================================================\n"))
  cat(sprintf("=== Rep %d (||alpha||=%.2f) ===\n", rep_i, sqrt(sum(alpha2^2))))
  cat(sprintf("============================================================\n"))

  cat(sprintf("\n--- r/pi - 1 distribution ---\n"))
  rpi_minus1 <- rpi - 1
  cat(sprintf("  Overall: mean=%.4f, sd=%.4f, range=[%.2f, %.2f]\n",
      mean(rpi_minus1), sd(rpi_minus1), min(rpi_minus1), max(rpi_minus1)))
  cat(sprintf("  r=0: mean=%.4f, sd=%.4f (all = -1 by definition)\n",
      mean(rpi_minus1[r_vec == 0]), sd(rpi_minus1[r_vec == 0])))
  cat(sprintf("  r=1: mean=%.4f, sd=%.4f, range=[%.2f, %.2f]\n",
      mean(rpi_minus1[r_vec == 1]), sd(rpi_minus1[r_vec == 1]),
      min(rpi_minus1[r_vec == 1]), max(rpi_minus1[r_vec == 1])))
  cat(sprintf("  r=1 & pi<0.06: n=%d, r/pi-1 range=[%.1f, %.1f]\n",
      sum(r_vec == 1 & pi2 < 0.06),
      ifelse(any(r_vec == 1 & pi2 < 0.06), min(rpi_minus1[r_vec == 1 & pi2 < 0.06]), NA),
      ifelse(any(r_vec == 1 & pi2 < 0.06), max(rpi_minus1[r_vec == 1 & pi2 < 0.06]), NA)))

  cat(sprintf("\n--- Moment conditions G = E[(r/pi-1)*h_x] ---\n"))
  cat(sprintf("  G: %s\n", paste(round(G_vec, 6), collapse = ", ")))
  cat(sprintf("  ||G||: %.6f\n", sqrt(sum(G_vec^2))))

  cat(sprintf("\n--- Variance of g: diag(g'g/n) ---\n"))
  var_g <- diag(crossprod(g_mat) / n)
  cat(sprintf("  diag(g'g/n): %s\n", paste(round(var_g, 2), collapse = ", ")))
  cat(sprintf("  max var: %.2f\n", max(var_g)))

  cat(sprintf("\n--- W = inv(g'g/n) ---\n"))
  cat(sprintf("  diag(W): %s\n", paste(round(diag(W_mat), 6), collapse = ", ")))
  cat(sprintf("  max |W|: %.6f\n", max(abs(W_mat))))

  cat(sprintf("\n--- GMM objective decomposition ---\n"))
  cat(sprintf("  G'WG = %.6f\n", obj))
  cat(sprintf("  ||G||^2 = %.6f\n", sum(G_vec^2)))
  cat(sprintf("  ||G||^2 * max_eig(W) <= %.6f\n", sum(G_vec^2) * max(eigen(W_mat)$values)))
  cat(sprintf("  eigenvalues of W: %s\n", paste(round(eigen(W_mat)$values, 6), collapse = ", ")))

  # Compare: what if we used W = I (unweighted)?
  obj_I <- sum(G_vec^2)
  cat(sprintf("\n  G'G (W=I): %.6f\n", obj_I))
  cat(sprintf("  G'WG / G'G ratio: %.6f\n", obj / obj_I))
}
