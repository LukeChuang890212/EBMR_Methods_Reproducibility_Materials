setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Data_Generation.r")
source("config/scenarios.R")
library(EBMRalgorithmFast4)
source("Basic_setup.r")

data_file <- misspecified_model_all_data_file.list$setting3$miss50[[1]]
all_data <- readRDS(data_file)
n_val <- 2000

ps_spec <- get_ps_spec("9-alt1")
subset_ps_spec <- list(
  formula.list = ps_spec$formula.list[2],
  h_alpha.list = ps_spec$h_alpha.list[2],
  inv_link = ps_spec$inv_link,
  outcome = ps_spec$outcome
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

for (rep_i in c(437, 15, 1, 4)) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]

  ebmr <- EBMRAlgorithmFast4$new("y", subset_ps_spec, dat, W_func)

  ps_fit <- ebmr$ps_fit.list[[1]]
  ps_fitted <- ps_fit$fitted.values
  r <- as.numeric(dat$r)
  y <- dat$y
  n <- length(y)

  gmm_fit <- ps_fit$gmm_fit
  psi_alpha <- gmm_fit$psi
  Gamma_hat <- gmm_fit$Gamma.hat
  W_hat <- gmm_fit$W.hat

  alpha_hat <- ps_fit$coefficients
  alpha_dim <- length(alpha_hat)
  h_dim <- nrow(psi_alpha)

  # h_alpha is design_matrix (the covariates used in PS model)
  h_alpha_mat <- as.matrix(ps_fit$design_matrix)

  cat(sprintf("========== Rep %d ==========\n", rep_i))
  cat(sprintf("alpha = [%s]\n", paste(round(alpha_hat, 4), collapse=", ")))
  cat(sprintf("alpha_dim=%d, h_dim=%d, dim(design_matrix)=%dx%d\n",
              alpha_dim, h_dim, nrow(h_alpha_mat), ncol(h_alpha_mat)))

  # dot_pi = d(ps)/d(alpha) = h_alpha * ps * (1-ps) for logistic
  dot_pi <- h_alpha_mat * ps_fitted * (1 - ps_fitted)  # n x alpha_dim

  # ry_ps_inv2 = r * y / ps^2
  ry_ps_inv2 <- as.vector(r * y * ps_fitted^(-2))

  # H_alpha_w = (1/n) sum_i dot_pi_i * ry_ps_inv2_i
  H_alpha_w <- colMeans(dot_pi * ry_ps_inv2)

  cat(sprintf("H_alpha_w = [%s]\n", paste(round(H_alpha_w, 4), collapse=", ")))
  cat(sprintf("||H_alpha_w|| = %.4f\n", sqrt(sum(H_alpha_w^2))))

  # alpha_adj_i = H_alpha_w' %*% psi_alpha_i
  alpha_adj <- as.vector(t(H_alpha_w) %*% psi_alpha)

  base_term <- r / ps_fitted * y
  mu_iid <- base_term - alpha_adj

  cat(sprintf("\nVariance decomposition (all /n):\n"))
  cat(sprintf("  var(base)/n     = %.6f  (se_naive = %.4f)\n", var(base_term)/n, sqrt(var(base_term)/n)))
  cat(sprintf("  var(adj)/n      = %.6f  (se_adj   = %.4f)\n", var(alpha_adj)/n, sqrt(var(alpha_adj)/n)))
  cat(sprintf("  -2cov(b,a)/n    = %.6f\n", -2*cov(base_term, alpha_adj)/n))
  cat(sprintf("  var(iid)/n      = %.6f  (se_total = %.4f)\n", var(mu_iid)/n, sqrt(var(mu_iid)/n)))
  cat(sprintf("  cor(base, adj)  = %.4f\n", cor(base_term, alpha_adj)))

  # Check psi_alpha structure
  cat(sprintf("\npsi_alpha row stats:\n"))
  for (j in 1:h_dim) {
    cat(sprintf("  psi[%d,]: sd=%.4f, max|.|=%.4f\n", j, sd(psi_alpha[j,]), max(abs(psi_alpha[j,]))))
  }

  # Top 5 |alpha_adj|
  top5_adj <- order(abs(alpha_adj), decreasing=TRUE)[1:5]
  cat(sprintf("\nTop 5 |alpha_adj|:\n"))
  for (k in top5_adj) {
    cat(sprintf("  obs %d: adj=%.2f, base=%.2f, ps=%.5f, r=%d, y=%.3f\n",
                k, alpha_adj[k], base_term[k], ps_fitted[k], r[k], y[k]))
  }

  # GMM sandwich condition
  GtW <- crossprod(Gamma_hat, W_hat)
  GtWG <- GtW %*% Gamma_hat
  GtWG_inv <- solve(GtWG)

  cat(sprintf("\ncond(GtWG) = %.2e\n", kappa(GtWG)))
  cat(sprintf("||GtWG_inv|| (Frob) = %.4f\n", sqrt(sum(GtWG_inv^2))))

  # Effective weight: H' (GtWG)^{-1} GtW — this is the linear combo applied to h_i
  eff_weight <- as.vector(t(H_alpha_w) %*% GtWG_inv %*% GtW)
  cat(sprintf("Effective weight = [%s]\n", paste(round(eff_weight, 4), collapse=", ")))
  cat(sprintf("||eff_weight|| = %.4f\n", sqrt(sum(eff_weight^2))))

  # What is the moment condition h_i for the top obs?
  cat(sprintf("\nMoment conditions h_i for top |adj| obs:\n"))
  # h_i from the GMM is: (r_i/ps_i - 1) * h_alpha_i
  for (k in top5_adj[1:3]) {
    h_i <- (r[k]/ps_fitted[k] - 1) * h_alpha_mat[k,]
    cat(sprintf("  obs %d: h_i=[%s], r/ps=%.2f\n",
                k, paste(round(h_i, 2), collapse=", "), r[k]/ps_fitted[k]))
  }

  cat("\n\n")
}

cat("Done!\n")
