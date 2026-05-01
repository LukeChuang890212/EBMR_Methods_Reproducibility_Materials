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

  cat(sprintf("\n=== Rep %d ===\n", rep_i))
  cat(sprintf("Response rate: %.3f\n", mean(dat$r)))

  # Model 2 formula: r ~ y + u2 + z2
  # Check if logistic regression has perfect separation
  glm_fit <- glm(r ~ y + u2 + z2, data = dat, family = binomial())
  cat(sprintf("GLM converged: %s\n", glm_fit$converged))
  cat(sprintf("GLM coefficients: %s\n", paste(round(coef(glm_fit), 4), collapse = ", ")))
  cat(sprintf("GLM max |coef|: %.4f\n", max(abs(coef(glm_fit)))))
  cat(sprintf("GLM deviance: %.2f, null deviance: %.2f\n", glm_fit$deviance, glm_fit$null.deviance))

  # Check predicted probabilities from GLM
  p_glm <- predict(glm_fit, type = "response")
  cat(sprintf("GLM p range: [%.6f, %.6f]\n", min(p_glm), max(p_glm)))
  cat(sprintf("GLM p < 0.01: %d, p > 0.99: %d\n", sum(p_glm < 0.01), sum(p_glm > 0.99)))
  cat(sprintf("GLM p < 0.05: %d, p > 0.95: %d\n", sum(p_glm < 0.05), sum(p_glm > 0.95)))

  # Check linear predictor range
  eta_glm <- predict(glm_fit, type = "link")
  cat(sprintf("GLM eta range: [%.2f, %.2f]\n", min(eta_glm), max(eta_glm)))

  # Now fit via the package (GMM with MNAR)
  ebmr <- EBMRAlgorithmFast4[["new"]]("y", ps_spec_13, dat, W_func)
  ps_fit2 <- ebmr[["ps_fit.list"]][[2]]
  alpha2 <- ps_fit2[["coefficients"]]
  pi2 <- ps_fit2[["fitted.values"]]

  cat(sprintf("\nGMM Model 2:\n"))
  cat(sprintf("  alpha: %s\n", paste(round(alpha2, 4), collapse = ", ")))
  cat(sprintf("  ||alpha||: %.4f\n", sqrt(sum(alpha2^2))))
  cat(sprintf("  pi range: [%.6f, %.6f]\n", min(pi2), max(pi2)))
  cat(sprintf("  pi < 0.06: %d, pi > 0.94: %d\n", sum(pi2 < 0.06), sum(pi2 > 0.94)))
  cat(sprintf("  GMM obj: %.6f, converged: %s\n",
      ps_fit2[["gmm_fit"]][["opt"]][["objective"]],
      ps_fit2[["gmm_fit"]][["opt"]][["converged"]]))

  # Check: does the design matrix for model 2 separate r perfectly?
  # Model 2 uses r ~ y + u2 + z2 (MNAR: y is in the model)
  # But y is only observed when r=1, so the design matrix uses y*r
  # Check the actual design matrix
  dm2 <- ps_fit2[["design_matrix"]]
  cat(sprintf("\n  Design matrix columns: %d\n", ncol(dm2)))

  # For r=0 observations, what does the linear predictor look like?
  r_vec <- dat$r
  eta2 <- as.vector(dm2 %*% alpha2)
  cat(sprintf("  eta (r=0): mean=%.2f, sd=%.2f, range=[%.2f, %.2f]\n",
      mean(eta2[r_vec == 0]), sd(eta2[r_vec == 0]),
      min(eta2[r_vec == 0]), max(eta2[r_vec == 0])))
  cat(sprintf("  eta (r=1): mean=%.2f, sd=%.2f, range=[%.2f, %.2f]\n",
      mean(eta2[r_vec == 1]), sd(eta2[r_vec == 1]),
      min(eta2[r_vec == 1]), max(eta2[r_vec == 1])))

  # Check overlap of eta between r=0 and r=1
  cat(sprintf("  eta overlap: r=0 max=%.2f vs r=1 min=%.2f (gap=%.2f)\n",
      max(eta2[r_vec == 0]), min(eta2[r_vec == 1]),
      min(eta2[r_vec == 1]) - max(eta2[r_vec == 0])))

  # Check if y*r creates a separation issue
  # When r=0, y*r = 0 for all obs; when r=1, y*r = y which varies
  # If alpha_y (coef on y) is large, this pushes r=1 obs to extreme pi
  cat(sprintf("\n  y distribution (r=1): mean=%.3f, sd=%.3f, range=[%.3f, %.3f]\n",
      mean(dat$y[r_vec == 1]), sd(dat$y[r_vec == 1]),
      min(dat$y[r_vec == 1]), max(dat$y[r_vec == 1])))

  # Separation check: can a hyperplane perfectly separate r=0 vs r=1 in design space?
  # Simple check: for each r=0 obs, is there always an r=1 obs with similar design matrix values?
  cat(sprintf("  pi at r=0: mean=%.4f, range=[%.6f, %.6f]\n",
      mean(pi2[r_vec == 0]), min(pi2[r_vec == 0]), max(pi2[r_vec == 0])))
  cat(sprintf("  pi at r=1: mean=%.4f, range=[%.6f, %.6f]\n",
      mean(pi2[r_vec == 1]), min(pi2[r_vec == 1]), max(pi2[r_vec == 1])))

  # Confusion matrix at threshold 0.5
  pred_class <- ifelse(pi2 > 0.5, 1, 0)
  acc <- mean(pred_class == r_vec)
  cat(sprintf("  Classification accuracy (pi>0.5): %.3f\n", acc))

  # AUC-like measure
  pi_r1 <- pi2[r_vec == 1]
  pi_r0 <- pi2[r_vec == 0]
  # Approximate AUC: P(pi(r=1) > pi(r=0))
  auc_approx <- mean(outer(pi_r1, pi_r0, ">")) + 0.5 * mean(outer(pi_r1, pi_r0, "=="))
  cat(sprintf("  Approximate AUC: %.4f\n", auc_approx))
}
