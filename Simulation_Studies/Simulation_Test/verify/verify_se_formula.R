## Verify the SE formula for mu_001 by comparing analytical SE vs numerical differentiation
## Focus on a REP where M3 converges well (not rep 1)
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)
library(numDeriv)

ps_spec_base <- get_ps_spec("9-alt1")
h_alpha_fn <- function(dat) cbind(
  u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2,
  u1u2 = dat$u1*dat$u2, z1z2 = dat$z1*dat$z2,
  u1z1 = dat$u1*dat$z1, z1u2 = dat$z1*dat$u2,
  u1z2 = dat$u1*dat$z2, u2z2 = dat$u2*dat$z2)
ps_spec <- list(
  formula.list = ps_spec_base$formula.list,
  h_alpha.list = list(h_alpha_fn, h_alpha_fn, h_alpha_fn),
  inv_link = ps_spec_base$inv_link,
  outcome = ps_spec_base$outcome
)
W_fn <- function(g.matrix) solve(var(g.matrix))
h_nu_fn <- h_alpha_fn

data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]]
all_data <- NULL
for (fi in seq_along(data_file)) {
  if (file.exists(data_file[[fi]])) {
    test <- readRDS(data_file[[fi]])
    if (nrow(test) / 1000 == 2000) { all_data <- test; break }
  }
}

# Use rep 2 which had reasonable SE
cat("=== Verifying SE formula for mu_001, rep 2 ===\n\n")
dat <- all_data[2001:4000, ]
n <- 2000

ebmr <- EBMRAlgorithmFast4$new("y", ps_spec, dat, W_fn)
ipw <- ebmr$EBMR_IPW(h_nu_fn, model_indices = 3, se.fit = TRUE)
cat(sprintf("Package SE: %.6f\n", ipw$se_ipw))
cat(sprintf("mu_ipw: %.6f\n\n", ipw$mu_ipw))

# Now manually compute SE using the influence function approach
# For single model (J=1), no ensemble step. SE comes purely from alpha GMM.
ps_fit3 <- ebmr$ps_fit.list[[3]]
alpha_hat <- ps_fit3$coefficients
gmm_fit <- ps_fit3$gmm_fit

cat("--- GMM fit details ---\n")
cat(sprintf("converged: %s, grad_norm: %.2e\n", gmm_fit$opt$converged, gmm_fit$opt$final_grad_norm))
cat(sprintf("alpha SE from GMM: %s\n\n", paste(round(gmm_fit$se, 6), collapse=", ")))

# Extract quantities from the GMM
g_mat <- gmm_fit$g.matrix      # n x h_dim
W_hat <- gmm_fit$W.hat         # h_dim x h_dim
Gamma_hat <- gmm_fit$Gamma.hat # h_dim x k
psi_alpha <- gmm_fit$psi       # k x n (influence function for alpha)

h_dim <- ncol(g_mat)
k <- length(alpha_hat)

cat(sprintf("Dimensions: h_dim=%d, k=%d, n=%d\n", h_dim, k, n))
cat(sprintf("dim(psi_alpha) = %d x %d\n", nrow(psi_alpha), ncol(psi_alpha)))
cat(sprintf("dim(Gamma_hat) = %d x %d\n", nrow(Gamma_hat), ncol(Gamma_hat)))
cat(sprintf("dim(W_hat) = %d x %d\n\n", nrow(W_hat), ncol(W_hat)))

# Verify: SE from psi should match gmm_fit$se
se_from_psi <- sqrt(diag(var(t(psi_alpha)) / n))
cat(sprintf("SE from psi: %s\n", paste(round(se_from_psi, 6), collapse=", ")))
cat(sprintf("SE from gmm: %s\n\n", paste(round(gmm_fit$se, 6), collapse=", ")))

# Now check the EBMR_IPW SE formula
# For J=1 (single model), the IPW SE formula uses:
# mu_ipw.iid = r/pi*y - H_alpha.w' %*% psi_alpha
# where H_alpha.w = E[r*y/pi^2 * dot_pi * w]  (but w=1 for J=1)

# Let's verify by numerical differentiation
# The IPW estimator as a function of alpha:
design_mat <- ps_fit3$design_matrix
r_vec <- as.vector(dat$r)
y_vec <- as.vector(dat$y)

mu_ipw_fn <- function(alpha) {
  eta <- as.vector(design_mat %*% alpha)
  eta <- pmin(pmax(eta, -20), 20)
  pi_vec <- 1 / (1 + exp(eta))  # logistic complement
  mean(r_vec * y_vec / pi_vec)
}

# Numerical gradient of mu_ipw w.r.t. alpha
dmu_dalpha_num <- numDeriv::grad(mu_ipw_fn, alpha_hat)
cat("dmu/dalpha (numerical):", round(dmu_dalpha_num, 6), "\n")

# The analytical version from the package:
# H_alpha.w = colMeans(r*y/pi^2 * dot_pi)
# For single model with logistic_complement: dot_pi = -design_mat * pi*(1-pi)
pi_hat <- ps_fit3$fitted.values
ry_pi2 <- as.vector(r_vec * y_vec / (pi_hat^2))
dot_pi <- -design_mat * (pi_hat * (1 - pi_hat))
H_alpha_w <- colMeans(ry_pi2 * dot_pi)
cat("H_alpha.w (analytical):", round(H_alpha_w, 6), "\n")
cat("These should match (H_alpha.w = -dmu/dalpha):", round(-H_alpha_w, 6), "\n\n")

# Now the full IPW SE:
# mu_ipw.iid[i] = r_i*y_i/pi_i - H_alpha.w' * psi_alpha[,i]
mu_base <- as.vector(r_vec * y_vec / pi_hat)
mu_iid_manual <- mu_base - as.vector(t(H_alpha_w) %*% psi_alpha)
se_manual <- sqrt(var(mu_iid_manual) / n)
cat(sprintf("SE manual (from influence function): %.6f\n", se_manual))
cat(sprintf("SE package: %.6f\n\n", ipw$se_ipw))

# Check: what if we DON'T account for alpha estimation uncertainty?
se_naive <- sqrt(var(mu_base) / n)
cat(sprintf("SE naive (no alpha uncertainty): %.6f\n", se_naive))

# Now do the same for REP 1 (the problematic one)
cat("\n\n=== Same check for rep 1 (problematic) ===\n\n")
dat1 <- all_data[1:2000, ]
ebmr1 <- EBMRAlgorithmFast4$new("y", ps_spec, dat1, W_fn)
ipw1 <- ebmr1$EBMR_IPW(h_nu_fn, model_indices = 3, se.fit = TRUE)
cat(sprintf("Package SE: %.6f\n", ipw1$se_ipw))

ps_fit3_1 <- ebmr1$ps_fit.list[[3]]
gmm_fit1 <- ps_fit3_1$gmm_fit
psi1 <- gmm_fit1$psi

cat(sprintf("converged: %s, grad_norm: %.2e\n", gmm_fit1$opt$converged, gmm_fit1$opt$final_grad_norm))
cat(sprintf("max|psi|: %.4f\n", max(abs(psi1))))
cat(sprintf("alpha SE from GMM: %s\n\n", paste(round(gmm_fit1$se, 6), collapse=", ")))

# Manual IPW SE for rep 1
pi1 <- ps_fit3_1$fitted.values
mu_base1 <- as.vector(dat1$r * dat1$y / pi1)
ry_pi2_1 <- as.vector(dat1$r * dat1$y / (pi1^2))
dot_pi1 <- -ps_fit3_1$design_matrix * (pi1 * (1 - pi1))
H_alpha_w1 <- colMeans(ry_pi2_1 * dot_pi1)

mu_iid1 <- mu_base1 - as.vector(t(H_alpha_w1) %*% psi1)
se_manual1 <- sqrt(var(mu_iid1) / n)
se_naive1 <- sqrt(var(mu_base1) / n)

cat(sprintf("SE manual (influence fn): %.6f\n", se_manual1))
cat(sprintf("SE naive (no alpha uncert): %.6f\n", se_naive1))
cat(sprintf("SE package: %.6f\n\n", ipw1$se_ipw))

# Decompose: what's the variance contribution from each term?
var_base1 <- var(mu_base1)
var_correction1 <- var(as.vector(t(H_alpha_w1) %*% psi1))
cov_term1 <- 2 * cov(mu_base1, as.vector(t(H_alpha_w1) %*% psi1))
cat("Variance decomposition (rep 1):\n")
cat(sprintf("  var(base) = %.6f\n", var_base1))
cat(sprintf("  var(H'psi) = %.6f\n", var_correction1))
cat(sprintf("  -2*cov(base, H'psi) = %.6f\n", -cov_term1))
cat(sprintf("  total var = %.6f\n", var_base1 + var_correction1 - cov_term1))
cat(sprintf("  SE = sqrt(total/n) = %.6f\n", sqrt((var_base1 + var_correction1 - cov_term1)/n)))
