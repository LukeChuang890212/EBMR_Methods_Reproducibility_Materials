setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
  library(EBMRalgorithmFast4)
})

W_func  <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_nu_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1_u2=dat$u1*dat$u2)
n_val   <- 2000

ps_spec <- get_ps_spec("9-alt1")

# ---------------------------------------------------------------
# CASE 1: Setting3, model 2 - inspect H_A vs H_B for one rep
# ---------------------------------------------------------------
cat("=== CASE 1: Setting3, Model 2 ===\n\n")
ps_sub3 <- list(formula.list = ps_spec[["formula.list"]][2],
                h_alpha.list = ps_spec[["h_alpha.list"]][2],
                inv_link     = ps_spec[["inv_link"]],
                outcome      = ps_spec[["outcome"]])
all_data3 <- readRDS("Simulation_Data/setting3.B1_n2000_replicate1000.RDS")

set.seed(12346)  # rep 1
dat3 <- all_data3[1:n_val, ]
ebmr3 <- EBMRAlgorithmFast4$new("y", ps_sub3, dat3, W_func)
res3  <- ebmr3$EBMR_IPW(h_nu=h_nu_fn, type="HT", se.fit=TRUE)
fit3  <- ebmr3$ps_fit.list[[1]]$gmm_fit

cat(sprintf("GMM objective: %.2e\n", fit3$opt$objective))
cat(sprintf("se1 (alpha): %s\n", paste(round(fit3$se1, 4), collapse=", ")))
cat(sprintf("se2 (alpha): %s\n", paste(round(fit3$se2, 4), collapse=", ")))

# Reconstruct H_A (se1) and H_B (se2) manually
n  <- n_val
k  <- length(fit3$estimates)
h  <- nrow(fit3$Gamma.hat)
Gm <- fit3$Gamma.hat   # h x k
W  <- fit3$W.hat       # h x h
G  <- as.vector(fit3$eta_s)  # h = G(theta_hat)
GtW <- crossprod(Gm, W)       # k x h
GtWG <- GtW %*% Gm            # k x k (H_A base)

cat(sprintf("\nG (eta_s) norm: %.2e  (how far from 0)\n", sqrt(sum(G^2))))
cat(sprintf("||Gamma'WGamma|| (k x k): eigenvalues = %s\n",
    paste(round(eigen(GtWG)$values, 4), collapse=", ")))

# What is se2's H_B doing?
# H_B = GtWG + (G'W ⊗ I_k)R - (G'W ⊗ Γ'W)S
# Check the S-correction term magnitude
g_mat <- fit3$g.matrix   # n x h

# Analytical S (no FD): S[:,j] = vec[(1/n) sum_i (dg_j_i g_i' + g_i dg_j_i')]
# Use t_Gamma_arr from the gmm. We can't access it directly; use FD S from se2.
# But instead let's look at what (G'W ⊗ Gamma'W) S actually contributes.

# Reconstruct using se2's FD S_mat
eps2 <- 1e-5
alpha_hat <- fit3$estimates
# Need g function - reconstruct from the gmm_fit internals via re-running
# Instead, directly check psi2 magnitudes

psi1_mat <- fit3$psi   # k x n
psi2_mat <- fit3$psi2  # k x n or NA

cat(sprintf("\nIs psi2 a matrix: %s\n", is.matrix(psi2_mat)))
if (is.matrix(psi2_mat)) {
  cat(sprintf("psi1 column norms: mean=%.4f, max=%.4f\n",
      mean(sqrt(colSums(psi1_mat^2))), max(sqrt(colSums(psi1_mat^2)))))
  cat(sprintf("psi2 column norms: mean=%.4f, max=%.4f\n",
      mean(sqrt(colSums(psi2_mat^2))), max(sqrt(colSums(psi2_mat^2)))))

  # Compare Q_i magnitudes
  # Q_i^(se1) = H_1n + H_2n_2 + H_2n_3 columns
  # We have psi1 = H_A^{-1} (Q^se1)  and  psi2 = H_B^{-1} Q_i^se2
  # H_A psi1 = Q^se1, H_B psi2 = -Q^se2
  # Check: var(psi1)/n should give cov for se1
  v1 <- var(t(psi1_mat)) / n
  v2 <- var(t(psi2_mat)) / n
  cat(sprintf("\nse1 from psi1: %s\n", paste(round(sqrt(diag(v1)), 4), collapse=", ")))
  cat(sprintf("se2 from psi2: %s\n", paste(round(sqrt(diag(v2)), 4), collapse=", ")))
} else {
  cat("psi2 is NA - se2 tryCatch caught an error. Trying to reproduce:\n")
}

# Show se_ipw for this rep
cat(sprintf("\nse_ipw (se1-based): %.4f\n", res3$se_ipw))
cat(sprintf("se2 (alpha se2):    %s\n", paste(round(fit3$se2, 4), collapse=", ")))

# ---------------------------------------------------------------
# CASE 3: Ensemble nu GMM - why psi2 is NA
# ---------------------------------------------------------------
cat("\n\n=== CASE 3: Ensemble nu GMM - diagnose psi2 = NA ===\n\n")
ps_sub4e <- list(formula.list = ps_spec[["formula.list"]][2:3],
                 h_alpha.list = ps_spec[["h_alpha.list"]][2:3],
                 inv_link     = ps_spec[["inv_link"]],
                 outcome      = ps_spec[["outcome"]])
all_data4 <- readRDS("Simulation_Data/setting4.B1_n2000_replicate1000.RDS")

set.seed(12346)
dat4  <- all_data4[1:n_val, ]
ebmr4 <- EBMRAlgorithmFast4$new("y", ps_sub4e, dat4, W_func)
res4  <- ebmr4$EBMR_IPW(h_nu=h_nu_fn, type="HT", se.fit=TRUE)

nu_fit  <- ebmr4$nu_fit
gmm_nu  <- nu_fit$gmm_fit

cat(sprintf("nu GMM objective: %.2e\n", gmm_nu$opt$objective))
cat(sprintf("nu estimates: %s\n", paste(round(gmm_nu$estimates, 4), collapse=", ")))
cat(sprintf("se1 (nu): %s\n", paste(round(gmm_nu$se1, 4), collapse=", ")))
cat(sprintf("psi2 is NA: %s\n", identical(gmm_nu$psi2, NA)))

# Check alpha gmm psi2
cat(sprintf("\nalpha1 psi2 is.matrix: %s\n", is.matrix(ebmr4$ps_fit.list[[1]]$gmm_fit$psi2)))
cat(sprintf("alpha2 psi2 is.matrix: %s\n", is.matrix(ebmr4$ps_fit.list[[2]]$gmm_fit$psi2)))

# ---------------------------------------------------------------
# CASE 2: Setting4, model 3 - se2 = 2x check
# ---------------------------------------------------------------
cat("\n\n=== CASE 2: Setting4, Model 3 (se2 ≈ 2x) ===\n\n")
ps_sub4m3 <- list(formula.list = ps_spec[["formula.list"]][3],
                  h_alpha.list = ps_spec[["h_alpha.list"]][3],
                  inv_link     = ps_spec[["inv_link"]],
                  outcome      = ps_spec[["outcome"]])

set.seed(12346)
dat4m3  <- all_data4[1:n_val, ]
ebmr4m3 <- EBMRAlgorithmFast4$new("y", ps_sub4m3, dat4m3, W_func)
res4m3  <- ebmr4m3$EBMR_IPW(h_nu=h_nu_fn, type="HT", se.fit=TRUE)
fit4m3  <- ebmr4m3$ps_fit.list[[1]]$gmm_fit

cat(sprintf("GMM objective: %.2e\n", fit4m3$opt$objective))
cat(sprintf("G (eta_s) norm: %.2e\n", sqrt(sum(fit4m3$eta_s^2))))
cat(sprintf("se1 (alpha): %s\n", paste(round(fit4m3$se1, 4), collapse=", ")))
cat(sprintf("se2 (alpha): %s\n", paste(round(fit4m3$se2, 4), collapse=", ")))

if (is.matrix(fit4m3$psi2)) {
  v1m3 <- var(t(fit4m3$psi)) / n
  v2m3 <- var(t(fit4m3$psi2)) / n
  cat(sprintf("psi1 norms: mean=%.4f, sd=%.4f\n",
      mean(sqrt(colSums(fit4m3$psi^2))), sd(sqrt(colSums(fit4m3$psi^2)))))
  cat(sprintf("psi2 norms: mean=%.4f, sd=%.4f\n",
      mean(sqrt(colSums(fit4m3$psi2^2))), sd(sqrt(colSums(fit4m3$psi2^2)))))
}
