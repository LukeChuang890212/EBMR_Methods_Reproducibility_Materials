## Diagnose why mu_001 SE explodes for rep 1 with W=solve(var(g))
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

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
W_var <- function(g.matrix) solve(var(g.matrix))
W_sm <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_nu_fn <- h_alpha_fn

data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]]
all_data <- NULL
for (fi in seq_along(data_file)) {
  if (file.exists(data_file[[fi]])) {
    test <- readRDS(data_file[[fi]])
    if (nrow(test) / 1000 == 2000) { all_data <- test; break }
  }
}

# Rep 1 is the problematic one
dat <- all_data[1:2000, ]

cat("=== Rep 1 diagnosis for mu_001 ===\n\n")

# Fit with var(g) W
ebmr_var <- EBMRAlgorithmFast4$new("y", ps_spec, dat, W_var)
ps_fit3_var <- ebmr_var$ps_fit.list[[3]]

# Fit with second moment W
ebmr_sm <- EBMRAlgorithmFast4$new("y", ps_spec, dat, W_sm)
ps_fit3_sm <- ebmr_sm$ps_fit.list[[3]]

cat("Alpha estimates (M3):\n")
cat("  var(g) W:", round(ps_fit3_var$coefficients, 4), "\n")
cat("  sm W:    ", round(ps_fit3_sm$coefficients, 4), "\n\n")

cat("PS distribution (M3):\n")
ps_var <- ps_fit3_var$fitted.values
ps_sm <- ps_fit3_sm$fitted.values
cat(sprintf("  var(g) W: min=%.4f, Q1=%.4f, med=%.4f, Q3=%.4f, max=%.4f\n",
    min(ps_var), quantile(ps_var,0.25), median(ps_var), quantile(ps_var,0.75), max(ps_var)))
cat(sprintf("  sm W:     min=%.4f, Q1=%.4f, med=%.4f, Q3=%.4f, max=%.4f\n\n",
    min(ps_sm), quantile(ps_sm,0.25), median(ps_sm), quantile(ps_sm,0.75), max(ps_sm)))

cat("GMM convergence (M3):\n")
cat(sprintf("  var(g) W: grad_norm=%.2e, obj=%.2e, cond=%.2e\n",
    ps_fit3_var$gmm_fit$opt$final_grad_norm, ps_fit3_var$gmm_fit$opt$objective,
    ps_fit3_var$gmm_fit$opt$solution_cond))
cat(sprintf("  sm W:     grad_norm=%.2e, obj=%.2e, cond=%.2e\n\n",
    ps_fit3_sm$gmm_fit$opt$final_grad_norm, ps_fit3_sm$gmm_fit$opt$objective,
    ps_fit3_sm$gmm_fit$opt$solution_cond))

# Check g matrix and W condition
g_mat_var <- ps_fit3_var$gmm_fit$g.matrix
g_mat_sm <- ps_fit3_sm$gmm_fit$g.matrix
W_hat_var <- ps_fit3_var$gmm_fit$W.hat
W_hat_sm <- ps_fit3_sm$gmm_fit$W.hat

cat("W matrix condition number:\n")
cat(sprintf("  var(g) W: cond(W)=%.2e\n", kappa(W_hat_var)))
cat(sprintf("  sm W:     cond(W)=%.2e\n\n", kappa(W_hat_sm)))

cat("var(g) vs g'g/n for M3 alpha g matrix:\n")
var_g <- var(g_mat_var)
gg_n <- crossprod(g_mat_var) / nrow(g_mat_var)
cat(sprintf("  cond(var(g))=%.2e\n", kappa(var_g)))
cat(sprintf("  cond(g'g/n)=%.2e\n\n", kappa(gg_n)))

# Check G (mean of g) — should be ~0 at solution
G_var <- colMeans(g_mat_var)
G_sm <- colMeans(g_mat_sm)
cat("G = mean(g) at solution:\n")
cat(sprintf("  var(g) W: max|G|=%.2e, norm(G)=%.2e\n", max(abs(G_var)), sqrt(sum(G_var^2))))
cat(sprintf("  sm W:     max|G|=%.2e, norm(G)=%.2e\n\n", max(abs(G_sm)), sqrt(sum(G_sm^2))))

# SE details
cat("SE from GMM fit (alpha SE):\n")
cat(sprintf("  var(g) W se: %s\n", paste(round(ps_fit3_var$gmm_fit$se, 6), collapse=", ")))
cat(sprintf("  sm W se:     %s\n\n", paste(round(ps_fit3_sm$gmm_fit$se, 6), collapse=", ")))

# Check psi (influence function) for alpha
psi_var <- ps_fit3_var$gmm_fit$psi
psi_sm <- ps_fit3_sm$gmm_fit$psi
cat("Influence function psi (alpha) diagnostics:\n")
cat(sprintf("  var(g) W: max|psi|=%.4f, max row var=%.6f\n",
    max(abs(psi_var)), max(apply(psi_var, 1, var))))
cat(sprintf("  sm W:     max|psi|=%.4f, max row var=%.6f\n\n",
    max(abs(psi_sm)), max(apply(psi_sm, 1, var))))
