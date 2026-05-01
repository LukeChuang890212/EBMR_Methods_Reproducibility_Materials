setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

n_val <- 2000
n_reps <- 200

# Compare: package SE2 (via psi2) vs standalone simplified SE2
# for setting4 model3 with constrained_nr

ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
mu_true <- get_mu_true("setting4")
model_idx <- 3

single_ps <- list(
  formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
  h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]],
  alpha_init.list = list(NULL),
  optimizer = "constrained_nr"
)
W <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

mu_vals <- rep(NA_real_, n_reps)
se_pkg1 <- rep(NA_real_, n_reps)   # Package psi1
se_pkg2 <- rep(NA_real_, n_reps)   # Package psi2
se_simple <- rep(NA_real_, n_reps) # Standalone simplified SE2

for (rep_i in 1:n_reps) {
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W)

    r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
    pi_hat <- ebmr$ps_fit.list[[1]]$fitted.values
    mu_hat <- mean(r_vec * y_vec / pi_hat)
    n <- n_val

    gmm_fit <- ebmr$ps_fit.list[[1]]$gmm_fit
    psi1 <- gmm_fit$psi
    psi2_mat <- gmm_fit$psi2

    design_mat <- ebmr$ps_fit.list[[1]]$design_matrix
    link_type <- ebmr$ps_fit.list[[1]]$link_type
    if (!is.null(link_type) && link_type == "logistic_complement") {
      dot_pi <- -design_mat * (pi_hat * (1 - pi_hat))
    } else {
      dot_pi <- design_mat * (pi_hat * (1 - pi_hat))
    }
    ry_ps_inv2 <- as.vector(r_vec * y_vec * (pi_hat^(-2)))
    H_alpha <- colMeans(dot_pi * ry_ps_inv2)

    # Package SE using psi1
    mu_iid1 <- as.vector(t(r_vec/pi_hat*y_vec) - t(H_alpha) %*% psi1)
    se_pkg1[rep_i] <- sqrt(var(mu_iid1)/n)

    # Package SE using psi2
    if (!is.null(psi2_mat) && is.matrix(psi2_mat) && !any(is.na(psi2_mat))) {
      mu_iid2 <- as.vector(t(r_vec/pi_hat*y_vec) - t(H_alpha) %*% psi2_mat)
      se_pkg2[rep_i] <- sqrt(var(mu_iid2)/n)
    }

    # Standalone simplified SE2 (what we tested before)
    # psi_simple = -(Gamma'WGamma)^{-1} Gamma'W g_i  (simple influence fn for alpha)
    # + CUE correction using simplified formula
    h_alpha_vars <- ps_spec[["h_alpha.list"]][[model_idx]]
    h_full <- cbind(1, as.matrix(dat[, h_alpha_vars, drop=FALSE]))
    p <- ncol(design_mat); h_dim <- ncol(h_full)
    alpha <- gmm_fit$estimates

    Gamma_hat <- gmm_fit$Gamma.hat
    g_mat <- gmm_fit$g.matrix
    W_hat <- gmm_fit$W.hat
    G_vec <- gmm_fit$eta_s

    GtW <- crossprod(Gamma_hat, W_hat)
    GtWG <- GtW %*% Gamma_hat
    H_inv <- tryCatch(solve(GtWG), error = function(e) NULL)
    if (!is.null(H_inv)) {
      psi1_simple <- H_inv %*% GtW %*% t(g_mat)
      dmu <- colMeans(-r_vec*y_vec*(1-pi_hat)/pi_hat * design_mat)

      # Simple CUE correction
      WG <- W_hat %*% G_vec
      gi_WG <- as.vector(g_mat %*% WG)
      Omega <- crossprod(g_mat)/n
      OmWG <- as.vector(Omega %*% WG)
      GtW_g <- GtW %*% t(g_mat)
      GtW_OmWG <- as.vector(GtW %*% OmWG)
      corr <- sweep(GtW_g * rep(gi_WG, each=p), 1, GtW_OmWG, "-")
      psi2_simple <- psi1_simple - H_inv %*% corr
      iid_simple <- r_vec*y_vec/pi_hat - mu_hat - as.vector(t(dmu) %*% psi2_simple)
      se_simple[rep_i] <- sqrt(mean(iid_simple^2)/n)
    }

    mu_vals[rep_i] <- mu_hat
  }, error = function(e) NULL)
}

valid <- !is.na(mu_vals)
v1 <- valid & !is.na(se_pkg1)
v2 <- valid & !is.na(se_pkg2)
vs <- valid & !is.na(se_simple)

esd <- sd(mu_vals[valid])
cat(sprintf("Setting4 model3, constrained_nr:\n"))
cat(sprintf("  Valid: mu=%d, pkg1=%d, pkg2=%d, simple=%d\n", sum(valid), sum(v1), sum(v2), sum(vs)))
cat(sprintf("  Bias=%.4f, ESD=%.4f\n", mean(mu_vals[valid]) - mu_true, esd))
cat(sprintf("  Package SE1 (psi1):     ESE/ESD = %.3f\n", mean(se_pkg1[v1])/sd(mu_vals[v1])))
cat(sprintf("  Package SE2 (psi2):     ESE/ESD = %.3f\n", mean(se_pkg2[v2])/sd(mu_vals[v2])))
cat(sprintf("  Standalone simplified:  ESE/ESD = %.3f\n", mean(se_simple[vs])/sd(mu_vals[vs])))

# Check: are psi1 from package the same as psi1_simple?
# Just check first rep
dat <- all_data[1:n_val, ]
ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W)
gmm_fit <- ebmr$ps_fit.list[[1]]$gmm_fit
psi1_pkg <- gmm_fit$psi  # k x n

h_alpha_vars <- ps_spec[["h_alpha.list"]][[model_idx]]
h_full <- cbind(1, as.matrix(dat[, h_alpha_vars, drop=FALSE]))
Gamma_hat <- gmm_fit$Gamma.hat
g_mat <- gmm_fit$g.matrix
W_hat <- gmm_fit$W.hat
GtW <- crossprod(Gamma_hat, W_hat)
GtWG <- GtW %*% Gamma_hat
H_inv <- solve(GtWG)
psi1_simple <- H_inv %*% GtW %*% t(g_mat)

cat(sprintf("\nPsi1 comparison (rep 1, first 3 obs):\n"))
cat(sprintf("  Package:    %s\n", paste(sprintf("%.6f", psi1_pkg[1,1:3]), collapse=", ")))
cat(sprintf("  Simple:     %s\n", paste(sprintf("%.6f", psi1_simple[1,1:3]), collapse=", ")))
cat(sprintf("  Max diff:   %.2e\n", max(abs(psi1_pkg - psi1_simple))))
