setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
source("Simulation.r")
library(EBMRalgorithmFast4)

# Setting 3, model 2 alone
data_file <- misspecified_model_all_data_file.list[["setting3"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
n_val <- 2000
mu_true <- get_mu_true("setting3")

ps_spec <- get_ps_spec("9-alt1")
# Model 2 only
single_ps_spec <- list(
  formula.list = ps_spec[["formula.list"]][2],
  h_alpha.list = ps_spec[["h_alpha.list"]][2],
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]]
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

n_reps <- 500
# Store: mu_ipw, mu_ipw.true(NA), se_ipw, se_ipw.true(NA), se_simple, se_H1only, se_noMn, se_noH2n3
results_raw <- matrix(NA, 8, n_reps)
rownames(results_raw) <- c("mu_ipw", "mu_ipw.true", "se_ipw", "se_ipw.true",
                           "se_simple", "se_H1only", "se_noMn", "se_noH2n3")

cat(sprintf("mu_true = %.4f\n", mu_true))
cat("=== RUNNING SE DIAGNOSTIC (500 reps) ===\n")
n_err <- 0

for (rep_i in 1:n_reps) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]
  tryCatch({
    ebmr <- EBMRAlgorithmFast4[["new"]]("y", single_ps_spec, dat, W_func)
    res <- ebmr[["EBMR_IPW"]](
      h_nu = function(dat) cbind(u1 = dat[["u1"]], u2 = dat[["u2"]], z1 = dat[["z1"]], z2 = dat[["z2"]]),
      se.fit = TRUE, type = "HT"
    )

    results_raw[1, rep_i] <- res[["mu_ipw"]]
    results_raw[2, rep_i] <- 0  # placeholder for mu_ipw.true (not used)
    results_raw[3, rep_i] <- res[["se_ipw"]]
    results_raw[4, rep_i] <- 0  # placeholder for se_ipw.true (not used)

    # Extract GMM quantities for SE variants
    gmm_fit <- ebmr[["ps_fit.list"]][[1]][["gmm_fit"]]
    Gamma.hat <- gmm_fit[["Gamma.hat"]]
    W.hat <- gmm_fit[["W.hat"]]
    g.matrix <- gmm_fit[["g.matrix"]]
    psi_alpha <- gmm_fit[["psi"]]
    eta_s <- gmm_fit[["eta_s"]]
    Q_full <- gmm_fit[["Q"]]

    n <- nrow(g.matrix)
    esteq_dim <- ncol(g.matrix)
    param_dim <- nrow(psi_alpha)

    GtW <- crossprod(Gamma.hat, W.hat)
    GtWG <- GtW %*% Gamma.hat
    H_0n <- tryCatch(-solve(GtWG), error = function(e) -solve(GtWG + 1e-8 * diag(nrow(GtWG))))
    t_g <- t(g.matrix)
    H_1n <- GtW %*% (t_g - eta_s)

    # Delta method: use the se_ipw influence function structure
    # mu_iid = base_term - H_alpha.w' %*% psi_alpha (from EBMR_IPW)
    # To compute SE variants, we replace psi_alpha with different versions
    # and recompute se = sqrt(var(mu_iid)/n)

    # Reconstruct base_term and H_alpha.w from the EBMR_IPW internals
    ps_fit <- ebmr[["ps_fit.list"]][[1]]
    pi_hat <- ps_fit[["fitted.values"]]
    r_vec <- dat[["r"]]
    y_vec <- dat[["y"]]
    design_mat <- ps_fit[["design_matrix"]]
    link_type <- ps_fit[["link_type"]]

    # Mixture parameters
    eps_mix <- 0.01
    scale_mix <- 0.99
    # Recover f from pi: f = (pi - eps) / scale
    f_hat <- (pi_hat - eps_mix) / scale_mix

    ry_ps_inv2 <- as.vector(r_vec * y_vec * (pi_hat^(-2)))
    # dot_pi = scale * f'(eta) * X  (matching the mixture model)
    if (link_type == "logistic_complement") {
      dot_pi <- -scale_mix * design_mat * (f_hat * (1 - f_hat))
    } else if (link_type == "logistic") {
      dot_pi <- scale_mix * design_mat * (f_hat * (1 - f_hat))
    }
    H_alpha.w <- colMeans(dot_pi * ry_ps_inv2)
    base_term <- as.vector(r_vec / pi_hat * y_vec)

    # SE variant: Simple GMM
    psi_simple <- H_0n %*% GtW %*% t_g
    mu_iid_simple <- base_term - as.vector(t(H_alpha.w) %*% psi_simple)
    results_raw[5, rep_i] <- sqrt(var(mu_iid_simple) / n)

    # SE variant: H_1n only (with M_n in bread)
    psi_H1only <- Q_full %*% H_1n
    mu_iid_H1 <- base_term - as.vector(t(H_alpha.w) %*% psi_H1only)
    results_raw[6, rep_i] <- sqrt(var(mu_iid_H1) / n)

    # Recover H_2n_2 and H_2n_3
    W_eta <- as.vector(W.hat %*% eta_s)
    tGamma_W_eta <- as.vector(GtW %*% eta_s)
    GtW_g <- GtW %*% t_g
    g_W_eta <- as.vector(g.matrix %*% W_eta)
    H_2n_3 <- matrix(tGamma_W_eta, param_dim, n) - GtW_g * matrix(g_W_eta, param_dim, n, byrow = TRUE)
    Q_inv <- solve(Q_full)
    H_sum <- Q_inv %*% psi_alpha
    H_2n_2 <- H_sum - H_1n - H_2n_3

    # SE variant: No M_n
    psi_noMn <- H_0n %*% (H_1n + H_2n_2 + H_2n_3)
    mu_iid_noMn <- base_term - as.vector(t(H_alpha.w) %*% psi_noMn)
    results_raw[7, rep_i] <- sqrt(var(mu_iid_noMn) / n)

    # SE variant: No H_2n_3
    psi_noH2n3 <- Q_full %*% (H_1n + H_2n_2)
    mu_iid_noH2n3 <- base_term - as.vector(t(H_alpha.w) %*% psi_noH2n3)
    results_raw[8, rep_i] <- sqrt(var(mu_iid_noH2n3) / n)

  }, error = function(e) {
    n_err <<- n_err + 1
    if (n_err <= 5) cat(sprintf("  Error rep %d: %s\n", rep_i, conditionMessage(e)))
  })
  if (rep_i %% 100 == 0) cat(sprintf("  %d / %d done (errors: %d)\n", rep_i, n_reps, n_err))
}

cat(sprintf("\nErrors: %d / %d\n", n_err, n_reps))

# --- Clean with clean_sim_result ---
cat("\n========== BEFORE CLEANING ==========\n")
valid_all <- !is.na(results_raw[1, ])
mu_all <- results_raw[1, valid_all]
esd_all <- sd(mu_all)
cat(sprintf("Valid: %d, Bias: %.4f, ESD: %.4f\n", sum(valid_all), mean(mu_all) - mu_true, esd_all))
cat(sprintf("%-15s %10s %10s %10s\n", "SE variant", "Mean ESE", "ESE/ESD", "Median ESE"))
for (nm in c("se_ipw", "se_simple", "se_H1only", "se_noMn", "se_noH2n3")) {
  idx <- which(rownames(results_raw) == nm)
  se_vals <- results_raw[idx, valid_all]
  cat(sprintf("%-15s %10.4f %10.4f %10.4f\n", nm, mean(se_vals), mean(se_vals)/esd_all, median(se_vals)))
}

cat("\n========== AFTER clean_sim_result ==========\n")
# Apply clean_sim_result to the first 4 rows (mu_ipw, mu_ipw.true, se_ipw, se_ipw.true)
cleaned <- clean_sim_result(results_raw[1:4, ], multiplier = 3, verbose = TRUE)
keep_mask <- !is.na(cleaned$result[1, ])
# But we need to map back to original columns — clean_sim_result removes columns
# So let's replicate its logic to get the kept indices
na_mask <- apply(results_raw[1:4, , drop = FALSE], 2, function(x) !any(is.na(x)))
results_nona <- results_raw[, na_mask, drop = FALSE]
mu_ipw <- results_nona[1, ]
Q1 <- quantile(mu_ipw, 0.25)
Q3 <- quantile(mu_ipw, 0.75)
IQR_value <- Q3 - Q1
is_outlier <- mu_ipw < (Q1 - 3 * IQR_value) | mu_ipw > (Q3 + 3 * IQR_value)
max_remove <- floor(0.01 * length(mu_ipw))
if (sum(is_outlier) > max_remove && max_remove > 0) {
  dist_from_median <- abs(mu_ipw - median(mu_ipw))
  outlier_idx <- which(is_outlier)
  keep_idx <- outlier_idx[order(dist_from_median[outlier_idx], decreasing = TRUE)[(max_remove + 1):length(outlier_idx)]]
  is_outlier[keep_idx] <- FALSE
}
results_clean <- results_nona[, !is_outlier, drop = FALSE]

cat(sprintf("Kept: %d (removed %d NA, %d outlier)\n",
    ncol(results_clean), sum(!na_mask), sum(is_outlier)))

mu_clean <- results_clean[1, ]
esd_clean <- sd(mu_clean)
cat(sprintf("Bias: %.4f, ESD: %.4f\n", mean(mu_clean) - mu_true, esd_clean))

cat(sprintf("\n%-15s %10s %10s %10s %10s\n", "SE variant", "Mean ESE", "ESE/ESD", "Median ESE", "Med/ESD"))
for (nm in c("se_ipw", "se_simple", "se_H1only", "se_noMn", "se_noH2n3")) {
  idx <- which(rownames(results_clean) == nm)
  se_vals <- results_clean[idx, ]
  cat(sprintf("%-15s %10.4f %10.4f %10.4f %10.4f\n",
      nm, mean(se_vals), mean(se_vals)/esd_clean, median(se_vals), median(se_vals)/esd_clean))
}

# Coverage probability
ci_l <- results_clean[1, ] - 1.96 * results_clean[3, ]
ci_u <- results_clean[1, ] + 1.96 * results_clean[3, ]
cp <- mean(mu_true >= ci_l & mu_true <= ci_u)
cat(sprintf("\nCP (se_ipw): %.4f\n", cp))

# Coverage with se_noH2n3
ci_l2 <- results_clean[1, ] - 1.96 * results_clean[8, ]
ci_u2 <- results_clean[1, ] + 1.96 * results_clean[8, ]
cp2 <- mean(mu_true >= ci_l2 & mu_true <= ci_u2)
cat(sprintf("CP (se_noH2n3): %.4f\n", cp2))
