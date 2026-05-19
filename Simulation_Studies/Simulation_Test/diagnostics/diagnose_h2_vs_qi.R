setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

n_val <- 2000
n_reps <- 200

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
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

mu_vals <- rep(NA_real_, n_reps)
# 4 variants:
# A: psi = -H2^{-1} Q_full  (package SE2)
# B: psi = -(GtWG)^{-1} Q_full  (simple Hessian + full Q)
# C: psi = -H2^{-1} (GtW g_i)  (full Hessian + simple Q)
# D: psi = -(GtWG)^{-1} GtW g_i  (fully simple)
se_A <- se_B <- se_C <- se_D <- rep(NA_real_, n_reps)

for (rep_i in 1:n_reps) {
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
    r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
    pi_hat <- ebmr$ps_fit.list[[1]]$fitted.values
    mu_hat <- mean(r_vec * y_vec / pi_hat)
    n <- n_val

    gmm_fit <- ebmr$ps_fit.list[[1]]$gmm_fit
    design_mat <- ebmr$ps_fit.list[[1]]$design_matrix
    link_type <- ebmr$ps_fit.list[[1]]$link_type

    if (!is.null(link_type) && link_type == "logistic_complement") {
      dot_pi <- -design_mat * (pi_hat * (1 - pi_hat))
    } else {
      dot_pi <- design_mat * (pi_hat * (1 - pi_hat))
    }
    ry_ps_inv2 <- as.vector(r_vec * y_vec * (pi_hat^(-2)))
    H_alpha <- colMeans(dot_pi * ry_ps_inv2)

    Gamma_hat <- gmm_fit$Gamma.hat
    g_mat <- gmm_fit$g.matrix
    W_hat <- gmm_fit$W.hat
    eta_s <- gmm_fit$eta_s
    GtW <- crossprod(Gamma_hat, W_hat)
    GtWG <- GtW %*% Gamma_hat
    GtWG_inv <- tryCatch(solve(GtWG), error = function(e) NULL)

    # Package psi2 (variant A) — directly from package
    psi2_pkg <- gmm_fit$psi2
    if (!is.null(psi2_pkg) && is.matrix(psi2_pkg) && !any(is.na(psi2_pkg))) {
      mu_iid_A <- as.vector(t(r_vec/pi_hat*y_vec) - t(H_alpha) %*% psi2_pkg)
      se_A[rep_i] <- sqrt(var(mu_iid_A)/n)
    }

    if (!is.null(GtWG_inv)) {
      # Reconstruct Q_full: need Gamma_i, W, G
      WG <- as.vector(W_hat %*% eta_s)
      GtW_g <- GtW %*% t(g_mat)  # k x n
      g_WG <- as.vector(g_mat %*% WG)  # n-vector

      # Gamma_i'WG: need individual Gamma_i
      # Use dg function to get t_Gamma_arr
      # Instead, use the package's t_Gamma construction
      # We need: for each i, Gamma_i' W G
      # Gamma_i[l, j] = dg_il/dalpha_j
      # Using the analytical formula
      if (link_type == "logistic_complement") {
        cf <- r_vec * (1 - pi_hat) / pi_hat
      } else {
        cf <- -r_vec * (1 - pi_hat) / pi_hat
      }
      h_alpha_vars <- ps_spec[["h_alpha.list"]][[model_idx]]
      h_x <- cbind(1, as.matrix(dat[, h_alpha_vars, drop = FALSE]))
      p <- ncol(design_mat)
      h <- ncol(h_x)

      # Gamma_i' W G for each i: sum over l of Gamma_i[l,j] * WG[l]
      # Gamma_i[l, j] = cf_i * h_x[i,l] * X[i,j]
      # sum_l Gamma_i[l,j] * WG[l] = cf_i * X[i,j] * sum_l h_x[i,l] * WG[l]
      hx_WG <- as.vector(h_x %*% WG)  # n-vector: sum_l h_x[i,l] * WG[l]
      cf_hxWG <- cf * hx_WG  # n-vector
      GammaT_WG <- t(design_mat * cf_hxWG)  # p x n

      Q_full <- GtW_g + GammaT_WG - sweep(GtW_g, 2, g_WG, `*`)

      # Simple Q: just GtW g_i
      Q_simple <- GtW %*% t(g_mat)  # k x n

      # Variant B: -(GtWG)^{-1} Q_full
      psi_B <- -(GtWG_inv %*% Q_full)
      mu_iid_B <- as.vector(t(r_vec/pi_hat*y_vec) - t(H_alpha) %*% psi_B)
      se_B[rep_i] <- sqrt(var(mu_iid_B)/n)

      # Variant C: need H2. We'll extract it from the package via psi2 and Q_full
      # Actually, we can compute: psi2_pkg = -H2^{-1} Q_full
      # So H2^{-1} Q_full = -psi2_pkg
      # Then H2^{-1} Q_simple = H2^{-1} (GtW g_i)
      # We can't easily extract H2, but we can compute H2^{-1} via the relationship
      # Actually, let me just compute H2 from the formula

      # Skip for now — focus on B and D
      # Variant D: simple psi = -(GtWG)^{-1} GtW g_i
      psi_D <- -(GtWG_inv %*% Q_simple)  # = -(GtWG)^{-1} GtW g_i
      dmu <- colMeans(-r_vec*y_vec*(1-pi_hat)/pi_hat * design_mat)
      iid_D <- r_vec*y_vec/pi_hat - mu_hat - as.vector(t(dmu) %*% (-psi_D))
      # Note: psi_D has negative sign convention, dmu = H_alpha
      # iid = r*y/pi - mu + H_alpha' * psi_D (where psi_D already has minus)
      # Actually: mu_iid_D = r/pi*y - H_alpha' * psi_D
      mu_iid_D <- as.vector(t(r_vec/pi_hat*y_vec) - t(H_alpha) %*% psi_D)
      se_D[rep_i] <- sqrt(var(mu_iid_D)/n)
    }

    mu_vals[rep_i] <- mu_hat
  }, error = function(e) { cat(sprintf("Rep %d error: %s\n", rep_i, e$message)) })
}

valid <- !is.na(mu_vals)
esd <- sd(mu_vals[valid])

cat("=== SE decomposition: Model 3, constrained_nr ===\n\n")
cat(sprintf("Bias=%.4f, ESD=%.4f\n\n", mean(mu_vals[valid]) - mu_true, esd))

for (nm in c("A", "B", "C", "D")) {
  sv <- get(paste0("se_", nm))
  v <- valid & !is.na(sv)
  label <- c(A="H2^{-1} Q_full (package)", B="(GtWG)^{-1} Q_full",
             C="H2^{-1} GtW_g (hybrid)", D="(GtWG)^{-1} GtW_g (simple)")[nm]
  if (sum(v) > 0) {
    cat(sprintf("  %-30s: ESE/ESD=%.3f  (valid=%d)\n", label, mean(sv[v])/sd(mu_vals[v]), sum(v)))
  }
}
