setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})

n_val <- 2000
n_reps <- 200
ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
mu_true <- get_mu_true("setting4")

formula_j <- ps_spec[["formula.list"]][[3]]
h_alpha_vars <- ps_spec[["h_alpha.list"]][[3]]

mu_vals <- rep(NA_real_, n_reps)
se1_vals <- rep(NA_real_, n_reps)
se2_vals <- rep(NA_real_, n_reps)
alpha_mat <- matrix(NA, n_reps, 4)
conv_vals <- rep(NA, n_reps)
grad_vals <- rep(NA_real_, n_reps)

for (rep_i in 1:n_reps) {
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
  tryCatch({
    r_vec <- dat[["r"]]
    y_vec <- dat[["y"]]
    design_mat <- model.matrix(formula_j, data=dat)
    h_alpha_mat <- as.matrix(dat[, h_alpha_vars, drop=FALSE])
    h_full <- cbind(1, h_alpha_mat)
    n <- n_val
    p <- ncol(design_mat)
    h_dim <- ncol(h_full)

    model_fn <- function(alpha) {
      eta <- as.vector(design_mat %*% alpha)
      1 / (1 + exp(-eta))
    }

    Phi_alpha <- function(alpha) {
      pi_hat <- model_fn(alpha)
      (r_vec / pi_hat - 1) * h_full
    }

    Gamma_fn <- function(alpha) {
      pi_hat <- model_fn(alpha)
      common <- -r_vec * (1 - pi_hat) / pi_hat
      crossprod(h_full * common, design_mat) / n
    }

    W_func <- function(g_mat) solve(crossprod(g_mat) / nrow(g_mat))

    bnd <- 5

    # Step 1: W=I with bounds
    obj_I <- function(alpha) {
      G <- colMeans(Phi_alpha(alpha))
      sum(G^2)
    }
    opt0 <- optim(rep(0, p), obj_I, method = "L-BFGS-B",
                  lower = rep(-bnd, p), upper = rep(bnd, p),
                  control = list(maxit = 5000))
    estimates <- opt0$par

    # Step 2+: iterative GMM with bounds
    best_grad <- Inf
    best_est <- estimates
    no_improve <- 0L
    converged <- FALSE

    for (t in 1:200) {
      g_mat <- Phi_alpha(estimates)
      G_vec <- colMeans(g_mat)
      W_hat <- tryCatch(W_func(g_mat), error = function(e) diag(h_dim))

      Gamma_hat <- Gamma_fn(estimates)
      conv_grad <- 2 * as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))
      grad_norm <- max(abs(conv_grad))

      if (grad_norm < 1e-6) { converged <- TRUE; break }

      if (grad_norm < best_grad) {
        best_grad <- grad_norm
        best_est <- estimates
        no_improve <- 0L
      } else {
        no_improve <- no_improve + 1L
        if (no_improve >= 20L) { estimates <- best_est; break }
      }

      obj_fn <- function(alpha) {
        G <- colMeans(Phi_alpha(alpha))
        as.numeric(t(G) %*% W_hat %*% G)
      }
      opt <- optim(estimates, obj_fn, method = "L-BFGS-B",
                   lower = rep(-bnd, p), upper = rep(bnd, p),
                   control = list(maxit = 1000))
      estimates <- opt$par
    }

    # === Compute mu and SE ===
    pi_hat <- model_fn(estimates)
    mu_hat <- mean(r_vec * y_vec / pi_hat)

    g_mat <- Phi_alpha(estimates)
    W_hat <- tryCatch(W_func(g_mat), error = function(e) diag(h_dim))
    Gamma_hat <- Gamma_fn(estimates)

    # H = Gamma' W Gamma
    H_mat <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
    H_inv <- tryCatch(solve(H_mat), error = function(e) NULL)

    if (!is.null(H_inv)) {
      # psi1 (standard GMM influence): H^{-1} Gamma' W g_i
      # psi1 is p x n
      psi1 <- H_inv %*% crossprod(Gamma_hat, W_hat) %*% t(g_mat)  # p x n

      # d(mu)/d(alpha): derivative of mean(r*y/pi) w.r.t. alpha
      dmu_dalpha <- colMeans(-r_vec * y_vec * (1 - pi_hat) / pi_hat * design_mat)

      # SE1: standard influence function
      # iid1_i = r_i*y_i/pi_i - mu_hat - dmu_dalpha' psi1_i
      iid1 <- r_vec * y_vec / pi_hat - mu_hat - as.vector(t(dmu_dalpha) %*% psi1)
      se1_vals[rep_i] <- sqrt(mean(iid1^2) / n)

      # SE2: use psi2 (CUE-corrected influence)
      # psi2 = H^{-1} Gamma' W g_i - H^{-1} [sum_l (e_l' H^{-1} Gamma' W) d2g_l' W g_i] ?
      # Simpler: for overidentified GMM with optimal W,
      # the CUE correction adds terms from dW/dalpha
      # psi2_i = psi1_i + H^{-1} * correction_i
      #
      # The correction for each i involves:
      # sum over moments l: d(W)/d(alpha) terms
      # This is complex. Use numerical approach:
      # psi2_i = -H^{-1} * [Gamma' W g_i + 0.5 * d(g_i' W g_i)/dalpha]
      # where the second term captures the W(alpha) dependence
      #
      # Actually, the standard formula (Newey & McFadden 1994):
      # For CUE: psi2 = H_cue^{-1} * Gamma' * Omega^{-1} * g_i
      # where H_cue = Gamma' Omega^{-1} Gamma (same as H when W=Omega^{-1})
      # Plus correction: - H^{-1} * P_i where
      # P_i = sum_k sum_l [H^{-1}]_{.,k} * [Gamma' W]_{k,.} * (g_i g_i'/n - Omega) * W * G
      #
      # For a simpler approach that matches the package:
      # Omega = g'g/n, W = Omega^{-1}
      # At optimal W, Omega^{-1} = W, so:
      # D_i = (g_i g_i' - Omega) * W * G_vec  (h_dim vector for each i)
      # correction_i = Gamma' * W * D_i  (p vector)
      # psi2_i = psi1_i - H^{-1} * correction_i

      G_vec <- colMeans(g_mat)
      Omega <- crossprod(g_mat) / n
      WG <- W_hat %*% G_vec

      # For each i: D_i = (g_i * g_i' - Omega) %*% WG
      # g_i is h_dim vector, g_i * g_i' is h_dim x h_dim
      # D_i = g_i * (g_i' WG) - Omega %*% WG
      gi_dot_WG <- as.vector(g_mat %*% WG)  # n-vector: g_i' * WG for each i
      OmegaWG <- as.vector(Omega %*% WG)     # h_dim vector (constant)

      # D_i = g_i * gi_dot_WG[i] - OmegaWG  (h_dim vector for each i)
      # correction_i = t(Gamma) %*% W %*% D_i
      GtW <- crossprod(Gamma_hat, W_hat)  # p x h_dim

      # correction for all i: GtW %*% (g_i * scalar - constant)
      # = scalar_i * (GtW %*% g_i) - GtW %*% OmegaWG
      GtW_g <- GtW %*% t(g_mat)  # p x n
      GtW_OmegaWG <- as.vector(GtW %*% OmegaWG)  # p vector (constant)

      # correction_i = gi_dot_WG[i] * GtW_g[,i] - GtW_OmegaWG
      correction <- sweep(GtW_g * rep(gi_dot_WG, each = p), 1, GtW_OmegaWG, "-")  # p x n

      psi2 <- psi1 - H_inv %*% correction  # p x n

      iid2 <- r_vec * y_vec / pi_hat - mu_hat - as.vector(t(dmu_dalpha) %*% psi2)
      se2_vals[rep_i] <- sqrt(mean(iid2^2) / n)
    }

    mu_vals[rep_i] <- mu_hat
    alpha_mat[rep_i, ] <- estimates
    conv_vals[rep_i] <- converged
    grad_vals[rep_i] <- best_grad
  }, error = function(e) NULL)
}

valid <- !is.na(mu_vals) & !is.na(se1_vals)
conv <- conv_vals[valid]
mu_v <- mu_vals[valid]
se1_v <- se1_vals[valid]
se2_v <- se2_vals[valid]

cat(sprintf("mu.true = %.4f\n", mu_true))
cat(sprintf("Total valid: %d, converged: %d, not converged: %d\n", sum(valid), sum(conv), sum(!conv)))
cat(sprintf("se2 available: %d\n\n", sum(!is.na(se2_v))))

esd <- sd(mu_v)
ese1 <- mean(se1_v)
ese2 <- mean(se2_v, na.rm=TRUE)
ci1_l <- mu_v - 1.96*se1_v; ci1_u <- mu_v + 1.96*se1_v
cp1 <- mean(mu_true >= ci1_l & mu_true <= ci1_u)
ci2_l <- mu_v - 1.96*se2_v; ci2_u <- mu_v + 1.96*se2_v
cp2 <- mean(mu_true >= ci2_l & mu_true <= ci2_u, na.rm=TRUE)
cat(sprintf("=== Overall ===\n"))
cat(sprintf("Bias=%.4f, ESD=%.4f\n", mean(mu_v)-mu_true, esd))
cat(sprintf("ESE1=%.4f, ESE1/ESD=%.3f, CP1=%.3f\n", ese1, ese1/esd, cp1))
cat(sprintf("ESE2=%.4f, ESE2/ESD=%.3f, CP2=%.3f\n\n", ese2, ese2/esd, cp2))

for (label in c("Converged", "NOT converged")) {
  sel <- if (label == "Converged") conv else !conv
  n_sel <- sum(sel)
  if (n_sel == 0) next
  mu_s <- mu_vals[valid][sel]
  se1_s <- se1_vals[valid][sel]
  se2_s <- se2_vals[valid][sel]
  esd_s <- sd(mu_s)
  ese1_s <- mean(se1_s)
  ese2_s <- mean(se2_s, na.rm=TRUE)
  ci1_l <- mu_s - 1.96*se1_s; ci1_u <- mu_s + 1.96*se1_s
  cp1_s <- mean(mu_true >= ci1_l & mu_true <= ci1_u)
  ci2_l <- mu_s - 1.96*se2_s; ci2_u <- mu_s + 1.96*se2_s
  cp2_s <- mean(mu_true >= ci2_l & mu_true <= ci2_u, na.rm=TRUE)
  cat(sprintf("=== %s reps (n=%d) ===\n", label, n_sel))
  cat(sprintf("Bias=%.4f, ESD=%.4f\n", mean(mu_s)-mu_true, esd_s))
  cat(sprintf("ESE1=%.4f, ESE1/ESD=%.3f, CP1=%.3f\n", ese1_s, ese1_s/esd_s, cp1_s))
  cat(sprintf("ESE2=%.4f, ESE2/ESD=%.3f, CP2=%.3f\n\n", ese2_s, ese2_s/esd_s, cp2_s))
}

cat(sprintf("=== alpha[2] ===\n"))
cat(sprintf("mean=%.4f, sd=%.4f, range=[%.4f, %.4f]\n",
    mean(alpha_mat[valid, 2]), sd(alpha_mat[valid, 2]),
    min(alpha_mat[valid, 2]), max(alpha_mat[valid, 2])))
