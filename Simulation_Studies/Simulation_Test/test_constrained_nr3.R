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

run_method <- function(method_name, inner_fn) {
  mu_vals <- rep(NA_real_, n_reps)
  se1_vals <- rep(NA_real_, n_reps)
  se2_vals <- rep(NA_real_, n_reps)
  alpha2_vals <- rep(NA_real_, n_reps)
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

      get_condH <- function(alpha, W_hat) {
        Gamma_hat <- Gamma_fn(alpha)
        H <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
        eig <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
        max(eig) / max(min(eig), 1e-15)
      }

      # Step 1: W = I
      estimates <- inner_fn(rep(0, p), diag(h_dim), p, h_dim,
                            Phi_alpha, Gamma_fn, W_func, get_condH)

      # Step 2+
      best_grad <- Inf; best_est <- estimates; no_improve <- 0L; converged <- FALSE
      for (t in 1:200) {
        g_mat <- Phi_alpha(estimates)
        G_vec <- colMeans(g_mat)
        W_hat <- tryCatch(W_func(g_mat), error = function(e) diag(h_dim))
        Gamma_hat <- Gamma_fn(estimates)
        conv_grad <- 2 * as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))
        grad_norm <- max(abs(conv_grad))
        if (grad_norm < 1e-6) { converged <- TRUE; break }
        if (grad_norm < best_grad) { best_grad <- grad_norm; best_est <- estimates; no_improve <- 0L
        } else { no_improve <- no_improve + 1L; if (no_improve >= 20L) { estimates <- best_est; break } }
        estimates <- inner_fn(estimates, W_hat, p, h_dim,
                              Phi_alpha, Gamma_fn, W_func, get_condH)
      }

      pi_hat <- model_fn(estimates)
      mu_hat <- mean(r_vec * y_vec / pi_hat)

      g_mat <- Phi_alpha(estimates)
      W_hat <- tryCatch(W_func(g_mat), error = function(e) diag(h_dim))
      Gamma_hat <- Gamma_fn(estimates)
      G_vec <- colMeans(g_mat)
      H_mat <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
      H_inv <- tryCatch(solve(H_mat), error = function(e) NULL)

      if (!is.null(H_inv)) {
        GtW <- crossprod(Gamma_hat, W_hat)
        psi1 <- H_inv %*% GtW %*% t(g_mat)
        dmu_dalpha <- colMeans(-r_vec * y_vec * (1 - pi_hat) / pi_hat * design_mat)
        iid1 <- r_vec * y_vec / pi_hat - mu_hat - as.vector(t(dmu_dalpha) %*% psi1)
        se1_vals[rep_i] <- sqrt(mean(iid1^2) / n)

        WG <- W_hat %*% G_vec
        gi_dot_WG <- as.vector(g_mat %*% WG)
        Omega <- crossprod(g_mat) / n
        OmegaWG <- as.vector(Omega %*% WG)
        GtW_g <- GtW %*% t(g_mat)
        GtW_OmegaWG <- as.vector(GtW %*% OmegaWG)
        correction <- sweep(GtW_g * rep(gi_dot_WG, each = p), 1, GtW_OmegaWG, "-")
        psi2 <- psi1 - H_inv %*% correction
        iid2 <- r_vec * y_vec / pi_hat - mu_hat - as.vector(t(dmu_dalpha) %*% psi2)
        se2_vals[rep_i] <- sqrt(mean(iid2^2) / n)
      }

      mu_vals[rep_i] <- mu_hat
      alpha2_vals[rep_i] <- estimates[2]
      conv_vals[rep_i] <- converged
      grad_vals[rep_i] <- best_grad
    }, error = function(e) NULL)
  }

  valid <- !is.na(mu_vals) & !is.na(se1_vals)
  conv <- conv_vals[valid]
  esd <- sd(mu_vals[valid])
  ese1 <- mean(se1_vals[valid])
  ese2 <- mean(se2_vals[valid], na.rm=TRUE)

  cat(sprintf("\n=== %s ===\n", method_name))
  cat(sprintf("Valid: %d, Converged: %d\n", sum(valid), sum(conv)))
  cat(sprintf("ESD=%.4f, ESE1/ESD=%.3f, ESE2/ESD=%.3f\n", esd, ese1/esd, ese2/esd))
  cat(sprintf("alpha[2]: mean=%.3f, sd=%.3f, range=[%.3f, %.3f]\n",
      mean(alpha2_vals[valid]), sd(alpha2_vals[valid]),
      min(alpha2_vals[valid]), max(alpha2_vals[valid])))
}

# ============================================================
# Method A: cond check + trust region (step norm cap = 2.0)
# ============================================================
inner_A <- function(start, W_hat, p, h_dim, Phi_alpha, Gamma_fn, W_func, get_condH) {
  alpha <- start; cond_thresh <- 1e7
  for (k in 1:500) {
    g_mat <- Phi_alpha(alpha); G_vec <- colMeans(g_mat)
    Gamma_hat <- Gamma_fn(alpha)
    g_vec <- 2 * as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))
    if (max(abs(g_vec)) < 1e-8) break
    H_gn <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
    direction <- tryCatch(solve(H_gn, -g_vec), error = function(e) -g_vec)
    sn <- sqrt(sum(direction^2))
    if (sn > 2.0) direction <- direction * (2.0 / sn)  # trust region = 2.0
    obj_cur <- as.numeric(t(G_vec) %*% W_hat %*% G_vec)
    step <- 1.0; accepted <- FALSE
    for (ls in 1:30) {
      cand <- alpha + step * direction
      G_new <- colMeans(Phi_alpha(cand))
      obj_new <- as.numeric(t(G_new) %*% W_hat %*% G_new)
      if (is.finite(obj_new) && obj_new < obj_cur - 1e-4 * step * sum(g_vec * direction)) {
        if (get_condH(cand, W_hat) < cond_thresh) { accepted <- TRUE; break }
      }
      step <- step * 0.5
    }
    if (accepted) { alpha <- cand } else { break }
  }
  alpha
}

# ============================================================
# Method B: adaptive threshold (reject if cond > 100x starting cond)
# ============================================================
inner_B <- function(start, W_hat, p, h_dim, Phi_alpha, Gamma_fn, W_func, get_condH) {
  alpha <- start
  cond_start <- get_condH(start, W_hat)
  cond_thresh <- max(cond_start * 100, 1e5)  # at least 1e5
  for (k in 1:500) {
    g_mat <- Phi_alpha(alpha); G_vec <- colMeans(g_mat)
    Gamma_hat <- Gamma_fn(alpha)
    g_vec <- 2 * as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))
    if (max(abs(g_vec)) < 1e-8) break
    H_gn <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
    direction <- tryCatch(solve(H_gn, -g_vec), error = function(e) -g_vec)
    obj_cur <- as.numeric(t(G_vec) %*% W_hat %*% G_vec)
    step <- 1.0; accepted <- FALSE
    for (ls in 1:30) {
      cand <- alpha + step * direction
      G_new <- colMeans(Phi_alpha(cand))
      obj_new <- as.numeric(t(G_new) %*% W_hat %*% G_new)
      if (is.finite(obj_new) && obj_new < obj_cur - 1e-4 * step * sum(g_vec * direction)) {
        if (get_condH(cand, W_hat) < cond_thresh) { accepted <- TRUE; break }
      }
      step <- step * 0.5
    }
    if (accepted) { alpha <- cand } else { break }
  }
  alpha
}

# ============================================================
# Method C: cond check + trust region 1.0
# ============================================================
inner_C <- function(start, W_hat, p, h_dim, Phi_alpha, Gamma_fn, W_func, get_condH) {
  alpha <- start; cond_thresh <- 1e7
  for (k in 1:500) {
    g_mat <- Phi_alpha(alpha); G_vec <- colMeans(g_mat)
    Gamma_hat <- Gamma_fn(alpha)
    g_vec <- 2 * as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))
    if (max(abs(g_vec)) < 1e-8) break
    H_gn <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
    direction <- tryCatch(solve(H_gn, -g_vec), error = function(e) -g_vec)
    sn <- sqrt(sum(direction^2))
    if (sn > 1.0) direction <- direction * (1.0 / sn)
    obj_cur <- as.numeric(t(G_vec) %*% W_hat %*% G_vec)
    step <- 1.0; accepted <- FALSE
    for (ls in 1:30) {
      cand <- alpha + step * direction
      G_new <- colMeans(Phi_alpha(cand))
      obj_new <- as.numeric(t(G_new) %*% W_hat %*% G_new)
      if (is.finite(obj_new) && obj_new < obj_cur - 1e-4 * step * sum(g_vec * direction)) {
        if (get_condH(cand, W_hat) < cond_thresh) { accepted <- TRUE; break }
      }
      step <- step * 0.5
    }
    if (accepted) { alpha <- cand } else { break }
  }
  alpha
}

# ============================================================
# Method D: min eigenvalue check (reject if min_eig(H) < threshold)
# Instead of cond number, directly check if H is near-singular
# ============================================================
inner_D <- function(start, W_hat, p, h_dim, Phi_alpha, Gamma_fn, W_func, get_condH) {
  alpha <- start; min_eig_thresh <- 1e-8
  for (k in 1:500) {
    g_mat <- Phi_alpha(alpha); G_vec <- colMeans(g_mat)
    Gamma_hat <- Gamma_fn(alpha)
    g_vec <- 2 * as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))
    if (max(abs(g_vec)) < 1e-8) break
    H_gn <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
    direction <- tryCatch(solve(H_gn, -g_vec), error = function(e) -g_vec)
    obj_cur <- as.numeric(t(G_vec) %*% W_hat %*% G_vec)
    step <- 1.0; accepted <- FALSE
    for (ls in 1:30) {
      cand <- alpha + step * direction
      G_new <- colMeans(Phi_alpha(cand))
      obj_new <- as.numeric(t(G_new) %*% W_hat %*% G_new)
      if (is.finite(obj_new) && obj_new < obj_cur - 1e-4 * step * sum(g_vec * direction)) {
        Gamma_cand <- Gamma_fn(cand)
        H_cand <- crossprod(Gamma_cand, W_hat %*% Gamma_cand)
        eig_cand <- eigen(H_cand, symmetric = TRUE, only.values = TRUE)$values
        if (min(eig_cand) > min_eig_thresh) { accepted <- TRUE; break }
      }
      step <- step * 0.5
    }
    if (accepted) { alpha <- cand } else { break }
  }
  alpha
}

# ============================================================
# Method E: adaptive + trust region 1.0
# ============================================================
inner_E <- function(start, W_hat, p, h_dim, Phi_alpha, Gamma_fn, W_func, get_condH) {
  alpha <- start
  cond_start <- get_condH(start, W_hat)
  cond_thresh <- max(cond_start * 10, 1e4)
  for (k in 1:500) {
    g_mat <- Phi_alpha(alpha); G_vec <- colMeans(g_mat)
    Gamma_hat <- Gamma_fn(alpha)
    g_vec <- 2 * as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))
    if (max(abs(g_vec)) < 1e-8) break
    H_gn <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
    direction <- tryCatch(solve(H_gn, -g_vec), error = function(e) -g_vec)
    sn <- sqrt(sum(direction^2))
    if (sn > 1.0) direction <- direction * (1.0 / sn)
    obj_cur <- as.numeric(t(G_vec) %*% W_hat %*% G_vec)
    step <- 1.0; accepted <- FALSE
    for (ls in 1:30) {
      cand <- alpha + step * direction
      G_new <- colMeans(Phi_alpha(cand))
      obj_new <- as.numeric(t(G_new) %*% W_hat %*% G_new)
      if (is.finite(obj_new) && obj_new < obj_cur - 1e-4 * step * sum(g_vec * direction)) {
        if (get_condH(cand, W_hat) < cond_thresh) { accepted <- TRUE; break }
      }
      step <- step * 0.5
    }
    if (accepted) { alpha <- cand } else { break }
  }
  alpha
}

run_method("A: cond<1e7 + trust=2.0", inner_A)
run_method("B: adaptive cond (100x start)", inner_B)
run_method("C: cond<1e7 + trust=1.0", inner_C)
run_method("D: min_eig(H) > 1e-8", inner_D)
run_method("E: adaptive cond (10x) + trust=1.0", inner_E)
