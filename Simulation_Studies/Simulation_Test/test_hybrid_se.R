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
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
compute_pi_fn <- function(eta) plogis(-eta)
tol <- 1e-6
max_outer <- 100
cond_limit <- 1e8

# ============================================================
# SE1 computation (one-step GMM SE) at given alpha with W_hat
# ============================================================
compute_se1_at_alpha <- function(dat, alpha_hat, design_mat, h_x, link_type, n, W_hat) {
  r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
  pi_hat <- compute_pi_fn(as.vector(design_mat %*% alpha_hat))
  mu_hat <- mean(r_vec * y_vec / pi_hat)
  k <- length(alpha_hat); h <- ncol(h_x)

  # g_i and Gamma
  g_mat <- (r_vec / pi_hat - 1) * h_x
  cf <- r_vec * (1 - pi_hat) / pi_hat
  Gamma_hat <- crossprod(h_x * cf, design_mat) / n  # h x k

  # t_Gamma_arr: n x k x h, t_Gamma_arr[i,j,l] = cf[i]*design_mat[i,j]*h_x[i,l]
  t_Gamma_arr <- array(0, dim = c(n, k, h))
  for (j in 1:k) {
    t_Gamma_arr[, j, ] <- cf * design_mat[, j] * h_x
  }

  eta_s <- .colMeans(g_mat, n, h)
  GtW <- crossprod(Gamma_hat, W_hat)  # k x h
  GtWG <- GtW %*% Gamma_hat            # k x k

  eig <- eigen(GtWG, symmetric = TRUE, only.values = TRUE)$values
  cond_val <- max(eig) / max(min(eig), 1e-15)

  H_0n <- tryCatch(-solve(GtWG), error = function(e) {
    -solve(GtWG + 1e-8 * diag(k))
  })
  H_1n <- GtW %*% (t(g_mat) - eta_s)  # k x n

  # Gamma_2: d vec(t(Gamma)) / d alpha, (k*h) x k
  # Numerical differentiation
  Gamma_fn_local <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    cf_v <- r_vec * (1 - pi_v) / pi_v
    crossprod(h_x * cf_v, design_mat) / n
  }
  eps2 <- 1e-5
  Gamma_2 <- matrix(0, k * h, k)
  for (j in 1:k) {
    ej <- rep(0, k); ej[j] <- eps2
    Gp <- Gamma_fn_local(alpha_hat + ej)
    Gm <- Gamma_fn_local(alpha_hat - ej)
    Gamma_2[, j] <- (as.vector(t(Gp)) - as.vector(t(Gm))) / (2 * eps2)
  }

  M_n <- kronecker(eta_s %*% W_hat, diag(k)) %*% Gamma_2  # k x k

  # H_2n_2: vectorized
  W_eta <- as.vector(W_hat %*% eta_s)
  tGamma_W_eta <- as.vector(GtW %*% eta_s)
  t_Gamma_mat <- matrix(t_Gamma_arr, nrow = n * k, ncol = h)
  H_2n_2_vec <- t_Gamma_mat %*% W_eta
  H_2n_2 <- matrix(H_2n_2_vec, nrow = n, ncol = k, byrow = FALSE)
  H_2n_2 <- t(H_2n_2) - tGamma_W_eta  # k x n

  # H_2n_3
  GtW_g <- GtW %*% t(g_mat)  # k x n
  g_W_eta <- as.vector(g_mat %*% W_eta)  # n-vector
  H_2n_3 <- matrix(tGamma_W_eta, k, n) - GtW_g * matrix(g_W_eta, k, n, byrow = TRUE)

  # Q = solve(I - H_0n M_n) %*% H_0n
  I_minus_H0M <- diag(k) - H_0n %*% M_n
  Q <- tryCatch(
    solve(I_minus_H0M) %*% H_0n,
    error = function(e) solve(I_minus_H0M + 1e-8 * diag(k)) %*% H_0n
  )

  psi1 <- Q %*% (H_1n + H_2n_2 + H_2n_3)  # k x n

  # IPW SE using psi1 influence function
  if (!is.null(link_type) && link_type == "logistic_complement") {
    dot_pi <- -design_mat * (pi_hat * (1 - pi_hat))
  } else {
    dot_pi <- design_mat * (pi_hat * (1 - pi_hat))
  }
  ry_ps_inv2 <- as.vector(r_vec * y_vec * (pi_hat^(-2)))
  H_alpha <- colMeans(dot_pi * ry_ps_inv2)
  mu_iid <- as.vector(t(r_vec / pi_hat * y_vec) - t(H_alpha) %*% psi1)
  se1 <- sqrt(var(mu_iid) / n)

  list(mu = mu_hat, se1 = se1, cond = cond_val)
}

# ============================================================
# SE2 computation (iterative/CUE GMM SE) at given alpha
# ============================================================
compute_se2_at_alpha <- function(dat, alpha_hat, design_mat, h_x, link_type, n) {
  r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
  pi_hat <- compute_pi_fn(as.vector(design_mat %*% alpha_hat))
  mu_hat <- mean(r_vec * y_vec / pi_hat)
  g_mat <- (r_vec / pi_hat - 1) * h_x
  W_hat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) NULL)
  if (is.null(W_hat)) return(list(mu = mu_hat, se2 = NA, cond = NA))
  cf <- r_vec * (1 - pi_hat) / pi_hat
  Gamma_hat <- crossprod(h_x * cf, design_mat) / n
  GtWG <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
  eig <- eigen(GtWG, symmetric = TRUE, only.values = TRUE)$values
  cond_val <- max(eig) / max(min(eig), 1e-15)
  k <- length(alpha_hat); h <- ncol(h_x)
  GtW <- crossprod(Gamma_hat, W_hat)
  G_vec <- colMeans(g_mat)
  dc_deta <- cf; dc_X <- dc_deta * design_mat
  R_mat <- matrix(0, h * k, k)
  for (j in 1:k) {
    Xj_h <- design_mat[, j] * h_x
    block <- crossprod(dc_X, Xj_h) / n
    for (l in 1:h) R_mat[(l-1)*k + j, ] <- block[, l]
  }
  eta_s <- matrix(G_vec, h, 1)
  WG <- as.vector(W_hat %*% eta_s)
  GW_row <- as.vector(t(eta_s) %*% W_hat)
  S_mat <- matrix(0, h^2, k)
  for (j in 1:k) {
    dg_j <- cf * design_mat[, j] * h_x
    cross <- (crossprod(dg_j, g_mat) + crossprod(g_mat, dg_j)) / n
    S_mat[, j] <- as.vector(cross)
  }
  GW_kron_Ik <- kronecker(matrix(GW_row, 1, h), diag(k))
  GW_kron_GtW <- kronecker(matrix(GW_row, 1, h), GtW)
  H2_mat <- GtWG + GW_kron_Ik %*% R_mat - GW_kron_GtW %*% S_mat
  H2_cond <- tryCatch({
    ev <- eigen(H2_mat, symmetric = FALSE, only.values = TRUE)$values
    max(Mod(ev)) / max(min(Mod(ev)), 1e-15)
  }, error = function(e) Inf)
  if (!is.finite(H2_cond) || H2_cond > 1e15)
    return(list(mu = mu_hat, se2 = NA, cond = cond_val))
  H2_inv <- tryCatch(solve(H2_mat), error = function(e) NULL)
  if (is.null(H2_inv)) return(list(mu = mu_hat, se2 = NA, cond = cond_val))
  GtW_g <- GtW %*% t(g_mat)
  g_WG <- as.vector(g_mat %*% WG)
  cf_hxWG <- cf * as.vector(h_x %*% WG)
  GammaT_WG <- t(design_mat * cf_hxWG)
  Q_mat <- GtW_g + GammaT_WG - sweep(GtW_g, 2, g_WG, `*`)
  psi2 <- -(H2_inv %*% Q_mat)
  if (!is.null(link_type) && link_type == "logistic_complement") {
    dot_pi <- -design_mat * (pi_hat * (1 - pi_hat))
  } else {
    dot_pi <- design_mat * (pi_hat * (1 - pi_hat))
  }
  ry_ps_inv2 <- as.vector(r_vec * y_vec * (pi_hat^(-2)))
  H_alpha <- colMeans(dot_pi * ry_ps_inv2)
  mu_iid <- as.vector(t(r_vec/pi_hat*y_vec) - t(H_alpha) %*% psi2)
  se2 <- sqrt(var(mu_iid) / n)
  list(mu = mu_hat, se2 = se2, cond = cond_val)
}

# ============================================================
# Guarded iterative GMM: reject W updates that cause degeneracy
# Returns: estimates, W used, whether W was ever updated, cond
# ============================================================
run_guarded_gmm <- function(dat, design_mat, h_x, link_type) {
  r_vec <- dat[["r"]]
  n <- nrow(dat); k <- ncol(design_mat); h_dim <- ncol(h_x)

  g_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    (r_vec / pi_v - 1) * h_x
  }
  Gamma_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    crossprod(h_x * (r_vec * (1 - pi_v) / pi_v), design_mat) / n
  }
  get_cond <- function(alpha, W_hat) {
    Gamma <- Gamma_fn(alpha)
    H <- crossprod(Gamma, W_hat %*% Gamma)
    eig <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    max(eig) / max(min(eig), 1e-15)
  }
  get_grad <- function(alpha, W_hat) {
    g_mat <- g_fn(alpha)
    G <- matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
    Gamma <- Gamma_fn(alpha)
    max(abs(2 * as.vector(crossprod(Gamma, W_hat %*% G))))
  }

  nr_inner <- function(alpha, W_hat) {
    for (iter in 1:50) {
      g_i <- g_fn(alpha)
      G_i <- matrix(.colMeans(g_i, n, h_dim), h_dim, 1)
      Gamma_i <- Gamma_fn(alpha)
      gv <- 2 * as.vector(crossprod(Gamma_i, W_hat %*% G_i))
      if (max(abs(gv)) < tol) break
      GtWG <- crossprod(Gamma_i, W_hat %*% Gamma_i)
      H_mat <- GtWG + GtWG
      direction <- tryCatch(-solve(H_mat, gv),
                            error = function(e) -gv * 2.0 / max(abs(gv)))
      if (sqrt(sum(direction^2)) > 2.0)
        direction <- direction * 2.0 / sqrt(sum(direction^2))
      obj_cur <- as.numeric(crossprod(G_i, W_hat %*% G_i))
      accepted <- FALSE; step <- 1.0
      for (ls in 1:20) {
        cand <- alpha + step * direction
        g_cand <- g_fn(cand)
        G_cand <- matrix(.colMeans(g_cand, n, h_dim), h_dim, 1)
        obj_new <- as.numeric(crossprod(G_cand, W_hat %*% G_cand))
        if (is.finite(obj_new) && obj_new < obj_cur - 1e-4 * step * sum(gv * direction)) {
          Gamma_c <- Gamma_fn(cand)
          H_c <- crossprod(Gamma_c, W_hat %*% Gamma_c)
          eig_c <- eigen(H_c, symmetric = TRUE, only.values = TRUE)$values
          cc <- max(eig_c) / max(min(eig_c), 1e-15)
          if (cc < cond_limit) { accepted <- TRUE; break }
        }
        step <- step * 0.5
      }
      if (accepted) alpha <- cand else break
    }
    alpha
  }

  # Step 1: W = I
  est <- nr_inner(rep(0, k), diag(h_dim))
  W_used <- diag(h_dim)
  n_W_updates <- 0  # number of successful W updates

  best_grad <- get_grad(est, diag(h_dim))
  best_est <- est
  last_good_est <- est
  last_good_W <- diag(h_dim)

  for (t in 1:max_outer) {
    g_mat <- g_fn(est)
    W_hat <- tryCatch(W_fn(g_mat), error = function(e) diag(h_dim))

    cur_cond <- get_cond(est, W_hat)
    if (cur_cond >= cond_limit) {
      est <- last_good_est
      break
    }

    Gamma <- Gamma_fn(est)
    G_vec <- matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
    grad <- max(abs(2 * as.vector(crossprod(Gamma, W_hat %*% G_vec))))

    if (grad < tol) {
      best_est <- est; best_grad <- grad
      last_good_est <- est; last_good_W <- W_hat
      n_W_updates <- n_W_updates + 1
      W_used <- W_hat
      break
    }

    last_good_est <- est; last_good_W <- W_hat
    n_W_updates <- n_W_updates + 1
    W_used <- W_hat

    new_est <- nr_inner(est, W_hat)
    g_new <- g_fn(new_est)
    W_new <- tryCatch(W_fn(g_new), error = function(e) diag(h_dim))
    new_cond <- get_cond(new_est, W_new)
    new_grad <- get_grad(new_est, W_new)

    if (new_cond < cond_limit) {
      est <- new_est
      if (new_grad < best_grad) {
        best_grad <- new_grad; best_est <- new_est
      }
    } else {
      est <- last_good_est
      break
    }
  }

  list(
    estimates = best_est,
    grad_norm = best_grad,
    W_used = last_good_W,
    n_W_updates = n_W_updates,
    is_onestep = (n_W_updates == 0)  # TRUE if never successfully updated W
  )
}

# ========== Run for M2 and M3 separately ==========
for (m_idx in c(2, 3)) {
  cat(sprintf("\n=== Model %d (guarded W + hybrid SE) ===\n", m_idx))
  mu_vals <- se_vals <- rep(NA_real_, n_reps)
  se_type <- rep(NA_character_, n_reps)
  n_onestep <- 0

  for (rep_i in 1:n_reps) {
    dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
    single_ps <- list(
      formula.list = list(ps_spec[["formula.list"]][[m_idx]]),
      h_alpha.list = list(ps_spec[["h_alpha.list"]][[m_idx]]),
      inv_link = ps_spec[["inv_link"]],
      outcome = ps_spec[["outcome"]],
      alpha_init.list = list(NULL),
      optimizer = "L-BFGS-B"
    )
    tryCatch({
      ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
      design_mat <- ebmr$ps_fit.list[[1]]$design_matrix
      h_x <- ebmr$ps_fit.list[[1]]$h_x
      link_type <- ebmr$ps_fit.list[[1]]$link_type

      res <- run_guarded_gmm(dat, design_mat, h_x, link_type)

      if (res$is_onestep) {
        # W update was rejected → one-step GMM → use SE1
        n_onestep <- n_onestep + 1
        se_res <- compute_se1_at_alpha(dat, res$estimates, design_mat, h_x, link_type, n_val, diag(ncol(h_x)))
        mu_vals[rep_i] <- se_res$mu
        se_vals[rep_i] <- se_res$se1
        se_type[rep_i] <- "SE1"
      } else {
        # Iterative GMM converged → use SE2
        se_res <- compute_se2_at_alpha(dat, res$estimates, design_mat, h_x, link_type, n_val)
        mu_vals[rep_i] <- se_res$mu
        se_vals[rep_i] <- se_res$se2
        se_type[rep_i] <- "SE2"
      }
    }, error = function(e) {
      cat(sprintf("  Rep %d ERROR: %s\n", rep_i, e$message))
    })
  }

  valid <- !is.na(mu_vals) & !is.na(se_vals)
  cat(sprintf("  Valid: %d/%d, One-step (SE1): %d, Iterative (SE2): %d\n",
      sum(valid), n_reps, sum(se_type == "SE1", na.rm = TRUE), sum(se_type == "SE2", na.rm = TRUE)))
  cat(sprintf("  Bias  = %.4f\n", mean(mu_vals[valid]) - mu_true))
  cat(sprintf("  ESD   = %.4f\n", sd(mu_vals[valid])))
  cat(sprintf("  ESE   = %.4f\n", mean(se_vals[valid])))
  cat(sprintf("  ESE/ESD = %.3f\n", mean(se_vals[valid]) / sd(mu_vals[valid])))

  # Breakdown by SE type
  se1_idx <- valid & se_type == "SE1"
  se2_idx <- valid & se_type == "SE2"
  if (sum(se1_idx) > 1) {
    cat(sprintf("  [SE1 reps] n=%d, Bias=%.4f, ESD=%.4f, ESE=%.4f, ESE/ESD=%.3f\n",
        sum(se1_idx), mean(mu_vals[se1_idx]) - mu_true, sd(mu_vals[se1_idx]),
        mean(se_vals[se1_idx]), mean(se_vals[se1_idx]) / sd(mu_vals[se1_idx])))
  }
  if (sum(se2_idx) > 1) {
    cat(sprintf("  [SE2 reps] n=%d, Bias=%.4f, ESD=%.4f, ESE=%.4f, ESE/ESD=%.3f\n",
        sum(se2_idx), mean(mu_vals[se2_idx]) - mu_true, sd(mu_vals[se2_idx]),
        mean(se_vals[se2_idx]), mean(se_vals[se2_idx]) / sd(mu_vals[se2_idx])))
  }
}
