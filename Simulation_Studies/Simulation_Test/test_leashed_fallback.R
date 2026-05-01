setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

n_val <- 2000
n_reps <- 1000
ps_spec <- get_ps_spec("9-alt1")
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
compute_pi_fn <- function(eta) plogis(-eta)
tol <- 1e-6
cond_limit <- 1e8
dist_max <- 6

# Leashed inner + L-BFGS-B fallback when leash doesn't converge
run_leashed_fallback_gmm <- function(dat, design_mat, h_x, link_type) {
  r_vec <- dat[["r"]]
  n <- nrow(dat); k <- ncol(design_mat); h_dim <- ncol(h_x)
  init <- rep(0, k)

  g_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    (r_vec / pi_v - 1) * h_x
  }
  Gamma_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    crossprod(h_x * (r_vec * (1 - pi_v) / pi_v), design_mat) / n
  }
  check_cond_W <- function(alpha, W_hat) {
    Gamma <- Gamma_fn(alpha)
    H <- crossprod(Gamma, W_hat %*% Gamma)
    eig <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    max(eig) / max(min(eig), 1e-15)
  }
  forward_cond <- function(alpha) {
    g_c <- g_fn(alpha)
    W_c <- tryCatch(solve(crossprod(g_c) / n), error = function(e) diag(h_dim))
    check_cond_W(alpha, W_c)
  }
  compute_grad_with_W <- function(alpha) {
    g_mat <- g_fn(alpha)
    G_vec <- matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
    W_hat <- tryCatch(W_fn(g_mat), error = function(e) diag(h_dim))
    Gamma <- Gamma_fn(alpha)
    max(abs(2 * as.vector(crossprod(Gamma, W_hat %*% G_vec))))
  }

  # ===== Phase 1: Leashed inner NR =====
  nr_inner <- function(start, W_hat) {
    alpha <- start
    for (iter in 1:500) {
      g_mat <- g_fn(alpha)
      G_vec <- matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
      Gamma <- Gamma_fn(alpha)
      gv <- 2 * as.vector(crossprod(Gamma, W_hat %*% G_vec))
      if (max(abs(gv)) < 1e-8) break
      H_gn <- crossprod(Gamma, W_hat %*% Gamma)
      direction <- tryCatch(solve(H_gn, -gv), error = function(e) -gv / max(abs(gv)))
      sn <- sqrt(sum(direction^2))
      if (sn > 2.0) direction <- direction * (2.0 / sn)
      obj_cur <- as.numeric(crossprod(G_vec, W_hat %*% G_vec))
      accepted <- FALSE; step <- 1.0
      for (ls in 1:30) {
        cand <- alpha + step * direction
        if (sqrt(sum((cand - init)^2)) > dist_max) { step <- step * 0.5; next }
        g_cand <- g_fn(cand)
        G_cand <- matrix(.colMeans(g_cand, n, h_dim), h_dim, 1)
        obj_new <- as.numeric(crossprod(G_cand, W_hat %*% G_cand))
        if (is.finite(obj_new) && obj_new < obj_cur - 1e-4 * step * sum(gv * direction)) {
          cc <- check_cond_W(cand, W_hat)
          if (cc < cond_limit) { accepted <- TRUE; break }
        }
        step <- step * 0.5
      }
      if (accepted) alpha <- cand else break
    }
    alpha
  }

  est <- init; best_grad <- Inf; best_est <- est
  for (t in 1:100) {
    g_mat <- g_fn(est)
    W_hat <- tryCatch(W_fn(g_mat), error = function(e) diag(h_dim))
    Gamma <- Gamma_fn(est)
    G_vec <- matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
    grad <- max(abs(2 * as.vector(crossprod(Gamma, W_hat %*% G_vec))))
    if (grad < best_grad) { best_grad <- grad; best_est <- est }
    if (grad < tol) break
    est <- nr_inner(est, W_hat)
  }

  leash_grad <- compute_grad_with_W(best_est)
  used_fallback <- FALSE

  # ===== Phase 2: L-BFGS-B fallback if leash didn't converge =====
  if (leash_grad >= tol) {
    run_lbfgsb <- function(start_point) {
      st <- new.env(parent = emptyenv())
      st$W.hat <- NULL; st$cache_param <- NULL
      st$cache_g <- NULL; st$cache_G <- NULL
      cache_fn <- function(param) {
        if (!is.null(st$cache_param) && identical(param, st$cache_param)) return()
        st$cache_param <- param
        st$cache_g <- g_fn(param)
        st$cache_G <- matrix(.colMeans(st$cache_g, n, h_dim), h_dim, 1)
      }
      obj_fn <- function(param) {
        cache_fn(param)
        v <- as.numeric(crossprod(st$cache_G, st$W.hat %*% st$cache_G))
        if (!is.finite(v)) 1e8 else v
      }
      grad_fn <- function(param) {
        cache_fn(param)
        2 * as.vector(crossprod(Gamma_fn(param), st$W.hat %*% st$cache_G))
      }
      est_fb <- start_point
      # Step 1: W = I
      st$W.hat <- diag(h_dim)
      opt0 <- optim(est_fb, obj_fn, gr = grad_fn, method = "L-BFGS-B",
                     control = list(maxit = 1000))
      est_fb <- opt0$par
      # Step 2: Iterative W updates
      for (t_fb in 1:100) {
        cache_fn(est_fb)
        st$W.hat <- tryCatch(W_fn(st$cache_g), error = function(e) diag(h_dim))
        fg <- max(abs(grad_fn(est_fb)))
        if (fg < tol) break
        st$cache_param <- NULL
        opt_t <- optim(est_fb, obj_fn, gr = grad_fn, method = "L-BFGS-B",
                       control = list(maxit = 1000))
        est_fb <- opt_t$par
      }
      g_f <- g_fn(est_fb)
      G_f <- matrix(.colMeans(g_f, n, h_dim), h_dim, 1)
      W_f <- tryCatch(W_fn(g_f), error = function(e) diag(h_dim))
      list(estimates = est_fb,
           grad_norm = max(abs(2 * as.vector(crossprod(Gamma_fn(est_fb), W_f %*% G_f)))))
    }

    # Try from two starting points
    fb1 <- run_lbfgsb(best_est)
    fb2 <- run_lbfgsb(init)
    fb <- if (fb1$grad_norm <= fb2$grad_norm) fb1 else fb2

    # Only accept fallback if it converges better AND is not degenerate
    fb_cond <- forward_cond(fb$estimates)
    if (fb$grad_norm < leash_grad && fb_cond < cond_limit) {
      best_est <- fb$estimates
      used_fallback <- TRUE
    }
  }

  final_cond <- forward_cond(best_est)
  final_grad <- compute_grad_with_W(best_est)
  pi_final <- compute_pi_fn(as.vector(design_mat %*% best_est))
  list(estimates = best_est, grad_norm = final_grad, cond = final_cond,
       pi_hat = pi_final, used_fallback = used_fallback)
}

# SE computation
compute_se <- function(dat, alpha_hat, design_mat, h_x, link_type, n) {
  r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
  pi_hat <- compute_pi_fn(as.vector(design_mat %*% alpha_hat))
  mu_hat <- mean(r_vec * y_vec / pi_hat)
  g_mat <- (r_vec / pi_hat - 1) * h_x
  W_hat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) NULL)
  if (is.null(W_hat)) return(list(mu = mu_hat, se = NA))
  cf <- r_vec * (1 - pi_hat) / pi_hat
  Gamma_hat <- crossprod(h_x * cf, design_mat) / n
  GtWG <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
  k <- length(alpha_hat); h <- ncol(h_x)
  GtW <- crossprod(Gamma_hat, W_hat)
  G_vec <- colMeans(g_mat)
  dc_X <- cf * design_mat
  R_mat <- matrix(0, h * k, k)
  for (j in 1:k) {
    block <- crossprod(dc_X, design_mat[, j] * h_x) / n
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
  H_mat <- GtWG + GW_kron_Ik %*% R_mat - GW_kron_GtW %*% S_mat
  H_cond <- tryCatch({
    ev <- eigen(H_mat, symmetric = FALSE, only.values = TRUE)$values
    max(Mod(ev)) / max(min(Mod(ev)), 1e-15)
  }, error = function(e) Inf)
  if (!is.finite(H_cond) || H_cond > 1e15) return(list(mu = mu_hat, se = NA))
  H_inv <- tryCatch(solve(H_mat), error = function(e) NULL)
  if (is.null(H_inv)) return(list(mu = mu_hat, se = NA))
  GtW_g <- GtW %*% t(g_mat)
  g_WG <- as.vector(g_mat %*% WG)
  cf_hxWG <- cf * as.vector(h_x %*% WG)
  GammaT_WG <- t(design_mat * cf_hxWG)
  Q_mat <- GtW_g + GammaT_WG - sweep(GtW_g, 2, g_WG, `*`)
  psi <- -(H_inv %*% Q_mat)
  if (!is.null(link_type) && link_type == "logistic_complement") {
    dot_pi <- -design_mat * (pi_hat * (1 - pi_hat))
  } else {
    dot_pi <- design_mat * (pi_hat * (1 - pi_hat))
  }
  ry_ps_inv2 <- as.vector(r_vec * y_vec * (pi_hat^(-2)))
  H_alpha <- colMeans(dot_pi * ry_ps_inv2)
  mu_iid <- as.vector(t(r_vec/pi_hat*y_vec) - t(H_alpha) %*% psi)
  se <- sqrt(var(mu_iid) / n)
  list(mu = mu_hat, se = se)
}

# ===== Run all three settings =====
settings <- list(
  list(setting = "setting4", m_idx = 3, label = "S4 M3"),
  list(setting = "setting4", m_idx = 2, label = "S4 M2"),
  list(setting = "setting3", m_idx = 2, label = "S3 M2")
)

for (s in settings) {
  cat(sprintf("\n========== %s, dist_max=%d, 1000 reps ==========\n", s$label, dist_max))
  data_file <- misspecified_model_all_data_file.list[[s$setting]][["miss50"]][[1]]
  all_data <- readRDS(data_file)
  mu_true <- get_mu_true(s$setting)

  # Pre-build design matrices
  cat("  Pre-building...\n")
  design_info <- vector("list", n_reps)
  for (rep_i in 1:n_reps) {
    dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
    single_ps <- list(
      formula.list = list(ps_spec[["formula.list"]][[s$m_idx]]),
      h_alpha.list = list(ps_spec[["h_alpha.list"]][[s$m_idx]]),
      inv_link = ps_spec[["inv_link"]],
      outcome = ps_spec[["outcome"]],
      alpha_init.list = list(NULL),
      optimizer = "L-BFGS-B"
    )
    tryCatch({
      ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
      design_info[[rep_i]] <- list(
        design_mat = ebmr$ps_fit.list[[1]]$design_matrix,
        h_x = ebmr$ps_fit.list[[1]]$h_x,
        link_type = ebmr$ps_fit.list[[1]]$link_type
      )
    }, error = function(e) NULL)
  }

  mu_vals <- se_vals <- cond_vals <- grad_vals <- rep(NA_real_, n_reps)
  n_fallback <- 0
  for (rep_i in 1:n_reps) {
    if (rep_i %% 200 == 0) cat(sprintf("  rep %d...\n", rep_i))
    if (is.null(design_info[[rep_i]])) next
    dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
    di <- design_info[[rep_i]]
    tryCatch({
      res <- run_leashed_fallback_gmm(dat, di$design_mat, di$h_x, di$link_type)
      cond_vals[rep_i] <- res$cond
      grad_vals[rep_i] <- res$grad_norm
      if (res$used_fallback) n_fallback <- n_fallback + 1
      se_res <- compute_se(dat, res$estimates, di$design_mat, di$h_x, di$link_type, n_val)
      mu_vals[rep_i] <- se_res$mu
      se_vals[rep_i] <- se_res$se
    }, error = function(e) {
      cat(sprintf("  Rep %d ERROR: %s\n", rep_i, e$message))
    })
  }

  valid <- !is.na(mu_vals) & !is.na(se_vals)
  n_degen <- sum(cond_vals[valid] >= cond_limit, na.rm = TRUE)
  n_conv <- sum(grad_vals[valid] < tol, na.rm = TRUE)
  cat(sprintf("  Valid: %d/%d, Degen: %d, Conv: %d, Fallback: %d\n",
      sum(valid), n_reps, n_degen, n_conv, n_fallback))
  cat(sprintf("  Bias     = %.4f\n", mean(mu_vals[valid]) - mu_true))
  cat(sprintf("  ESD      = %.4f\n", sd(mu_vals[valid])))
  cat(sprintf("  ESE      = %.4f\n", mean(se_vals[valid])))
  cat(sprintf("  ESE/ESD  = %.3f\n", mean(se_vals[valid]) / sd(mu_vals[valid])))
  ci_lo <- mu_vals[valid] - 1.96 * se_vals[valid]
  ci_hi <- mu_vals[valid] + 1.96 * se_vals[valid]
  cat(sprintf("  CP       = %.3f\n", mean(mu_true >= ci_lo & mu_true <= ci_hi)))
}
