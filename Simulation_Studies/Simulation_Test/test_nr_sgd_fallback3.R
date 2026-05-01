setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

n_val <- 2000
ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
mu_true <- get_mu_true("setting4")
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
compute_pi_fn <- function(eta) plogis(-eta)
tol <- 1e-6
cond_limit <- 1e6  # Tighter threshold

# Stubborn reps that had cond ~1e7 in fallback2
stubborn_reps <- c(35, 38, 143, 155, 156, 190)
# Also test a few well-behaved fallback reps for comparison
good_fallback_reps <- c(34, 68, 109)
all_test_reps <- c(stubborn_reps, good_fallback_reps)
m_idx <- 3

run_sgd_nr_fallback_v3 <- function(dat, design_mat, h_x, link_type, n_restarts = 10) {
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
  obj_fn <- function(alpha, W_hat) {
    g_mat <- g_fn(alpha)
    G <- matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
    as.numeric(crossprod(G, W_hat %*% G))
  }
  grad_fn <- function(alpha, W_hat) {
    g_mat <- g_fn(alpha)
    G <- matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
    Gamma <- Gamma_fn(alpha)
    2 * as.vector(crossprod(Gamma, W_hat %*% G))
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
  full_grad_norm <- function(alpha) {
    g_mat <- g_fn(alpha)
    G <- matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
    W <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(h_dim))
    Gamma <- Gamma_fn(alpha)
    max(abs(2 * as.vector(crossprod(Gamma, W %*% G))))
  }

  sgd_grad <- function(alpha, W_hat, idx) {
    nb <- length(idx)
    dm_b <- design_mat[idx, , drop = FALSE]
    hx_b <- h_x[idx, , drop = FALSE]
    r_b <- r_vec[idx]
    pi_v <- compute_pi_fn(as.vector(dm_b %*% alpha))
    g_b <- (r_b / pi_v - 1) * hx_b
    G_b <- matrix(colMeans(g_b), h_dim, 1)
    cf_b <- r_b * (1 - pi_v) / pi_v
    Gamma_b <- crossprod(hx_b * cf_b, dm_b) / nb
    2 * as.vector(crossprod(Gamma_b, W_hat %*% G_b))
  }

  sgd_inner <- function(start, W_hat, batch_size = 256, max_iter = 3000, lr0 = 0.05) {
    alpha <- start; best_alpha <- alpha; best_obj <- Inf
    for (iter in 1:max_iter) {
      idx <- sample.int(n, min(batch_size, n))
      sg <- sgd_grad(alpha, W_hat, idx)
      lr <- lr0 / sqrt(iter)
      step_vec <- -lr * sg
      sn <- sqrt(sum(step_vec^2))
      if (sn > 1.0) step_vec <- step_vec * (1.0 / sn)
      alpha <- alpha + step_vec
      if (iter %% 100 == 0) {
        obj <- obj_fn(alpha, W_hat)
        if (is.finite(obj) && obj < best_obj) { best_obj <- obj; best_alpha <- alpha }
      }
    }
    best_alpha
  }

  nr_inner <- function(start, W_hat, trust_r = 2.0) {
    alpha <- start
    for (iter in 1:500) {
      g_mat <- g_fn(alpha)
      G_vec <- matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
      Gamma <- Gamma_fn(alpha)
      gv <- 2 * as.vector(crossprod(Gamma, W_hat %*% G_vec))
      if (max(abs(gv)) < 1e-8) break
      H_gn <- crossprod(Gamma, W_hat %*% Gamma)
      direction <- tryCatch(solve(H_gn, -gv), error = function(e) -gv)
      sn <- sqrt(sum(direction^2))
      if (sn > trust_r) direction <- direction * (trust_r / sn)
      obj_cur <- as.numeric(crossprod(G_vec, W_hat %*% G_vec))
      accepted <- FALSE; step <- 1.0
      for (ls in 1:30) {
        cand <- alpha + step * direction
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

  best_overall_est <- NULL
  best_overall_grad <- Inf

  for (restart in 1:n_restarts) {
    if (restart == 1) {
      init <- rep(0, k)
    } else {
      init <- rnorm(k) * runif(1, 0.1, 1.0)  # Wider diversity
    }

    # SGD phase
    est <- sgd_inner(init, diag(h_dim), max_iter = 3000, lr0 = 0.05)
    for (t in 1:5) {
      g_mat <- g_fn(est)
      W_hat <- tryCatch(W_fn(g_mat), error = function(e) diag(h_dim))
      est <- sgd_inner(est, W_hat, max_iter = 1500, lr0 = 0.03)
    }

    # NR phase
    prev_est <- est; prev_W <- diag(h_dim)
    best_grad <- Inf; best_est <- est

    for (t in 1:50) {
      g_mat <- g_fn(est)
      W_hat <- tryCatch(W_fn(g_mat), error = function(e) diag(h_dim))
      cc <- check_cond_W(est, W_hat)

      if (cc >= cond_limit) {
        retry_ok <- FALSE
        for (tr in c(1.0, 0.5, 0.2, 0.1)) {
          retry <- nr_inner(prev_est, prev_W, trust_r = tr)
          if (forward_cond(retry) < cond_limit) {
            est <- retry; retry_ok <- TRUE; break
          }
        }
        if (!retry_ok) { est <- best_est; break }
        next
      }

      Gamma <- Gamma_fn(est)
      G_vec <- matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
      grad <- max(abs(2 * as.vector(crossprod(Gamma, W_hat %*% G_vec))))
      if (grad < best_grad) { best_grad <- grad; best_est <- est }
      if (grad < tol) break

      prev_est <- est; prev_W <- W_hat
      est <- nr_inner(est, W_hat)
    }

    fc <- forward_cond(best_est)
    fg <- full_grad_norm(best_est)

    cat(sprintf("    restart %2d: grad=%.2e, cond=%.2e, alpha[2]=%.3f\n",
        restart, fg, fc, best_est[2]))

    if (fc < cond_limit && fg < best_overall_grad) {
      best_overall_grad <- fg
      best_overall_est <- best_est
    }
  }

  if (is.null(best_overall_est)) {
    best_overall_est <- rep(0, k)
    best_overall_grad <- Inf
  }

  final_cond <- forward_cond(best_overall_est)
  list(estimates = best_overall_est, grad_norm = best_overall_grad, cond = final_cond)
}

# SE2
compute_se2 <- function(dat, alpha_hat, design_mat, h_x, link_type, n) {
  r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
  pi_hat <- compute_pi_fn(as.vector(design_mat %*% alpha_hat))
  mu_hat <- mean(r_vec * y_vec / pi_hat)
  g_mat <- (r_vec / pi_hat - 1) * h_x
  W_hat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) NULL)
  if (is.null(W_hat)) return(list(mu = mu_hat, se2 = NA))
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
  H2_mat <- GtWG + GW_kron_Ik %*% R_mat - GW_kron_GtW %*% S_mat
  H2_cond <- tryCatch({
    ev <- eigen(H2_mat, symmetric = FALSE, only.values = TRUE)$values
    max(Mod(ev)) / max(min(Mod(ev)), 1e-15)
  }, error = function(e) Inf)
  if (!is.finite(H2_cond) || H2_cond > 1e15) return(list(mu = mu_hat, se2 = NA))
  H2_inv <- tryCatch(solve(H2_mat), error = function(e) NULL)
  if (is.null(H2_inv)) return(list(mu = mu_hat, se2 = NA))
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
  list(mu = mu_hat, se2 = se2)
}

# Test on stubborn + good reps
cat("=== SGD->NR v3: cond_limit=1e6, 10 restarts, 3000 SGD iters ===\n\n")
for (rep_i in all_test_reps) {
  cat(sprintf("--- Rep %d ---\n", rep_i))
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

    set.seed(rep_i * 200)
    res <- run_sgd_nr_fallback_v3(dat, design_mat, h_x, link_type, n_restarts = 10)

    se_res <- compute_se2(dat, res$estimates, design_mat, h_x, link_type, n_val)
    cat(sprintf("  BEST: mu=%.4f, se2=%.4f, grad=%.2e, cond=%.2e, alpha[2]=%.3f\n\n",
        se_res$mu, ifelse(is.na(se_res$se2), NA, se_res$se2),
        res$grad_norm, res$cond, res$estimates[2]))
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n\n", e$message))
  })
}
