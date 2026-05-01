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
cond_limit <- 1e8
m_idx <- 3
delta <- 0.01  # cap: pi <= 1 - delta

# GMM with capped pi INSIDE the moment function
# g_i = (r_i / min(pi_i, 1-delta) - 1) * h(X_i)
# This prevents g_i from going to zero when pi -> 1
run_capped_gmm <- function(dat, design_mat, h_x, link_type, delta) {
  r_vec <- dat[["r"]]
  n <- nrow(dat); k <- ncol(design_mat); h_dim <- ncol(h_x)

  # Raw pi (uncapped)
  pi_raw_fn <- function(alpha) compute_pi_fn(as.vector(design_mat %*% alpha))

  # Capped pi for use in moment function
  pi_cap_fn <- function(alpha) pmin(pi_raw_fn(alpha), 1 - delta)

  # Moment function with capped pi
  g_fn <- function(alpha) {
    pi_c <- pi_cap_fn(alpha)
    (r_vec / pi_c - 1) * h_x
  }

  # Gamma with cap: derivative is 0 for observations where pi > 1-delta
  Gamma_fn <- function(alpha) {
    pi_raw <- pi_raw_fn(alpha)
    pi_c <- pmin(pi_raw, 1 - delta)
    # For capped obs (pi_raw > 1-delta), the derivative of g w.r.t. alpha is 0
    # For uncapped obs, dg/dalpha = -(r / pi^2) * dpi/dalpha * h_x
    uncapped <- pi_raw <= (1 - delta)
    cf <- ifelse(uncapped, r_vec * (1 - pi_c) / pi_c, 0)
    crossprod(h_x * cf, design_mat) / n
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

  nr_inner <- function(start, W_hat, trust_r = 2.0) {
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

  # Iterative GMM
  est <- rep(0, k)
  best_grad <- Inf; best_est <- est

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

  # Return: use RAW pi (uncapped) for the final IPW estimate
  final_cond <- forward_cond(best_est)
  final_grad <- max(abs(grad_fn(best_est,
    tryCatch(W_fn(g_fn(best_est)), error = function(e) diag(h_dim)))))
  pi_final <- pi_raw_fn(best_est)
  list(estimates = best_est, grad_norm = final_grad, cond = final_cond,
       pi_hat = pi_final)
}

# SE computation using raw (uncapped) pi
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

# ===== Stubborn reps =====
stubborn_reps <- c(35, 38, 143, 155, 156, 190)

for (d in c(0.01, 0.05)) {
  cat(sprintf("\n=== Stubborn reps: pi cap at 1-%.2f = %.2f ===\n", d, 1-d))
  for (rep_i in stubborn_reps) {
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

      res <- run_capped_gmm(dat, design_mat, h_x, link_type, delta = d)
      r_vec <- dat$r; y_vec <- dat$y
      mu <- mean(r_vec * y_vec / res$pi_hat)
      cat(sprintf("  Rep %3d: grad=%.2e, cond=%.2e, mu=%.4f, alpha[2]=%.3f, #pi>%.2f=%d, #pi>0.99=%d\n",
          rep_i, res$grad_norm, res$cond, mu, res$estimates[2],
          1-d, sum(res$pi_hat > 1-d), sum(res$pi_hat > 0.99)))
    }, error = function(e) {
      cat(sprintf("  Rep %3d: ERROR %s\n", rep_i, e$message))
    })
  }
}

# ===== All 200 with delta=0.05 =====
cat("\n\n=== All 200 reps: M3, pi capped at 0.95 in moment fn ===\n")
mu_vals <- se_vals <- cond_vals <- grad_vals <- rep(NA_real_, n_reps)

for (rep_i in 1:n_reps) {
  if (rep_i %% 50 == 0) cat(sprintf("  rep %d...\n", rep_i))
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

    res <- run_capped_gmm(dat, design_mat, h_x, link_type, delta = 0.05)
    cond_vals[rep_i] <- res$cond
    grad_vals[rep_i] <- res$grad_norm

    se_res <- compute_se(dat, res$estimates, design_mat, h_x, link_type, n_val)
    mu_vals[rep_i] <- se_res$mu
    se_vals[rep_i] <- se_res$se
  }, error = function(e) {
    cat(sprintf("  Rep %d ERROR: %s\n", rep_i, e$message))
  })
}

valid <- !is.na(mu_vals) & !is.na(se_vals)
n_degen <- sum(cond_vals[valid] >= cond_limit, na.rm = TRUE)
n_conv <- sum(grad_vals[valid] < tol, na.rm = TRUE)
cat(sprintf("\n  Valid: %d/%d, Degenerate: %d, Converged: %d\n",
    sum(valid), n_reps, n_degen, n_conv))
cat(sprintf("  Bias     = %.4f\n", mean(mu_vals[valid]) - mu_true))
cat(sprintf("  ESD      = %.4f\n", sd(mu_vals[valid])))
cat(sprintf("  ESE      = %.4f\n", mean(se_vals[valid])))
cat(sprintf("  ESE/ESD  = %.3f\n", mean(se_vals[valid]) / sd(mu_vals[valid])))
ci_lo <- mu_vals[valid] - 1.96 * se_vals[valid]
ci_hi <- mu_vals[valid] + 1.96 * se_vals[valid]
cp <- mean(mu_true >= ci_lo & mu_true <= ci_hi)
cat(sprintf("  CP       = %.3f\n", cp))

well_cond <- valid & cond_vals < cond_limit
if (sum(well_cond) > 1 && sum(well_cond) < sum(valid)) {
  cat(sprintf("\n  [Cond<1e8] n=%d, Bias=%.4f, ESD=%.4f, ESE=%.4f, ESE/ESD=%.3f\n",
      sum(well_cond), mean(mu_vals[well_cond]) - mu_true, sd(mu_vals[well_cond]),
      mean(se_vals[well_cond]), mean(se_vals[well_cond]) / sd(mu_vals[well_cond])))
}
