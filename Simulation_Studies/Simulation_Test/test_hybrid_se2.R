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
cond_limit <- 1e8

# ============================================================
# One-step GMM at W=I + SE1
# ============================================================
run_onestep_se1 <- function(dat, design_mat, h_x, link_type, n) {
  r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
  k <- ncol(design_mat); h <- ncol(h_x)
  tol <- 1e-6

  # g and Gamma functions
  g_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    (r_vec / pi_v - 1) * h_x
  }
  Gamma_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    crossprod(h_x * (r_vec * (1 - pi_v) / pi_v), design_mat) / n
  }

  # Minimize G'I*G = ||G||^2 using NR with cond guard
  W_I <- diag(h)
  est <- rep(0, k)
  for (iter in 1:100) {
    g_mat <- g_fn(est)
    G_i <- matrix(.colMeans(g_mat, n, h), h, 1)
    Gamma_i <- Gamma_fn(est)
    gv <- 2 * as.vector(crossprod(Gamma_i, G_i))
    if (max(abs(gv)) < tol) break
    GtG <- crossprod(Gamma_i)
    direction <- tryCatch(-solve(GtG + GtG, gv),
                          error = function(e) -gv * 2.0 / max(abs(gv)))
    if (sqrt(sum(direction^2)) > 2.0)
      direction <- direction * 2.0 / sqrt(sum(direction^2))
    obj_cur <- sum(G_i^2)
    accepted <- FALSE; step <- 1.0
    for (ls in 1:20) {
      cand <- est + step * direction
      G_cand <- matrix(.colMeans(g_fn(cand), n, h), h, 1)
      obj_new <- sum(G_cand^2)
      if (is.finite(obj_new) && obj_new < obj_cur - 1e-4 * step * sum(gv * direction)) {
        accepted <- TRUE; break
      }
      step <- step * 0.5
    }
    if (accepted) est <- cand else break
  }

  # Compute SE1 at est with W = I
  pi_hat <- compute_pi_fn(as.vector(design_mat %*% est))
  mu_hat <- mean(r_vec * y_vec / pi_hat)
  g_mat <- g_fn(est)
  cf <- r_vec * (1 - pi_hat) / pi_hat
  Gamma_hat <- Gamma_fn(est)  # h x k

  # t_Gamma_arr: n x k x h
  t_Gamma_arr <- array(0, dim = c(n, k, h))
  for (j in 1:k) {
    t_Gamma_arr[, j, ] <- cf * design_mat[, j] * h_x
  }

  W_hat <- W_I
  eta_s <- .colMeans(g_mat, n, h)
  GtW <- crossprod(Gamma_hat, W_hat)  # k x h (= t(Gamma) since W=I)
  GtWG <- GtW %*% Gamma_hat

  H_0n <- tryCatch(-solve(GtWG), error = function(e) {
    -solve(GtWG + 1e-8 * diag(k))
  })
  H_1n <- GtW %*% (t(g_mat) - eta_s)

  # Gamma_2 via numerical differentiation
  eps2 <- 1e-5
  Gamma_2 <- matrix(0, k * h, k)
  for (j in 1:k) {
    ej <- rep(0, k); ej[j] <- eps2
    Gp <- Gamma_fn(est + ej)
    Gm <- Gamma_fn(est - ej)
    Gamma_2[, j] <- (as.vector(t(Gp)) - as.vector(t(Gm))) / (2 * eps2)
  }

  M_n <- kronecker(eta_s %*% W_hat, diag(k)) %*% Gamma_2

  # H_2n_2
  W_eta <- as.vector(W_hat %*% eta_s)
  tGamma_W_eta <- as.vector(GtW %*% eta_s)
  t_Gamma_mat <- matrix(t_Gamma_arr, nrow = n * k, ncol = h)
  H_2n_2_vec <- t_Gamma_mat %*% W_eta
  H_2n_2 <- matrix(H_2n_2_vec, nrow = n, ncol = k, byrow = FALSE)
  H_2n_2 <- t(H_2n_2) - tGamma_W_eta

  # H_2n_3
  GtW_g <- GtW %*% t(g_mat)
  g_W_eta <- as.vector(g_mat %*% W_eta)
  H_2n_3 <- matrix(tGamma_W_eta, k, n) - GtW_g * matrix(g_W_eta, k, n, byrow = TRUE)

  # Q
  I_minus_H0M <- diag(k) - H_0n %*% M_n
  Q <- tryCatch(
    solve(I_minus_H0M) %*% H_0n,
    error = function(e) solve(I_minus_H0M + 1e-8 * diag(k)) %*% H_0n
  )

  psi1 <- Q %*% (H_1n + H_2n_2 + H_2n_3)

  # IPW SE
  if (!is.null(link_type) && link_type == "logistic_complement") {
    dot_pi <- -design_mat * (pi_hat * (1 - pi_hat))
  } else {
    dot_pi <- design_mat * (pi_hat * (1 - pi_hat))
  }
  ry_ps_inv2 <- as.vector(r_vec * y_vec * (pi_hat^(-2)))
  H_alpha <- colMeans(dot_pi * ry_ps_inv2)
  mu_iid <- as.vector(t(r_vec / pi_hat * y_vec) - t(H_alpha) %*% psi1)
  se1 <- sqrt(var(mu_iid) / n)

  list(mu = mu_hat, se = se1)
}

# ========== Run for M2 and M3 ==========
for (m_idx in c(2, 3)) {
  cat(sprintf("\n=== Model %d: Package + hybrid SE ===\n", m_idx))
  mu_vals <- se_vals <- rep(NA_real_, n_reps)
  cond_vals <- rep(NA_real_, n_reps)
  se_type <- rep(NA_character_, n_reps)

  for (rep_i in 1:n_reps) {
    dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
    single_ps <- list(
      formula.list = list(ps_spec[["formula.list"]][[m_idx]]),
      h_alpha.list = list(ps_spec[["h_alpha.list"]][[m_idx]]),
      inv_link = ps_spec[["inv_link"]],
      outcome = ps_spec[["outcome"]],
      alpha_init.list = list(NULL),
      optimizer = "constrained_nr"
    )
    tryCatch({
      ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
      gmm_fit <- ebmr$ps_fit.list[[1]]$gmm_fit
      sol_cond <- gmm_fit$opt$solution_cond
      cond_vals[rep_i] <- sol_cond

      if (sol_cond < cond_limit) {
        # Good solution → use SE2 from psi2
        pi_hat <- ebmr$ps_fit.list[[1]]$fitted.values
        design_mat <- ebmr$ps_fit.list[[1]]$design_matrix
        link_type <- ebmr$ps_fit.list[[1]]$link_type
        r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
        mu <- mean(r_vec * y_vec / pi_hat)

        psi2_mat <- gmm_fit$psi2
        if (!is.null(psi2_mat) && is.matrix(psi2_mat) && !any(is.na(psi2_mat))) {
          if (!is.null(link_type) && link_type == "logistic_complement") {
            dot_pi <- -design_mat * (pi_hat * (1 - pi_hat))
          } else {
            dot_pi <- design_mat * (pi_hat * (1 - pi_hat))
          }
          ry_ps_inv2 <- as.vector(r_vec * y_vec * (pi_hat^(-2)))
          H_alpha <- colMeans(dot_pi * ry_ps_inv2)
          mu_iid <- as.vector(t(r_vec/pi_hat*y_vec) - t(H_alpha) %*% psi2_mat)
          se2 <- sqrt(var(mu_iid) / n_val)
          mu_vals[rep_i] <- mu
          se_vals[rep_i] <- se2
          se_type[rep_i] <- "SE2"
        }
      } else {
        # Degenerate solution → fall back to one-step GMM + SE1
        design_mat <- ebmr$ps_fit.list[[1]]$design_matrix
        h_x <- ebmr$ps_fit.list[[1]]$h_x
        link_type <- ebmr$ps_fit.list[[1]]$link_type
        res <- run_onestep_se1(dat, design_mat, h_x, link_type, n_val)
        mu_vals[rep_i] <- res$mu
        se_vals[rep_i] <- res$se
        se_type[rep_i] <- "SE1"
      }
    }, error = function(e) {
      cat(sprintf("  Rep %d ERROR: %s\n", rep_i, e$message))
    })
  }

  valid <- !is.na(mu_vals) & !is.na(se_vals)
  n_se1 <- sum(se_type == "SE1", na.rm = TRUE)
  n_se2 <- sum(se_type == "SE2", na.rm = TRUE)
  cat(sprintf("  Valid: %d/%d, SE1 (one-step): %d, SE2 (iterative): %d\n",
      sum(valid), n_reps, n_se1, n_se2))
  cat(sprintf("  Bias  = %.4f\n", mean(mu_vals[valid]) - mu_true))
  cat(sprintf("  ESD   = %.4f\n", sd(mu_vals[valid])))
  cat(sprintf("  ESE   = %.4f\n", mean(se_vals[valid])))
  cat(sprintf("  ESE/ESD = %.3f\n", mean(se_vals[valid]) / sd(mu_vals[valid])))

  se1_idx <- valid & se_type == "SE1"
  se2_idx <- valid & se_type == "SE2"
  if (sum(se1_idx) > 1) {
    cat(sprintf("  [SE1 reps] n=%d, mean_SE=%.4f, mu_range=[%.3f,%.3f]\n",
        sum(se1_idx), mean(se_vals[se1_idx]),
        min(mu_vals[se1_idx]), max(mu_vals[se1_idx])))
  }
  if (sum(se2_idx) > 1) {
    cat(sprintf("  [SE2 reps] n=%d, Bias=%.4f, ESD=%.4f, ESE=%.4f, ESE/ESD=%.3f\n",
        sum(se2_idx), mean(mu_vals[se2_idx]) - mu_true, sd(mu_vals[se2_idx]),
        mean(se_vals[se2_idx]), mean(se_vals[se2_idx]) / sd(mu_vals[se2_idx])))
  }
}
