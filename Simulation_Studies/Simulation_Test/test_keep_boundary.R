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

model_idx <- 3

# Strategy: run CNR with NO fallback (just use constrained_nr inner logic),
# then separately run the full package (with fallback) for comparison.
# To simulate "no fallback", we can't easily modify the package,
# but we CAN reconstruct SE2 at the CNR boundary solution manually.

compute_se2_at_alpha <- function(dat, alpha_hat, design_mat, h_x, link_type, n) {
  r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
  compute_pi <- function(eta) plogis(-eta)
  pi_hat <- compute_pi(as.vector(design_mat %*% alpha_hat))
  mu_hat <- mean(r_vec * y_vec / pi_hat)

  # Build g, Gamma, W at this alpha
  g_mat <- (r_vec / pi_hat - 1) * h_x
  G_vec <- colMeans(g_mat)
  W_hat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) NULL)
  if (is.null(W_hat)) return(list(mu = mu_hat, se2 = NA, cond = NA, foc = NA))

  cf <- r_vec * (1 - pi_hat) / pi_hat
  Gamma_hat <- crossprod(h_x * cf, design_mat) / n
  GtW <- crossprod(Gamma_hat, W_hat)
  GtWG <- GtW %*% Gamma_hat
  foc <- max(abs(GtW %*% G_vec))

  eig <- eigen(GtWG, symmetric = TRUE, only.values = TRUE)$values
  cond_val <- max(eig) / max(min(eig), 1e-15)

  k <- length(alpha_hat)
  h <- ncol(h_x)

  # H2 computation
  dc_deta <- cf
  dc_X <- dc_deta * design_mat
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

  if (!is.finite(H2_cond) || H2_cond > 1e15) return(list(mu = mu_hat, se2 = NA, cond = cond_val, foc = foc))

  H2_inv <- tryCatch(solve(H2_mat), error = function(e) NULL)
  if (is.null(H2_inv)) return(list(mu = mu_hat, se2 = NA, cond = cond_val, foc = foc))

  # Q_i
  GtW_g <- GtW %*% t(g_mat)
  g_WG <- as.vector(g_mat %*% WG)
  cf_hxWG <- cf * as.vector(h_x %*% WG)
  GammaT_WG <- t(design_mat * cf_hxWG)
  Q_mat <- GtW_g + GammaT_WG - sweep(GtW_g, 2, g_WG, `*`)

  psi2 <- -(H2_inv %*% Q_mat)

  # SE2 for mu
  if (!is.null(link_type) && link_type == "logistic_complement") {
    dot_pi <- -design_mat * (pi_hat * (1 - pi_hat))
  } else {
    dot_pi <- design_mat * (pi_hat * (1 - pi_hat))
  }
  ry_ps_inv2 <- as.vector(r_vec * y_vec * (pi_hat^(-2)))
  H_alpha <- colMeans(dot_pi * ry_ps_inv2)
  mu_iid <- as.vector(t(r_vec/pi_hat*y_vec) - t(H_alpha) %*% psi2)
  se2 <- sqrt(var(mu_iid) / n)

  list(mu = mu_hat, se2 = se2, cond = cond_val, foc = foc, H2_cond = H2_cond)
}

# For each rep, get both:
# (a) Package result (with L-BFGS-B fallback)
# (b) CNR-only boundary result (reconstructed from trace data)

# We need the CNR boundary alpha. Run the manual CNR loop from trace_cnr_path.R
# but simplified: just run the package CNR and capture the boundary solution
# before the fallback overwrites it.
# Easiest: run with optimizer="constrained_nr", but also run the manual CNR
# to get the pre-fallback solution.

# Actually, simpler approach: the trace showed the best_estimates for each
# degenerate rep. Let me just compute SE2 at those points.
# But we need the actual alpha values. Let me re-run the manual CNR for all reps.

mu_pkg <- se2_pkg <- cond_pkg <- rep(NA_real_, n_reps)
mu_bnd <- se2_bnd <- cond_bnd <- foc_bnd <- rep(NA_real_, n_reps)
is_boundary <- rep(FALSE, n_reps)

compute_pi_fn <- function(eta) plogis(-eta)

for (rep_i in 1:n_reps) {
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
  r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
  n <- n_val

  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL),
    optimizer = "constrained_nr"
  )

  tryCatch({
    # Package result (includes fallback)
    ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
    gmm_fit <- ebmr$ps_fit.list[[1]]$gmm_fit
    design_mat <- ebmr$ps_fit.list[[1]]$design_matrix
    h_x <- ebmr$ps_fit.list[[1]]$h_x
    link_type <- ebmr$ps_fit.list[[1]]$link_type

    pi_hat <- ebmr$ps_fit.list[[1]]$fitted.values
    mu_pkg[rep_i] <- mean(r_vec * y_vec / pi_hat)
    cond_pkg[rep_i] <- gmm_fit$opt$solution_cond

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
      se2_pkg[rep_i] <- sqrt(var(mu_iid) / n)
    }

    # Check if this is a degenerate solution (fallback was used)
    if (cond_pkg[rep_i] >= 1e8) {
      is_boundary[rep_i] <- TRUE

      # Run manual CNR without fallback to get boundary solution
      k <- ncol(design_mat)
      h_dim <- ncol(h_x)

      g_fn <- function(alpha) (r_vec / compute_pi_fn(as.vector(design_mat %*% alpha)) - 1) * h_x
      Gamma_fn <- function(alpha) {
        pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
        crossprod(h_x * (r_vec * (1 - pi_v) / pi_v), design_mat) / n
      }

      nr_inner_manual <- function(alpha, W_hat, cond_limit = 1e8) {
        for (iter in 1:50) {
          g_i <- g_fn(alpha)
          G_i <- matrix(colMeans(g_i), h_dim, 1)
          Gamma_i <- Gamma_fn(alpha)
          g_vec <- 2 * as.vector(crossprod(Gamma_i, W_hat %*% G_i))
          if (max(abs(g_vec)) < 1e-6) break
          GtWG <- crossprod(Gamma_i, W_hat %*% Gamma_i)
          H_mat <- GtWG + GtWG
          direction <- tryCatch(-solve(H_mat, g_vec),
                                error = function(e) -g_vec * 2.0 / max(abs(g_vec)))
          if (sqrt(sum(direction^2)) > 2.0)
            direction <- direction * 2.0 / sqrt(sum(direction^2))
          obj_cur <- as.numeric(crossprod(G_i, W_hat %*% G_i))
          accepted <- FALSE
          step <- 1.0
          for (ls in 1:20) {
            cand <- alpha + step * direction
            g_cand <- g_fn(cand)
            G_cand <- matrix(colMeans(g_cand), h_dim, 1)
            obj_new <- as.numeric(crossprod(G_cand, W_hat %*% G_cand))
            if (is.finite(obj_new) && obj_new < obj_cur - 1e-4 * step * sum(g_vec * direction)) {
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

      estimates <- nr_inner_manual(rep(0, k), diag(h_dim))
      best_grad <- Inf; best_est <- estimates
      for (t in 1:100) {
        g_mat <- g_fn(estimates)
        W_hat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(h_dim))
        Gamma_hat <- Gamma_fn(estimates)
        G_vec <- matrix(colMeans(g_mat), h_dim, 1)
        grad <- max(abs(2 * as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))))
        if (grad < 1e-6) { best_est <- estimates; best_grad <- grad; break }
        if (grad < best_grad) { best_grad <- grad; best_est <- estimates }
        estimates <- nr_inner_manual(estimates, W_hat)
      }

      # Compute SE2 at boundary solution
      res <- compute_se2_at_alpha(dat, best_est, design_mat, h_x, link_type, n)
      mu_bnd[rep_i] <- res$mu
      se2_bnd[rep_i] <- res$se2
      cond_bnd[rep_i] <- res$cond
      foc_bnd[rep_i] <- res$foc
    }
  }, error = function(e) NULL)
}

# Report
valid_pkg <- !is.na(mu_pkg) & !is.na(se2_pkg)
cat(sprintf("=== Package (with fallback): all valid ===\n"))
cat(sprintf("  n=%d, Bias=%.4f, ESD=%.4f, ESE2/ESD=%.3f\n",
    sum(valid_pkg), mean(mu_pkg[valid_pkg]) - mu_true, sd(mu_pkg[valid_pkg]),
    mean(se2_pkg[valid_pkg]) / sd(mu_pkg[valid_pkg])))

# Hybrid: use package for non-degenerate, boundary solution for degenerate
mu_hyb <- mu_pkg
se2_hyb <- se2_pkg
for (i in which(is_boundary)) {
  if (!is.na(mu_bnd[i])) mu_hyb[i] <- mu_bnd[i]
  if (!is.na(se2_bnd[i])) se2_hyb[i] <- se2_bnd[i]
}
valid_hyb <- !is.na(mu_hyb) & !is.na(se2_hyb)
cat(sprintf("\n=== Hybrid (boundary for degenerate): all valid ===\n"))
cat(sprintf("  n=%d, Bias=%.4f, ESD=%.4f, ESE2/ESD=%.3f\n",
    sum(valid_hyb), mean(mu_hyb[valid_hyb]) - mu_true, sd(mu_hyb[valid_hyb]),
    mean(se2_hyb[valid_hyb]) / sd(mu_hyb[valid_hyb])))

# Detail for degenerate reps
cat(sprintf("\n=== Degenerate reps detail ===\n"))
cat(sprintf("  %4s  %10s  %10s  %10s  %10s  %10s  %10s\n",
    "rep", "mu_pkg", "se2_pkg", "cond_pkg", "mu_bnd", "se2_bnd", "cond_bnd"))
for (i in which(is_boundary)) {
  cat(sprintf("  %4d  %10.4f  %10.4f  %10.2e  %10.4f  %10s  %10.2e\n",
      i, mu_pkg[i], ifelse(is.na(se2_pkg[i]), NA, se2_pkg[i]), cond_pkg[i],
      mu_bnd[i], ifelse(is.na(se2_bnd[i]), "NA", sprintf("%.4f", se2_bnd[i])),
      cond_bnd[i]))
}
