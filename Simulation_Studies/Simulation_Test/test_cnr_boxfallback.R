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
compute_pi_fn <- function(eta) plogis(-eta)
tol <- 1e-6
max_outer <- 100

# SE2 computation at a given alpha
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

# Simulate constrained_nr + box-constrained L-BFGS-B fallback externally
# CNR part: run the package with constrained_nr, extract cnr boundary solution
# Box fallback: run L-BFGS-B with bounds derived from cnr solution

run_cnr_box_fallback <- function(dat, design_mat, h_x, link_type,
                                  cnr_alpha, cnr_grad, box_radius = 5.0) {
  r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
  n <- nrow(dat); k <- ncol(design_mat); h_dim <- ncol(h_x)

  g_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    (r_vec / pi_v - 1) * h_x
  }
  Gamma_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    crossprod(h_x * (r_vec * (1 - pi_v) / pi_v), design_mat) / n
  }

  # Box bounds centered on cnr solution
  lower <- cnr_alpha - box_radius
  upper <- cnr_alpha + box_radius

  st <- new.env(parent = emptyenv())
  st$W <- NULL; st$cp <- NULL; st$cg <- NULL; st$cG <- NULL
  cache <- function(p) {
    if (!is.null(st$cp) && identical(p, st$cp)) return()
    st$cp <- p; st$cg <- g_fn(p)
    st$cG <- matrix(.colMeans(st$cg, n, h_dim), h_dim, 1)
  }
  of <- function(p) {
    cache(p)
    v <- as.numeric(crossprod(st$cG, st$W %*% st$cG))
    if (!is.finite(v)) 1e8 else v
  }
  gf <- function(p) {
    cache(p)
    2 * as.vector(crossprod(Gamma_fn(p), st$W %*% st$cG))
  }

  run_one <- function(start) {
    est <- pmax(pmin(start, upper), lower)
    st$W <<- diag(h_dim)
    st$cp <<- NULL
    opt0 <- tryCatch(
      optim(est, of, gr = gf, method = "L-BFGS-B",
            lower = lower, upper = upper, control = list(maxit = 1000)),
      error = function(e) list(par = est))
    est <- opt0$par
    for (t in 1:max_outer) {
      cache(est)
      st$W <<- tryCatch(W_fn(st$cg), error = function(e) diag(h_dim))
      fg <- max(abs(gf(est)))
      if (fg < tol) break
      st$cp <<- NULL
      opt_t <- tryCatch(
        optim(est, of, gr = gf, method = "L-BFGS-B",
              lower = lower, upper = upper, control = list(maxit = 1000)),
        error = function(e) list(par = est))
      est <- opt_t$par
    }
    g_f <- g_fn(est)
    G_f <- matrix(.colMeans(g_f, n, h_dim), h_dim, 1)
    W_f <- tryCatch(W_fn(g_f), error = function(e) diag(h_dim))
    gn <- max(abs(2 * as.vector(crossprod(Gamma_fn(est), W_f %*% G_f))))
    list(estimates = est, grad_norm = gn)
  }

  fb1 <- run_one(cnr_alpha)
  fb2 <- run_one(rep(0, k))
  fb <- if (fb1$grad_norm <= fb2$grad_norm) fb1 else fb2
  if (fb$grad_norm < cnr_grad) fb else list(estimates = cnr_alpha, grad_norm = cnr_grad)
}

# ========== Main loop ==========
mu_res <- se2_res <- cond_res <- rep(NA_real_, n_reps)
used_box <- rep(FALSE, n_reps)

for (rep_i in 1:n_reps) {
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]

  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL),
    optimizer = "constrained_nr"
  )
  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
    gmm_fit <- ebmr$ps_fit.list[[1]]$gmm_fit
    design_mat <- ebmr$ps_fit.list[[1]]$design_matrix
    h_x <- ebmr$ps_fit.list[[1]]$h_x
    link_type <- ebmr$ps_fit.list[[1]]$link_type
    pkg_cond <- gmm_fit$opt$solution_cond
    pkg_alpha <- gmm_fit$estimates
    pkg_grad <- gmm_fit$opt$final_grad_norm

    if (pkg_cond < 1e8) {
      # Normal rep: use package result directly
      pi_hat <- ebmr$ps_fit.list[[1]]$fitted.values
      r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
      mu_res[rep_i] <- mean(r_vec * y_vec / pi_hat)
      cond_res[rep_i] <- pkg_cond
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
        se2_res[rep_i] <- sqrt(var(mu_iid) / n_val)
      }
    } else {
      # Degenerate rep: use box-constrained L-BFGS-B fallback
      # Need CNR boundary alpha (before fallback overwrote it).
      # Re-run manual CNR to get boundary solution.
      r_vec <- dat[["r"]]; n <- n_val; k <- ncol(design_mat); h_dim <- ncol(h_x)
      g_fn <- function(alpha) {
        pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
        (r_vec / pi_v - 1) * h_x
      }
      Gamma_fn <- function(alpha) {
        pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
        crossprod(h_x * (r_vec * (1 - pi_v) / pi_v), design_mat) / n
      }
      nr_inner <- function(alpha, W_hat, cond_limit = 1e8) {
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

      est <- nr_inner(rep(0, k), diag(h_dim))
      best_grad <- Inf; best_est <- est
      for (t in 1:max_outer) {
        g_mat <- g_fn(est)
        W_hat <- tryCatch(W_fn(g_mat), error = function(e) diag(h_dim))
        Gamma_hat <- Gamma_fn(est)
        G_vec <- matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
        grad <- max(abs(2 * as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))))
        if (grad < tol) { best_est <- est; best_grad <- grad; break }
        if (grad < best_grad) { best_grad <- grad; best_est <- est }
        est <- nr_inner(est, W_hat)
      }
      cnr_alpha <- best_est; cnr_grad <- best_grad

      # Now run box-constrained L-BFGS-B fallback
      fb <- run_cnr_box_fallback(dat, design_mat, h_x, link_type,
                                  cnr_alpha, cnr_grad, box_radius = 2.0)
      used_box[rep_i] <- TRUE

      se_res <- compute_se2_at_alpha(dat, fb$estimates, design_mat, h_x, link_type, n_val)
      mu_res[rep_i] <- se_res$mu
      se2_res[rep_i] <- se_res$se2
      cond_res[rep_i] <- se_res$cond
    }
  }, error = function(e) NULL)
}

# ========== Report ==========
valid <- !is.na(mu_res) & !is.na(se2_res)
cat(sprintf("\n=== CNR + box-constrained fallback (radius=2) ===\n"))
cat(sprintf("  Valid: %d/%d, box used: %d\n", sum(valid), n_reps, sum(used_box)))
cat(sprintf("  Bias = %.4f\n", mean(mu_res[valid]) - mu_true))
cat(sprintf("  ESD  = %.4f\n", sd(mu_res[valid])))
cat(sprintf("  ESE2 = %.4f\n", mean(se2_res[valid])))
cat(sprintf("  ESE2/ESD = %.3f\n", mean(se2_res[valid]) / sd(mu_res[valid])))

# Detail on box-used reps
cat(sprintf("\n=== Detail on box-fallback reps ===\n"))
cat(sprintf("  %4s  %8s  %8s  %10s\n", "rep", "mu", "se2", "cond"))
for (i in which(used_box)) {
  cat(sprintf("  %4d  %8.4f  %8s  %10.2e\n",
      i, mu_res[i],
      ifelse(is.na(se2_res[i]), "NA", sprintf("%.4f", se2_res[i])),
      cond_res[i]))
}

# Also report excluding degenerate for comparison
normal_idx <- valid & !used_box
cat(sprintf("\n=== Normal reps only (for reference) ===\n"))
cat(sprintf("  n=%d, Bias=%.4f, ESD=%.4f, ESE2=%.4f, ESE2/ESD=%.3f\n",
    sum(normal_idx), mean(mu_res[normal_idx]) - mu_true,
    sd(mu_res[normal_idx]), mean(se2_res[normal_idx]),
    mean(se2_res[normal_idx]) / sd(mu_res[normal_idx])))
