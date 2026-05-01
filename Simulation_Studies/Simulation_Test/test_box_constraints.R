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

# ========== Step 1: Collect normal-rep alpha range ==========
cat("=== Collecting alpha range from normal reps ===\n")
alpha_all <- matrix(NA, n_reps, 4)
cond_all <- rep(NA, n_reps)

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
    alpha_all[rep_i, ] <- ebmr$ps_fit.list[[1]]$gmm_fit$estimates
    cond_all[rep_i] <- ebmr$ps_fit.list[[1]]$gmm_fit$opt$solution_cond
  }, error = function(e) NULL)
}

normal <- which(!is.na(cond_all) & cond_all < 1e8)
cat(sprintf("Normal reps: %d\n", length(normal)))
for (j in 1:4) {
  cat(sprintf("  alpha[%d]: range [%.3f, %.3f], mean=%.3f, sd=%.3f\n",
      j, min(alpha_all[normal, j]), max(alpha_all[normal, j]),
      mean(alpha_all[normal, j]), sd(alpha_all[normal, j])))
}

# Set box bounds: mean +/- 5*sd for each alpha component
alpha_mean <- colMeans(alpha_all[normal, ])
alpha_sd <- apply(alpha_all[normal, ], 2, sd)
lower_bounds <- alpha_mean - 5 * alpha_sd
upper_bounds <- alpha_mean + 5 * alpha_sd
cat(sprintf("\nBox bounds (mean +/- 5*sd):\n"))
for (j in 1:4) {
  cat(sprintf("  alpha[%d]: [%.3f, %.3f]\n", j, lower_bounds[j], upper_bounds[j]))
}

# ========== Step 2: L-BFGS-B with box constraints ==========
cat("\n=== Running box-constrained L-BFGS-B for all reps ===\n")

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

  k <- length(alpha_hat)
  h <- ncol(h_x)

  # H2 computation
  GtW <- crossprod(Gamma_hat, W_hat)
  G_vec <- colMeans(g_mat)
  dc_deta <- cf
  dc_X <- dc_deta * design_mat
  R_mat <- matrix(0, h * k, k)
  for (j_idx in 1:k) {
    Xj_h <- design_mat[, j_idx] * h_x
    block <- crossprod(dc_X, Xj_h) / n
    for (l in 1:h) R_mat[(l-1)*k + j_idx, ] <- block[, l]
  }

  eta_s <- matrix(G_vec, h, 1)
  WG <- as.vector(W_hat %*% eta_s)
  GW_row <- as.vector(t(eta_s) %*% W_hat)

  S_mat <- matrix(0, h^2, k)
  for (j_idx in 1:k) {
    dg_j <- cf * design_mat[, j_idx] * h_x
    cross <- (crossprod(dg_j, g_mat) + crossprod(g_mat, dg_j)) / n
    S_mat[, j_idx] <- as.vector(cross)
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

  # Q_i and psi2
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

  list(mu = mu_hat, se2 = se2, cond = cond_val)
}

# Box-constrained L-BFGS-B iterative GMM
run_box_lbfgsb <- function(dat, design_mat, h_x, link_type, lower, upper) {
  r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
  n <- nrow(dat)
  k <- ncol(design_mat)
  h_dim <- ncol(h_x)

  g_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    (r_vec / pi_v - 1) * h_x
  }
  Gamma_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    crossprod(h_x * (r_vec * (1 - pi_v) / pi_v), design_mat) / n
  }

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

  # Clamp init to bounds
  init <- pmax(pmin(rep(0, k), upper), lower)

  # Step 1: W = I
  est <- init
  st$W <- diag(h_dim)
  opt0 <- tryCatch(
    optim(est, of, gr = gf, method = "L-BFGS-B",
          lower = lower, upper = upper, control = list(maxit = 1000)),
    error = function(e) list(par = est))
  est <- opt0$par

  # Step 2: Iterative W updates
  for (t in 1:100) {
    cache(est)
    st$W <- tryCatch(W_fn(st$cg), error = function(e) diag(h_dim))
    fg <- max(abs(gf(est)))
    if (fg < 1e-6) break
    st$cp <- NULL
    opt_t <- tryCatch(
      optim(est, of, gr = gf, method = "L-BFGS-B",
            lower = lower, upper = upper, control = list(maxit = 1000)),
      error = function(e) list(par = est))
    est <- opt_t$par
  }

  g_f <- g_fn(est)
  G_f <- matrix(.colMeans(g_f, n, h_dim), h_dim, 1)
  W_f <- tryCatch(W_fn(g_f), error = function(e) diag(h_dim))
  grad_norm <- max(abs(2 * as.vector(crossprod(Gamma_fn(est), W_f %*% G_f))))

  list(estimates = est, grad_norm = grad_norm)
}

# Run for all reps
mu_box <- se2_box <- cond_box <- rep(NA_real_, n_reps)
grad_box <- rep(NA_real_, n_reps)
at_bound <- rep(FALSE, n_reps)

for (rep_i in 1:n_reps) {
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]

  # Get design_mat and h_x from package
  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
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

    res <- run_box_lbfgsb(dat, design_mat, h_x, link_type, lower_bounds, upper_bounds)
    grad_box[rep_i] <- res$grad_norm

    # Check if at boundary
    at_bound[rep_i] <- any(abs(res$estimates - lower_bounds) < 1e-4 |
                           abs(res$estimates - upper_bounds) < 1e-4)

    # Compute SE2
    se_res <- compute_se2_at_alpha(dat, res$estimates, design_mat, h_x, link_type, n_val)
    mu_box[rep_i] <- se_res$mu
    se2_box[rep_i] <- se_res$se2
    cond_box[rep_i] <- se_res$cond
  }, error = function(e) NULL)
}

# ========== Results ==========
valid <- !is.na(mu_box) & !is.na(se2_box)
cat(sprintf("\n=== Box-constrained L-BFGS-B (mean +/- 5*sd) ===\n"))
cat(sprintf("  Valid: %d/%d\n", sum(valid), n_reps))
cat(sprintf("  At boundary: %d\n", sum(at_bound, na.rm = TRUE)))
cat(sprintf("  Bias = %.4f\n", mean(mu_box[valid]) - mu_true))
cat(sprintf("  ESD  = %.4f\n", sd(mu_box[valid])))
cat(sprintf("  ESE2 = %.4f\n", mean(se2_box[valid])))
cat(sprintf("  ESE2/ESD = %.3f\n", mean(se2_box[valid]) / sd(mu_box[valid])))

# Compare with package (unconstrained)
cat(sprintf("\n=== For reference: package constrained_nr ===\n"))
mu_pkg <- se2_pkg <- cond_pkg_v <- rep(NA_real_, n_reps)
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
    pi_hat <- ebmr$ps_fit.list[[1]]$fitted.values
    design_mat <- ebmr$ps_fit.list[[1]]$design_matrix
    link_type <- ebmr$ps_fit.list[[1]]$link_type
    r_vec <- dat[["r"]]; y_vec <- dat[["y"]]

    mu_pkg[rep_i] <- mean(r_vec * y_vec / pi_hat)
    cond_pkg_v[rep_i] <- gmm_fit$opt$solution_cond

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
      se2_pkg[rep_i] <- sqrt(var(mu_iid) / n_val)
    }
  }, error = function(e) NULL)
}

valid_pkg <- !is.na(mu_pkg) & !is.na(se2_pkg)
cat(sprintf("  Valid: %d/%d\n", sum(valid_pkg), n_reps))
cat(sprintf("  Bias = %.4f\n", mean(mu_pkg[valid_pkg]) - mu_true))
cat(sprintf("  ESD  = %.4f\n", sd(mu_pkg[valid_pkg])))
cat(sprintf("  ESE2 = %.4f\n", mean(se2_pkg[valid_pkg])))
cat(sprintf("  ESE2/ESD = %.3f\n", mean(se2_pkg[valid_pkg]) / sd(mu_pkg[valid_pkg])))

# Detail on formerly-degenerate reps
degen <- which(!is.na(cond_pkg_v) & cond_pkg_v >= 1e8)
cat(sprintf("\n=== Detail on %d formerly-degenerate reps ===\n", length(degen)))
cat(sprintf("  %4s  %8s  %8s  %8s  %8s  %8s  %5s\n",
    "rep", "mu_box", "se2_box", "cond_box", "mu_pkg", "se2_pkg", "bound"))
for (i in degen) {
  cat(sprintf("  %4d  %8.4f  %8.4f  %8.2e  %8.4f  %8.4f  %5s\n",
      i, mu_box[i], ifelse(is.na(se2_box[i]), NA, se2_box[i]),
      cond_box[i], mu_pkg[i],
      ifelse(is.na(se2_pkg[i]), NA, se2_pkg[i]),
      ifelse(at_bound[i], "YES", "no")))
}
