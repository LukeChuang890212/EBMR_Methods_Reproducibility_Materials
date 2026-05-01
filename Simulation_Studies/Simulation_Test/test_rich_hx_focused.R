## Focused test: Enriched h(X) on problematic reps only
## Step 1: Identify degenerate reps with standard h_x (h_dim=5)
## Step 2: Re-run those reps with enriched h_x (h_dim=10) and diagnose SE
## Focus on S4 M3 (the main problem case)
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

n_val <- 2000; n_reps <- 200  # only first 200 reps for speed
ps_spec <- get_ps_spec("9-alt1")
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
compute_pi_fn <- function(eta) plogis(-eta)
cond_limit <- 1e8

## GMM optimizer (same as before)
run_gmm <- function(dat, design_mat, h_x, link_type) {
  r_vec <- dat[["r"]]
  n <- nrow(dat); k <- ncol(design_mat); h_dim <- ncol(h_x)
  init <- rep(0, k)
  tol <- 1e-6; max_outer <- 1000L; trust_radius <- 0.5

  g_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    (r_vec / pi_v - 1) * h_x
  }
  Gamma_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    crossprod(h_x * (r_vec * (1 - pi_v) / pi_v), design_mat) / n
  }

  nr_inner <- function(start, W_hat) {
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
      if (sn > trust_radius) direction <- direction * (trust_radius / sn)
      obj_cur <- as.numeric(crossprod(G_vec, W_hat %*% G_vec))
      step <- 1.0; accepted <- FALSE
      for (ls in 1:30) {
        cand <- alpha + step * direction
        g_cand <- g_fn(cand)
        G_cand <- matrix(.colMeans(g_cand, n, h_dim), h_dim, 1)
        obj_new <- as.numeric(crossprod(G_cand, W_hat %*% G_cand))
        if (is.finite(obj_new) && obj_new < obj_cur - 1e-4 * step * sum(gv * direction)) {
          accepted <- TRUE; break
        }
        step <- step * 0.5
      }
      if (accepted) alpha <- cand else break
    }
    alpha
  }

  estimates <- nr_inner(init, diag(h_dim))
  for (t in 1:max_outer) {
    g_mat <- g_fn(estimates)
    G_vec <- matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
    W_hat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(h_dim))
    Gamma <- Gamma_fn(estimates)
    cg <- 2 * as.vector(crossprod(Gamma, W_hat %*% G_vec))
    if (max(abs(cg)) < tol) break
    H_gn <- crossprod(Gamma, W_hat %*% Gamma)
    direction <- tryCatch(solve(H_gn, -cg), error = function(e) -cg)
    sn <- sqrt(sum(direction^2))
    if (sn > trust_radius) direction <- direction * (trust_radius / sn)
    obj_cur <- as.numeric(crossprod(G_vec, W_hat %*% G_vec))
    step <- 1.0; accepted <- FALSE
    for (ls in 1:30) {
      cand <- estimates + step * direction
      g_cand <- g_fn(cand)
      G_cand <- matrix(.colMeans(g_cand, n, h_dim), h_dim, 1)
      obj_new <- as.numeric(crossprod(G_cand, W_hat %*% G_cand))
      if (is.finite(obj_new) && obj_new < obj_cur - 1e-4 * step * sum(cg * direction)) {
        accepted <- TRUE; break
      }
      step <- step * 0.5
    }
    if (accepted) estimates <- cand
  }

  g_mat_f <- g_fn(estimates)
  W_f <- tryCatch(solve(crossprod(g_mat_f) / n), error = function(e) diag(h_dim))
  pi_final <- compute_pi_fn(as.vector(design_mat %*% estimates))
  cond_val <- tryCatch({
    Gamma_f <- Gamma_fn(estimates)
    H <- crossprod(Gamma_f, W_f %*% Gamma_f)
    ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    max(ev) / max(min(ev), 1e-15)
  }, error = function(e) Inf)
  list(estimates = estimates, pi_hat = pi_final, cond = cond_val)
}

## SE with detailed diagnostics for debugging
compute_se_debug <- function(dat, alpha_hat, design_mat, h_x, link_type, n) {
  r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
  pi_hat <- compute_pi_fn(as.vector(design_mat %*% alpha_hat))
  mu_hat <- mean(r_vec * y_vec / pi_hat)
  g_mat <- (r_vec / pi_hat - 1) * h_x

  # W matrix
  W_cross <- crossprod(g_mat) / n
  W_cond <- tryCatch({
    ev <- eigen(W_cross, symmetric = TRUE, only.values = TRUE)$values
    max(ev) / max(min(ev), 1e-15)
  }, error = function(e) Inf)
  W_hat <- tryCatch(solve(W_cross), error = function(e) NULL)
  if (is.null(W_hat)) return(list(mu = mu_hat, se = NA,
    diag = list(fail = "W_hat singular", W_cond = W_cond)))

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

  # Detailed H_mat diagnostics
  H_ev <- tryCatch(eigen(H_mat, symmetric = FALSE, only.values = TRUE)$values,
                   error = function(e) rep(NA, k))
  H_cond <- max(Mod(H_ev)) / max(min(Mod(H_ev)), 1e-15)

  # Check individual components
  GtWG_cond <- tryCatch({
    ev <- eigen(GtWG, symmetric = TRUE, only.values = TRUE)$values
    max(ev) / max(min(ev), 1e-15)
  }, error = function(e) Inf)

  # Norms of the correction terms
  correction_R <- GW_kron_Ik %*% R_mat
  correction_S <- GW_kron_GtW %*% S_mat
  norm_GtWG <- max(abs(GtWG))
  norm_corrR <- max(abs(correction_R))
  norm_corrS <- max(abs(correction_S))
  norm_G <- max(abs(G_vec))

  diag_info <- list(
    W_cond = W_cond,
    GtWG_cond = GtWG_cond,
    H_cond = H_cond,
    H_ev = H_ev,
    norm_GtWG = norm_GtWG,
    norm_corrR = norm_corrR,
    norm_corrS = norm_corrS,
    norm_G = norm_G,
    fail = "none"
  )

  if (!is.finite(H_cond) || H_cond > 1e15) {
    diag_info$fail <- paste0("H_cond=", format(H_cond, scientific=TRUE))
    return(list(mu = mu_hat, se = NA, diag = diag_info))
  }

  H_inv <- tryCatch(solve(H_mat), error = function(e) NULL)
  if (is.null(H_inv)) {
    diag_info$fail <- "H_inv singular"
    return(list(mu = mu_hat, se = NA, diag = diag_info))
  }

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
  diag_info$fail <- "none"
  list(mu = mu_hat, se = se, diag = diag_info)
}

## ============================================================
## STEP 1: Find degenerate reps with standard h_x
## ============================================================
cat("=== Step 1: Identify degenerate reps with standard h_x (S4 M3) ===\n\n")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
mu_true <- get_mu_true("setting4")
cat(sprintf("True mean: %.6f\n", mu_true))

m_idx <- 3  # Model 3: r ~ y + u2 + z2
degen_reps <- c()

for (rep_i in 1:n_reps) {
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[m_idx]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[m_idx]]),
    inv_link = ps_spec[["inv_link"]], outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL), optimizer = "L-BFGS-B"
  )
  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
    dm <- ebmr$ps_fit.list[[1]]$design_matrix
    hx <- ebmr$ps_fit.list[[1]]$h_x
    lt <- ebmr$ps_fit.list[[1]]$link_type
    res <- run_gmm(dat, dm, hx, lt)
    if (res$cond >= cond_limit) degen_reps <- c(degen_reps, rep_i)
  }, error = function(e) NULL)
}
cat(sprintf("Found %d degenerate reps out of %d: %s\n\n",
    length(degen_reps), n_reps,
    paste(head(degen_reps, 20), collapse=", ")))

## ============================================================
## STEP 2: Run degenerate reps with enriched h_x
## h_rich = (1, u1, u2, z1, z2, u1^2, u2^2, z1^2, z2^2, u1*u2)
## ============================================================
cat("=== Step 2: Re-run degenerate reps with enriched h_x ===\n\n")

# Also run some non-degenerate reps as control
non_degen_reps <- setdiff(1:n_reps, degen_reps)
control_reps <- head(non_degen_reps, min(20, length(non_degen_reps)))
test_reps <- c(degen_reps, control_reps)

cat(sprintf("Testing %d degenerate + %d control = %d reps\n",
    length(degen_reps), length(control_reps), length(test_reps)))

for (rep_i in test_reps) {
  is_degen <- rep_i %in% degen_reps
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]

  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[m_idx]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[m_idx]]),
    inv_link = ps_spec[["inv_link"]], outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL), optimizer = "L-BFGS-B"
  )

  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
    dm <- ebmr$ps_fit.list[[1]]$design_matrix  # (Intercept, y, u2, z2)
    lt <- ebmr$ps_fit.list[[1]]$link_type

    # Build enriched h_x: (1, u1, u2, z1, z2, u2^2, z2^2, u1*u2)
    # Note: u1 and z1 are binary so u1^2=u1, z1^2=z1 — skip those
    h_rich <- cbind(
      1,
      dat$u1, dat$u2, dat$z1, dat$z2,
      dat$u2^2, dat$z2^2,
      dat$u1 * dat$u2
    )
    colnames(h_rich) <- c("1", "u1", "u2", "z1", "z2",
                          "u2sq", "z2sq", "u1u2")

    # Run GMM with enriched h_x
    res <- run_gmm(dat, dm, h_rich, lt)

    # Compute SE with diagnostics
    se_res <- compute_se_debug(dat, res$estimates, dm, h_rich, lt, n_val)

    tag <- if(is_degen) "DEGEN" else "CTRL"
    cat(sprintf("[%s] Rep %3d: cond=%.1e, alpha=(%.2f,%.2f,%.2f,%.2f), pi=[%.4f,%.4f], mu=%.4f, se=%s\n",
        tag, rep_i, res$cond,
        res$estimates[1], res$estimates[2], res$estimates[3], res$estimates[4],
        min(res$pi_hat), max(res$pi_hat),
        se_res$mu,
        if(is.na(se_res$se)) paste0("NA(", se_res$diag$fail, ")") else sprintf("%.4f", se_res$se)))

    # Print diagnostics for first few
    if (rep_i %in% head(test_reps, 5)) {
      d <- se_res$diag
      cat(sprintf("        W_cond=%.1e, GtWG_cond=%.1e, H_cond=%.1e\n",
          d$W_cond, d$GtWG_cond, d$H_cond))
      cat(sprintf("        |GtWG|=%.1e, |corrR|=%.1e, |corrS|=%.1e, |G|=%.1e\n",
          d$norm_GtWG, d$norm_corrR, d$norm_corrS, d$norm_G))
      if (!all(is.na(d$H_ev))) {
        cat(sprintf("        H eigenvalues: %s\n",
            paste(format(d$H_ev, digits=3, scientific=TRUE), collapse=", ")))
      }
    }
  }, error = function(e) {
    cat(sprintf("[%s] Rep %3d: ERROR: %s\n",
        if(is_degen) "DEGEN" else "CTRL", rep_i, e$message))
  })
}

## ============================================================
## STEP 3: If SE fails, try relaxed threshold or pseudoinverse
## ============================================================
cat("\n=== Step 3: Try fixes for SE computation ===\n\n")

# Pick first degenerate rep for detailed diagnosis
if (length(degen_reps) > 0) {
  rep_i <- degen_reps[1]
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]

  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[m_idx]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[m_idx]]),
    inv_link = ps_spec[["inv_link"]], outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL), optimizer = "L-BFGS-B"
  )

  ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
  dm <- ebmr$ps_fit.list[[1]]$design_matrix
  lt <- ebmr$ps_fit.list[[1]]$link_type

  h_rich <- cbind(1, dat$u1, dat$u2, dat$z1, dat$z2,
                  dat$u2^2, dat$z2^2, dat$u1 * dat$u2)

  res <- run_gmm(dat, dm, h_rich, lt)
  cat(sprintf("Rep %d with enriched h_x: cond=%.1e, still_degen=%s\n",
      rep_i, res$cond, res$cond >= cond_limit))
  cat(sprintf("  alpha = (%s)\n", paste(round(res$estimates, 3), collapse=", ")))
  cat(sprintf("  pi range = [%.6f, %.6f]\n", min(res$pi_hat), max(res$pi_hat)))

  # Try SE with different H_cond thresholds
  for (thresh in c(1e15, 1e20, 1e25, Inf)) {
    r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
    pi_hat <- res$pi_hat
    mu_hat <- mean(r_vec * y_vec / pi_hat)
    g_mat <- (r_vec / pi_hat - 1) * h_rich
    n <- n_val; k <- ncol(dm); h <- ncol(h_rich)

    W_hat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) NULL)
    if (is.null(W_hat)) { cat(sprintf("  thresh=%.0e: W singular\n", thresh)); next }

    cf <- r_vec * (1 - pi_hat) / pi_hat
    Gamma_hat <- crossprod(h_rich * cf, dm) / n
    GtWG <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
    GtW <- crossprod(Gamma_hat, W_hat)
    G_vec <- colMeans(g_mat)
    dc_X <- cf * dm
    R_mat <- matrix(0, h * k, k)
    for (j in 1:k) {
      block <- crossprod(dc_X, dm[, j] * h_rich) / n
      for (l in 1:h) R_mat[(l-1)*k + j, ] <- block[, l]
    }
    eta_s <- matrix(G_vec, h, 1)
    WG <- as.vector(W_hat %*% eta_s)
    GW_row <- as.vector(t(eta_s) %*% W_hat)
    S_mat <- matrix(0, h^2, k)
    for (j in 1:k) {
      dg_j <- cf * dm[, j] * h_rich
      cross <- (crossprod(dg_j, g_mat) + crossprod(g_mat, dg_j)) / n
      S_mat[, j] <- as.vector(cross)
    }
    GW_kron_Ik <- kronecker(matrix(GW_row, 1, h), diag(k))
    GW_kron_GtW <- kronecker(matrix(GW_row, 1, h), GtW)
    H_mat <- GtWG + GW_kron_Ik %*% R_mat - GW_kron_GtW %*% S_mat

    H_ev <- tryCatch(eigen(H_mat, symmetric = FALSE, only.values = TRUE)$values,
                     error = function(e) rep(NA, k))
    H_cond <- max(Mod(H_ev)) / max(min(Mod(H_ev)), 1e-15)

    if (H_cond > thresh) {
      cat(sprintf("  thresh=%.0e: H_cond=%.1e > threshold, skip\n", thresh, H_cond))
      next
    }

    H_inv <- tryCatch(solve(H_mat), error = function(e) NULL)
    if (is.null(H_inv)) {
      # Try pseudoinverse
      svd_H <- svd(H_mat)
      tol_sv <- max(dim(H_mat)) * max(svd_H$d) * .Machine$double.eps
      d_inv <- ifelse(svd_H$d > tol_sv, 1/svd_H$d, 0)
      H_inv <- svd_H$v %*% diag(d_inv) %*% t(svd_H$u)
      cat(sprintf("  thresh=%.0e: H_cond=%.1e, used pseudoinverse\n", thresh, H_cond))
    }

    GtW_g <- GtW %*% t(g_mat)
    g_WG <- as.vector(g_mat %*% WG)
    cf_hxWG <- cf * as.vector(h_rich %*% WG)
    GammaT_WG <- t(dm * cf_hxWG)
    Q_mat <- GtW_g + GammaT_WG - sweep(GtW_g, 2, g_WG, `*`)
    psi <- -(H_inv %*% Q_mat)

    dot_pi <- -dm * (pi_hat * (1 - pi_hat))
    ry_ps_inv2 <- as.vector(r_vec * y_vec * (pi_hat^(-2)))
    H_alpha <- colMeans(dot_pi * ry_ps_inv2)
    mu_iid <- as.vector(t(r_vec/pi_hat*y_vec) - t(H_alpha) %*% psi)
    se <- sqrt(var(mu_iid) / n)

    cat(sprintf("  thresh=%.0e: H_cond=%.1e, SE=%.4f, mu=%.4f\n",
        thresh, H_cond, se, mu_hat))
  }
}

cat("\nDone.\n")
