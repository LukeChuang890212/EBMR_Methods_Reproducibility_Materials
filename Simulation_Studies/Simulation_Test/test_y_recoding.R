## Test: Recode binary y from {0,1} to {-1,1} in design matrix
## Rationale: With y in {0,1}, extreme alpha_y only affects y=1 obs (alpha_y*0=0 for y=0),
## allowing degeneracy where pi->1 for y=1 group without penalty.
## With y in {-1,1}, extreme alpha_y pushes pi in opposite directions for BOTH groups,
## so pi->1 for one group forces pi->0 for the other, making weights explode.
##
## We test: S4 M3 (binary, misspecified overid - main problem case),
##          S4 M2 (binary, misspecified overid),
##          S3 M2 (continuous, misspecified overid - should be unaffected)
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)
library(parallel)

n_val <- 2000; n_reps <- 1000
ps_spec <- get_ps_spec("9-alt1")
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
compute_pi_fn <- function(eta) plogis(-eta)

## Use the same constrained_nr optimizer from the package
## (standard approach, just with recoded design matrix)

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

  # Iterative GMM: W=I start, then iterate
  estimates <- nr_inner(init, diag(h_dim))

  for (t in 1:max_outer) {
    g_mat <- g_fn(estimates)
    G_vec <- matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
    W_hat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(h_dim))
    Gamma <- Gamma_fn(estimates)
    cg <- 2 * as.vector(crossprod(Gamma, W_hat %*% G_vec))
    if (max(abs(cg)) < tol) break

    # One Gauss-Newton step per outer iteration
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

  # Final diagnostics
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

## Recode binary y columns in design matrix from {0,1} to {-1,1}
## n_y = number of y-related columns (starts at column 2 after intercept)
recode_y_in_design <- function(design_mat, n_y, y_binary) {
  if (n_y == 0 || !y_binary) return(design_mat)
  dm <- design_mat
  # The y columns are columns 2 through (1 + n_y)
  for (j in 2:(1 + n_y)) {
    # For main y column: y -> 2*y - 1
    # For interaction y*x column: y*x -> (2*y-1)*x
    # Both cases: just replace wherever y=0 with -1, y=1 stays 1
    # Since these columns are y * (something), where y in {0,1}:
    #   when y=0: column value = 0 -> should become -1 * (something)
    #   when y=1: column value = 1 * (something) -> stays same
    # But we need to know which obs have y=0 vs y=1
    # Simpler: for pure y column, 2*col - 1. For interaction y*x: 2*col - x
    # Actually simplest: replace column with (2*y_01 - 1) * (original_x_part)
    # For the main y column (col value = y): new = 2*y - 1
    # For interaction y*z: col value = y*z. new = (2*y-1)*z
    # We can detect: for y=1 obs, the column value = x_part. For y=0 obs, value = 0.
    # So: x_part[y=1] = col[y=1], x_part[y=0] = need to figure out...
    # Actually this is tricky for interactions. Let's just handle the case we know:
    # In our settings, n_y = 1 and the y column is just y itself (no interactions).
    NULL
  }
  # For our specific case: column 2 is pure y (no interactions in u1_z2, u2_z2 formulas)
  # Recode: 0 -> -1, 1 -> 1
  dm[, 2] <- 2 * design_mat[, 2] - 1
  return(dm)
}

settings <- list(
  list(name = "S4 M3", setting = "setting4", m_idx = 3, binary_y = TRUE),
  list(name = "S4 M2", setting = "setting4", m_idx = 2, binary_y = TRUE),
  list(name = "S3 M2", setting = "setting3", m_idx = 2, binary_y = FALSE)
)

cond_limit <- 1e8

cat("=== Test: y recoding {0,1} -> {-1,1} in design matrix ===\n\n")

for (cfg in settings) {
  cat(sprintf("========== %s ==========\n", cfg$name))
  data_file <- misspecified_model_all_data_file.list[[cfg$setting]][["miss50"]][[1]]
  all_data <- readRDS(data_file)
  mu_true <- get_mu_true(cfg$setting)
  cat(sprintf("True mean: %.6f\n", mu_true))

  cat("Pre-building design matrices...\n")
  design_info <- vector("list", n_reps)
  for (rep_i in 1:n_reps) {
    dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
    single_ps <- list(
      formula.list = list(ps_spec[["formula.list"]][[cfg$m_idx]]),
      h_alpha.list = list(ps_spec[["h_alpha.list"]][[cfg$m_idx]]),
      inv_link = ps_spec[["inv_link"]], outcome = ps_spec[["outcome"]],
      alpha_init.list = list(NULL), optimizer = "L-BFGS-B"
    )
    tryCatch({
      ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
      dm <- ebmr$ps_fit.list[[1]]$design_matrix
      hx <- ebmr$ps_fit.list[[1]]$h_x
      lt <- ebmr$ps_fit.list[[1]]$link_type
      n_y <- ebmr$ps_fit.list[[1]]$n_y
      # Apply recoding for binary y
      if (cfg$binary_y) {
        dm[, 2] <- 2 * dm[, 2] - 1
      }
      design_info[[rep_i]] <- list(design_mat = dm, h_x = hx, link_type = lt, n_y = n_y)
    }, error = function(e) NULL)
  }
  cat("Done.\n")

  # Show design matrix info for first rep
  if (!is.null(design_info[[1]])) {
    cat(sprintf("  Design matrix cols: %s\n", paste(colnames(design_info[[1]]$design_mat), collapse=", ")))
    cat(sprintf("  k=%d, h_dim=%d, n_y=%d\n",
        ncol(design_info[[1]]$design_mat), ncol(design_info[[1]]$h_x),
        design_info[[1]]$n_y))
    if (cfg$binary_y) {
      cat(sprintf("  y column (recoded) range: [%g, %g]\n",
          min(design_info[[1]]$design_mat[, 2]), max(design_info[[1]]$design_mat[, 2])))
    }
  }

  # Run in parallel
  n_cores <- min(detectCores() - 1, 10)
  cl <- makeCluster(n_cores)
  clusterExport(cl, c("run_gmm", "compute_se", "compute_pi_fn", "W_fn",
                       "all_data", "design_info", "n_val"),
                envir = environment())

  results <- parLapply(cl, 1:n_reps, function(rep_i) {
    if (is.null(design_info[[rep_i]])) return(list(mu = NA, se = NA, cond = NA))
    dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
    di <- design_info[[rep_i]]
    tryCatch({
      res <- run_gmm(dat, di$design_mat, di$h_x, di$link_type)
      se_res <- compute_se(dat, res$estimates, di$design_mat, di$h_x, di$link_type, n_val)
      list(mu = se_res$mu, se = se_res$se, cond = res$cond,
           alpha = res$estimates, pi_min = min(res$pi_hat), pi_max = max(res$pi_hat))
    }, error = function(e) list(mu = NA, se = NA, cond = NA))
  })
  stopCluster(cl)

  mu_vals <- sapply(results, `[[`, "mu")
  se_vals <- sapply(results, `[[`, "se")
  cond_vals <- sapply(results, `[[`, "cond")

  valid <- !is.na(mu_vals) & !is.na(se_vals)
  n_degen <- sum(cond_vals[valid] >= cond_limit, na.rm = TRUE)

  # Check pi distributions for degenerate reps
  if (cfg$binary_y) {
    pi_mins <- sapply(results, function(r) if(!is.null(r$pi_min)) r$pi_min else NA)
    pi_maxs <- sapply(results, function(r) if(!is.null(r$pi_max)) r$pi_max else NA)
    degen_mask <- valid & (cond_vals >= cond_limit)
    if (any(degen_mask, na.rm = TRUE)) {
      cat(sprintf("  Degenerate reps pi range: [%.6f, %.6f]\n",
          min(pi_mins[degen_mask], na.rm=TRUE), max(pi_maxs[degen_mask], na.rm=TRUE)))
    }
  }

  # Check for extreme alpha_y in degenerate reps
  if (cfg$binary_y && any(valid & (cond_vals >= cond_limit), na.rm = TRUE)) {
    degen_idx <- which(valid & (cond_vals >= cond_limit))
    alpha_y_vals <- sapply(results[degen_idx], function(r) r$alpha[2])
    cat(sprintf("  Degenerate alpha_y: median=%.2f, range=[%.2f, %.2f]\n",
        median(alpha_y_vals, na.rm=TRUE),
        min(alpha_y_vals, na.rm=TRUE), max(alpha_y_vals, na.rm=TRUE)))
  }

  cat(sprintf("  Valid: %d/%d, Degen: %d\n", sum(valid), n_reps, n_degen))
  cat(sprintf("  Bias     = %.4f\n", mean(mu_vals[valid]) - mu_true))
  cat(sprintf("  ESD      = %.4f\n", sd(mu_vals[valid])))
  cat(sprintf("  ESE      = %.4f\n", mean(se_vals[valid])))
  cat(sprintf("  ESE/ESD  = %.3f\n", mean(se_vals[valid]) / sd(mu_vals[valid])))
  ci_lo <- mu_vals[valid] - 1.96 * se_vals[valid]
  ci_hi <- mu_vals[valid] + 1.96 * se_vals[valid]
  cat(sprintf("  CP       = %.3f\n", mean(mu_true >= ci_lo & mu_true <= ci_hi)))
  cat("\n")
}
