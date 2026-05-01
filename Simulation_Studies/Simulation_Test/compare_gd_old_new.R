setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
source("Data_Generation.r")
source("config/scenarios.R")

old_wd <- getwd()
setwd("../EBMRalgorithmFast4")
source("R/EBMRAlgorithm.r")
setwd(old_wd)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
inv_link_fn <- function(eta) 1 / (1 + exp(eta))

ps_spec_m2 <- list(
  formula.list = list(FORMULAS$u1_z2),
  h_alpha.list = list(H_ALPHA$full),
  inv_link = inv_link_fn,
  outcome = "y",
  alpha_init.list = list(NULL)
)

h_nu_func <- function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1,
                                  z2 = dat$z2, u1_u2 = dat$u1 * dat$u2)

# True mean
set.seed(999)
big <- setting3.A1(1e6)
mu_true <- mean(big$y)
rm(big); gc()
cat(sprintf("mu_true = %.6f\n\n", mu_true))

# Generate datasets (misspecified: setting3.B1)
set.seed(42)
n_reps <- 200
n_val <- 2000
datasets <- vector("list", n_reps)
for (i in 1:n_reps) datasets[[i]] <- setting3.B1(n_val)

# ============================================================================
# GD (new): L-BFGS-B with frequent W updates, bounds = Inf
# ============================================================================
cat("\n========== GD (new): L-BFGS-B inner steps, bounds=Inf ==========\n")
mu_new <- numeric(n_reps); se_new <- numeric(n_reps)
alpha_new_mat <- vector("list", n_reps)
conv_new <- logical(n_reps); iter_new <- integer(n_reps)
obj_new <- numeric(n_reps); gnorm_new <- numeric(n_reps)
min_ps_new <- numeric(n_reps); max_ipw_new <- numeric(n_reps)

t_start <- proc.time()
for (i in 1:n_reps) {
  dat <- datasets[[i]]
  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", ps_spec_m2, dat, W_func, bounds = Inf, gmm_method = "gd")
    ipw_result <- ebmr$EBMR_IPW(h_nu = h_nu_func)
    fit <- ebmr$ps_fit.list[[1]]
    ps_vals <- fit$fitted.values
    gmm_fit <- fit$gmm_fit

    mu_new[i] <- ipw_result$mu_ipw
    se_new[i] <- ipw_result$se_ipw
    alpha_new_mat[[i]] <- fit$coefficients
    conv_new[i] <- gmm_fit$opt$converged
    iter_new[i] <- gmm_fit$opt$iterations
    obj_new[i] <- gmm_fit$opt$objective
    gnorm_new[i] <- ifelse(is.null(gmm_fit$opt$final_grad_norm), NA, gmm_fit$opt$final_grad_norm)
    min_ps_new[i] <- min(ps_vals)
    max_ipw_new[i] <- max(dat$r / ps_vals)
  }, error = function(e) {
    mu_new[i] <<- NA; se_new[i] <<- NA
    alpha_new_mat[[i]] <<- rep(NA, 4); conv_new[i] <<- NA
    cat(sprintf("Rep %d: ERROR - %s\n", i, conditionMessage(e)))
  })
  if (i %% 50 == 0) cat(sprintf("  %d/%d done\n", i, n_reps))
}
elapsed_new <- (proc.time() - t_start)[3]
cat(sprintf("  Elapsed: %.1f sec\n", elapsed_new))

# ============================================================================
# GD (old): Pure steepest descent with Armijo line search, bounds = Inf
# Uses a custom gmm_old function embedded here
# ============================================================================
cat("\n========== GD (old): Pure steepest descent, bounds=Inf ==========\n")

# Standalone old GD GMM implementation
gmm_old_gd <- function(g, W, n, esteq_dim, param_dim, init, se.fit = TRUE, dg = NULL, d2g = NULL) {

  G = function(param) {
    g.matrix = g(param)
    return(matrix(.colMeans(g.matrix, n, esteq_dim), esteq_dim, 1))
  }

  t_Gamma_i = function(param) {
    if (!is.null(dg)) return(dg(param))
    Gamma.arr = array(NA, dim = c(n, param_dim, esteq_dim))
    for(l in 1:esteq_dim) {
      Gamma.arr[,,l] = numDeriv::jacobian(function(param) g(param)[, l], param)
    }
    return(Gamma.arr)
  }

  Gamma = function(param, t_Gamma_arr = NULL) {
    if (is.null(t_Gamma_arr)) t_Gamma_arr = t_Gamma_i(param)
    Gamma_mat = matrix(0, esteq_dim, param_dim)
    for (j in 1:param_dim) Gamma_mat[, j] = colMeans(t_Gamma_arr[, j, ])
    return(Gamma_mat)
  }

  Gamma_2 = function(param) {
    if (!is.null(d2g)) return(d2g(param))
    Gamma_from_G = function(p) numDeriv::jacobian(function(pp) as.vector(G(pp)), p)
    return(numDeriv::jacobian(function(p) as.vector(t(Gamma_from_G(p))), param))
  }

  if(is.null(init)) init = rep(0, param_dim)

  # Step 1: Initialize with W=I via L-BFGS-B (same as new GD)
  obj_init <- function(param) {
    g_mat <- g(param)
    G_hat <- matrix(.colMeans(g_mat, n, esteq_dim), esteq_dim, 1)
    as.numeric(crossprod(G_hat))
  }
  grad_init <- if (!is.null(dg)) {
    function(param) {
      g_mat <- g(param)
      G_hat <- matrix(.colMeans(g_mat, n, esteq_dim), esteq_dim, 1)
      dg_arr <- dg(param)
      Gamma_mat <- matrix(0, esteq_dim, param_dim)
      for (j in 1:param_dim) Gamma_mat[, j] <- .colMeans(dg_arr[, j, , drop = FALSE], n, esteq_dim)
      2 * as.vector(crossprod(Gamma_mat, G_hat))
    }
  } else { NULL }

  opt_init <- optim(init, obj_init, gr = grad_init, method = "L-BFGS-B",
                    control = list(maxit = 1000))
  estimates <- opt_init$par

  # Step 2: Pure steepest descent with Armijo backtracking
  max_gd <- 2000L
  gd_tol <- 1e-6
  final_grad_norm <- NA

  for (t in 1:max_gd) {
    g_mat <- g(estimates)
    G_hat <- matrix(.colMeans(g_mat, n, esteq_dim), esteq_dim, 1)
    W_t <- tryCatch(W(g_mat), error = function(e) diag(esteq_dim))
    Q_current <- as.numeric(crossprod(G_hat, W_t %*% G_hat))

    # Frozen-W gradient: 2 * Gamma' * W_t * G
    if (!is.null(dg)) {
      dg_arr <- dg(estimates)
    } else {
      dg_arr <- t_Gamma_i(estimates)
    }
    Gamma_mat <- matrix(0, esteq_dim, param_dim)
    for (j in 1:param_dim) {
      Gamma_mat[, j] <- .colMeans(dg_arr[, j, , drop = FALSE], n, esteq_dim)
    }
    grad_val <- 2 * as.vector(crossprod(Gamma_mat, W_t %*% G_hat))
    grad_norm <- max(abs(grad_val))
    final_grad_norm <- grad_norm
    if (grad_norm < gd_tol) break

    # Armijo backtracking line search
    step <- 1.0; c1 <- 1e-4; direction <- -grad_val
    slope <- sum(grad_val * direction)
    for (ls in 1:30) {
      alpha_try <- estimates + step * direction
      g_try <- g(alpha_try)
      G_try <- matrix(.colMeans(g_try, n, esteq_dim), esteq_dim, 1)
      Q_try <- as.numeric(crossprod(G_try, W_t %*% G_try))
      if (is.finite(Q_try) && Q_try <= Q_current + c1 * step * slope) break
      step <- step * 0.5
    }
    alpha_new_val <- estimates + step * direction
    conv_err <- max(abs(alpha_new_val - estimates))
    estimates <- alpha_new_val
    if (conv_err < gd_tol) break
  }

  iter_total <- t

  # Compute objective at final estimates
  g_mat_final <- g(estimates)
  G_final <- matrix(.colMeans(g_mat_final, n, esteq_dim), esteq_dim, 1)
  W_final <- tryCatch(W(g_mat_final), error = function(e) diag(esteq_dim))
  objective <- as.numeric(crossprod(G_final, W_final %*% G_final))
  opt <- list(
    converged = !is.na(final_grad_norm) && final_grad_norm < 1e-4,
    objective = objective,
    iterations = iter_total,
    method = "gd_old",
    final_grad_norm = final_grad_norm
  )

  # SE computation (same as standard gmm)
  Gamma.hat = W.hat = g.matrix = eta_s = Q = h_x = psi = se = NA
  if(se.fit) {
    t_Gamma_arr = t_Gamma_i(estimates)
    Gamma.hat = Gamma(estimates, t_Gamma_arr)
    g.matrix = g(estimates)
    W.hat = W(g.matrix)
    t_g.matrix = t(g.matrix)
    eta_s = .colMeans(g.matrix, n, esteq_dim)
    M_n = kronecker(eta_s %*% W.hat, diag(param_dim)) %*% Gamma_2(estimates)
    GtW = crossprod(Gamma.hat, W.hat)
    GtWG = GtW %*% Gamma.hat
    H_0n = tryCatch(-solve(GtWG), error = function(e) -solve(GtWG + 1e-8 * diag(nrow(GtWG))))
    H_1n = GtW %*% (t_g.matrix - eta_s)
    W_eta = as.vector(W.hat %*% eta_s)
    tGamma_W_eta = as.vector(GtW %*% eta_s)
    dim_arr = dim(t_Gamma_arr)
    t_Gamma_mat = matrix(t_Gamma_arr, nrow = dim_arr[1] * dim_arr[2], ncol = dim_arr[3])
    H_2n_2_vec = t_Gamma_mat %*% W_eta
    H_2n_2 = matrix(H_2n_2_vec, nrow = n, ncol = param_dim, byrow = FALSE)
    H_2n_2 = t(H_2n_2) - tGamma_W_eta
    GtW_eta = tGamma_W_eta
    GtW_g = GtW %*% t_g.matrix
    g_W_eta = as.vector(g.matrix %*% W_eta)
    H_2n_3 = matrix(GtW_eta, param_dim, n) - GtW_g * matrix(g_W_eta, param_dim, n, byrow = TRUE)
    I_minus_H0M = diag(param_dim) - H_0n %*% M_n
    Q = tryCatch(solve(I_minus_H0M) %*% H_0n, error = function(e) solve(I_minus_H0M + 1e-8 * diag(param_dim)) %*% H_0n)
    psi = Q %*% (H_1n + H_2n_2 + H_2n_3)
    cov.hat = var(t(psi)) / n
    se = sqrt(diag(cov.hat))
  }

  return(list(estimates = estimates, se = se, Gamma.hat = Gamma.hat, W.hat = W.hat,
              g.matrix = g.matrix, eta_s = eta_s, Q = Q, h_x = h_x, psi = psi, opt = opt))
}

# We need to run the old GD through the EBMRAlgorithm framework.
# The simplest approach: temporarily patch the gmm function.
# Save original, replace, run, restore.

# Actually, let's just source the EBMRAlgorithm code with a modified gmm.
# We'll create a wrapper class that uses gmm_old_gd.

# Simplest: monkey-patch the private gmm method
# Store the original
orig_gmm <- EBMRAlgorithmFast4$private_methods$gmm

# Patch with old GD
EBMRAlgorithmFast4$set("private", "gmm_patched", gmm_old_gd, overwrite = TRUE)

# We can't easily swap private methods in R6, so let's use a different approach:
# Create a modified copy of the GMM function file and source it.

# Actually the cleanest way: just temporarily replace the gmm in the class.
# R6 allows overwriting private methods with $set if overwrite=TRUE

# Let's wrap gmm_old_gd to match the expected signature
gmm_old_wrapper <- function(g, W, n, esteq_dim, param_dim, init, se.fit = TRUE,
                            dg = NULL, d2g = NULL, bounds = 5, method = c("iterated", "gd")) {
  gmm_old_gd(g, W, n, esteq_dim, param_dim, init, se.fit, dg, d2g)
}

# Replace the private gmm method
EBMRAlgorithmFast4$set("private", "gmm", gmm_old_wrapper, overwrite = TRUE)

mu_old <- numeric(n_reps); se_old <- numeric(n_reps)
alpha_old_mat <- vector("list", n_reps)
conv_old <- logical(n_reps); iter_old <- integer(n_reps)
obj_old <- numeric(n_reps); gnorm_old <- numeric(n_reps)
min_ps_old <- numeric(n_reps); max_ipw_old <- numeric(n_reps)

t_start <- proc.time()
for (i in 1:n_reps) {
  dat <- datasets[[i]]
  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", ps_spec_m2, dat, W_func, bounds = Inf, gmm_method = "gd")
    ipw_result <- ebmr$EBMR_IPW(h_nu = h_nu_func)
    fit <- ebmr$ps_fit.list[[1]]
    ps_vals <- fit$fitted.values
    gmm_fit <- fit$gmm_fit

    mu_old[i] <- ipw_result$mu_ipw
    se_old[i] <- ipw_result$se_ipw
    alpha_old_mat[[i]] <- fit$coefficients
    conv_old[i] <- gmm_fit$opt$converged
    iter_old[i] <- gmm_fit$opt$iterations
    obj_old[i] <- gmm_fit$opt$objective
    gnorm_old[i] <- ifelse(is.null(gmm_fit$opt$final_grad_norm), NA, gmm_fit$opt$final_grad_norm)
    min_ps_old[i] <- min(ps_vals)
    max_ipw_old[i] <- max(dat$r / ps_vals)
  }, error = function(e) {
    mu_old[i] <<- NA; se_old[i] <<- NA
    alpha_old_mat[[i]] <<- rep(NA, 4); conv_old[i] <<- NA
    cat(sprintf("Rep %d: ERROR - %s\n", i, conditionMessage(e)))
  })
  if (i %% 50 == 0) cat(sprintf("  %d/%d done\n", i, n_reps))
}
elapsed_old <- (proc.time() - t_start)[3]
cat(sprintf("  Elapsed: %.1f sec\n", elapsed_old))

# Restore original gmm
EBMRAlgorithmFast4$set("private", "gmm", orig_gmm, overwrite = TRUE)

# ============================================================================
# Report
# ============================================================================
report <- function(label, mu_hat, se_hat, alpha_mat, converged, iterations,
                   objective, grad_norm_vec, min_ps, max_ipw, elapsed) {
  cat(sprintf("\n\n########## %s (elapsed: %.1f sec) ##########\n", label, elapsed))
  valid <- !is.na(mu_hat) & !is.na(se_hat)
  cat(sprintf("Valid reps: %d/%d\n", sum(valid), n_reps))
  if (sum(valid) == 0) return()

  mu_v <- mu_hat[valid]; se_v <- se_hat[valid]
  q1 <- quantile(mu_v, 0.25); q3 <- quantile(mu_v, 0.75); iqr <- q3 - q1
  is_outlier <- mu_v < (q1 - 3*iqr) | mu_v > (q3 + 3*iqr)
  cat(sprintf("Outliers removed: %d/%d (%.1f%%)\n", sum(is_outlier), sum(valid), 100*mean(is_outlier)))

  mu_v <- mu_v[!is_outlier]; se_v <- se_v[!is_outlier]; n_kept <- length(mu_v)
  empirical_sd <- sd(mu_v)
  cat(sprintf("Bias: %.4f, RMSE: %.4f, SD: %.4f\n",
              mean(mu_v) - mu_true, sqrt(mean((mu_v - mu_true)^2)), empirical_sd))
  cat(sprintf("  mean(SE): %.4f, mean(SE)/SD: %.4f\n", mean(se_v), mean(se_v)/empirical_sd))
  cat(sprintf("  median(SE): %.4f, median(SE)/SD: %.4f\n", median(se_v), median(se_v)/empirical_sd))

  ci_lower <- mu_v - 1.96*se_v; ci_upper <- mu_v + 1.96*se_v
  covers <- (ci_lower <= mu_true) & (mu_true <= ci_upper)
  cat(sprintf("  95%% CI coverage: %d/%d (%.1f%%)\n", sum(covers), n_kept, 100*mean(covers)))

  valid_idx <- which(valid); kept_idx <- valid_idx[!is_outlier]
  cat(sprintf("  Convergence: %d/%d (%.0f%%)\n",
              sum(converged[kept_idx]), n_kept, 100*mean(converged[kept_idx])))
  cat(sprintf("  Iterations: mean=%.1f, range=[%d, %d]\n",
              mean(iterations[kept_idx]), min(iterations[kept_idx]), max(iterations[kept_idx])))
  fgn <- grad_norm_vec[kept_idx]; fgn <- fgn[!is.na(fgn)]
  if (length(fgn) > 0) cat(sprintf("  Final grad norm: mean=%.6f, median=%.6f, max=%.6f\n",
                                     mean(fgn), median(fgn), max(fgn)))

  alpha_all <- do.call(rbind, alpha_mat[kept_idx])
  anames <- names(alpha_mat[[kept_idx[1]]])
  if (is.null(anames)) anames <- paste0("alpha", 1:ncol(alpha_all))
  cat(sprintf("  Alpha coefficients:\n"))
  for (k in 1:ncol(alpha_all)) {
    cat(sprintf("    %s: mean=%7.3f, sd=%.3f, range=[%.3f, %.3f]\n",
                anames[k], mean(alpha_all[,k]), sd(alpha_all[,k]), min(alpha_all[,k]), max(alpha_all[,k])))
  }
  cat(sprintf("  min(PS): mean=%.6f, min=%.6f\n", mean(min_ps[kept_idx]), min(min_ps[kept_idx])))
  cat(sprintf("  max(r/PS): mean=%.1f, max=%.1f\n", mean(max_ipw[kept_idx]), max(max_ipw[kept_idx])))
}

report("GD (new): L-BFGS-B inner steps, bounds=Inf", mu_new, se_new, alpha_new_mat,
       conv_new, iter_new, obj_new, gnorm_new, min_ps_new, max_ipw_new, elapsed_new)

report("GD (old): Pure steepest descent, bounds=Inf", mu_old, se_old, alpha_old_mat,
       conv_old, iter_old, obj_old, gnorm_old, min_ps_old, max_ipw_old, elapsed_old)

cat("\nDone!\n")
