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

model_idx <- 3

# First identify all 17 degenerate reps
degen_reps <- c()
for (rep_i in 1:200) {
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
    sc <- ebmr$ps_fit.list[[1]]$gmm_fit$opt$solution_cond
    if (!is.na(sc) && sc >= 1e8) degen_reps <- c(degen_reps, rep_i)
  }, error = function(e) { degen_reps <<- c(degen_reps, rep_i) })
}
cat(sprintf("Degenerate reps (%d): %s\n\n", length(degen_reps), paste(degen_reps, collapse=", ")))

# Now manually run constrained_nr for each, tracing grad/cond at every outer iteration
for (rep_i in degen_reps) {
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
  r_vec <- dat[["r"]]
  n <- n_val

  # Build components from package
  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL),
    optimizer = "L-BFGS-B"  # just to build the components
  )
  ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
  design_mat <- ebmr$ps_fit.list[[1]]$design_matrix
  h_x <- ebmr$ps_fit.list[[1]]$h_x

  compute_pi <- function(eta) plogis(-eta)
  k <- ncol(design_mat)
  h <- ncol(h_x)

  g_fn <- function(alpha) {
    pi_v <- compute_pi(as.vector(design_mat %*% alpha))
    (r_vec / pi_v - 1) * h_x
  }
  Gamma_fn <- function(alpha) {
    pi_v <- compute_pi(as.vector(design_mat %*% alpha))
    cf <- r_vec * (1 - pi_v) / pi_v
    crossprod(h_x * cf, design_mat) / n
  }

  # Manual constrained NR with full tracing
  cond_threshold <- 1e8
  trust_radius <- 2.0
  max_inner <- 50
  max_outer <- 200
  tol <- 1e-6
  stall_limit <- 5

  nr_inner <- function(alpha, W_hat) {
    for (iter in 1:max_inner) {
      g_mat <- g_fn(alpha)
      G_vec <- matrix(colMeans(g_mat), h, 1)
      Gamma_hat <- Gamma_fn(alpha)
      g_vec <- 2 * as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))
      if (max(abs(g_vec)) < tol) break
      GtWG <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
      H_mat <- GtWG + GtWG
      direction <- tryCatch(-solve(H_mat, g_vec),
                            error = function(e) -g_vec * trust_radius / max(abs(g_vec)))
      if (sqrt(sum(direction^2)) > trust_radius)
        direction <- direction * trust_radius / sqrt(sum(direction^2))
      obj_cur <- as.numeric(crossprod(G_vec, W_hat %*% G_vec))
      accepted <- FALSE
      step <- 1.0
      for (ls in 1:20) {
        cand <- alpha + step * direction
        g_cand <- g_fn(cand)
        G_cand <- matrix(colMeans(g_cand), h, 1)
        obj_new <- as.numeric(crossprod(G_cand, W_hat %*% G_cand))
        if (is.finite(obj_new) && obj_new < obj_cur - 1e-4 * step * sum(g_vec * direction)) {
          Gamma_c <- Gamma_fn(cand)
          H_c <- crossprod(Gamma_c, W_hat %*% Gamma_c)
          eig_c <- eigen(H_c, symmetric = TRUE, only.values = TRUE)$values
          cc <- max(eig_c) / max(min(eig_c), 1e-15)
          if (cc < cond_threshold) { accepted <- TRUE; break }
        }
        step <- step * 0.5
      }
      if (accepted) { alpha <- cand } else { break }
    }
    alpha
  }

  # Step 1: W = I
  estimates <- nr_inner(rep(0, k), diag(h))

  # Step 2: Iterative with full trace
  best_grad_norm <- Inf
  best_estimates <- estimates
  best_cond <- NA
  no_improve_count <- 0L
  n_restarts <- 0L
  max_restarts <- 5L

  cat(sprintf("=== Rep %d ===\n", rep_i))
  cat(sprintf("  %4s  %12s  %12s  %12s  %s\n", "iter", "grad", "cond", "obj", "alpha[2]"))

  for (t in 1:max_outer) {
    g_mat <- g_fn(estimates)
    G_vec <- matrix(colMeans(g_mat), h, 1)
    W_hat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(h))
    Gamma_hat <- Gamma_fn(estimates)
    cg <- 2 * as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))
    grad_norm <- max(abs(cg))
    H_cur <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
    eig_cur <- eigen(H_cur, symmetric = TRUE, only.values = TRUE)$values
    cond_cur <- max(eig_cur) / max(min(eig_cur), 1e-15)
    obj_cur <- as.numeric(crossprod(G_vec, W_hat %*% G_vec))

    # Log every iteration for first 20, then every 10th
    if (t <= 20 || t %% 10 == 0 || grad_norm < tol) {
      cat(sprintf("  %4d  %12.4e  %12.2e  %12.6f  %.4f\n",
          t, grad_norm, cond_cur, obj_cur, estimates[2]))
    }

    if (grad_norm < tol) {
      cat(sprintf("  ** CONVERGED at iter %d **\n", t))
      break
    }

    if (grad_norm < best_grad_norm) {
      best_grad_norm <- grad_norm
      best_estimates <- estimates
      best_cond <- cond_cur
      no_improve_count <- 0L
    } else {
      no_improve_count <- no_improve_count + 1L
      if (no_improve_count >= stall_limit) {
        if (n_restarts < max_restarts) {
          n_restarts <- n_restarts + 1L
          cat(sprintf("  ** STALL at iter %d, restart #%d (best_grad=%.2e, best_cond=%.2e) **\n",
              t, n_restarts, best_grad_norm, best_cond))
          if (n_restarts <= 2L) {
            scale <- 0.1 * n_restarts
            perturb <- rnorm(k) * scale * (abs(best_estimates) + 0.1)
            restart_init <- best_estimates + perturb
          } else if (n_restarts == 3L) {
            g_b <- g_fn(best_estimates)
            G_b <- matrix(colMeans(g_b), h, 1)
            W_b <- tryCatch(solve(crossprod(g_b) / n), error = function(e) diag(h))
            grad_b <- 2 * as.vector(crossprod(Gamma_fn(best_estimates), W_b %*% G_b))
            restart_init <- best_estimates - trust_radius * grad_b / max(sqrt(sum(grad_b^2)), 1e-10)
          } else {
            restart_init <- rnorm(k) * 0.5
          }
          estimates <- nr_inner(restart_init, W_hat)
          g_rc <- g_fn(estimates)
          G_rc <- matrix(colMeans(g_rc), h, 1)
          W_rc <- tryCatch(solve(crossprod(g_rc) / n), error = function(e) diag(h))
          rg <- max(abs(2 * as.vector(crossprod(Gamma_fn(estimates), W_rc %*% G_rc))))
          if (rg < best_grad_norm) {
            best_grad_norm <- rg
            best_estimates <- estimates
          }
          no_improve_count <- 0L
          next
        } else {
          cat(sprintf("  ** EXHAUSTED restarts at iter %d (best_grad=%.2e, best_cond=%.2e) **\n",
              t, best_grad_norm, best_cond))
          estimates <- best_estimates
          break
        }
      }
    }
    estimates <- nr_inner(estimates, W_hat)
  }

  cat(sprintf("  Final: grad=%.2e, cond=%.2e, alpha=%s\n\n",
      best_grad_norm, best_cond, paste(round(best_estimates, 4), collapse=", ")))
}
