## Test: Replicate GMM inline with controllable cond_threshold
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Mental_Health_Data_Application")
library(Matrix)
library(numDeriv)
library(tidyverse)
source("MHD_functions.R")

# Also load package for design matrix building only
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

original_data <- read.csv("data_application.csv")
percent <- original_data$Percentage
n_total <- 2486
class_n <- round(n_total * percent / 100)
dat <- gen_data(original_data, class_n, n_total)
dat$y <- dat$teacher_report
n <- nrow(dat)

W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
compute_pi_fn <- function(eta) plogis(-eta)  # logistic_complement: 1/(1+exp(eta))

# Use the package just to build design matrices and h_x
ps_specifications <- list(
  formula.list = list(
    r ~ teacher_report + health + father + teacher_report:health + teacher_report:father,
    r ~ teacher_report + health + parent_report + teacher_report:health + teacher_report:parent_report,
    r ~ teacher_report + parent_report + father + teacher_report:parent_report + teacher_report:father
  ),
  h_alpha.list = list(
    c("health", "father", "parent_report", "fp", "fh", "hp"),
    c("health", "father", "parent_report", "fp", "fh", "hp"),
    c("health", "father", "parent_report", "fp", "fh", "hp")
  ),
  outcome = "teacher_report",
  inv_link = function(eta) 1 / (1 + exp(eta))
)

# Build design matrices via package
ebmr_ref <- EBMRAlgorithmFast4$new("teacher_report", ps_specifications, dat, W_fn)
model_info <- lapply(1:3, function(i) {
  pf <- ebmr_ref$ps_fit.list[[i]]
  list(design_mat = pf$design_matrix, h_x = pf$h_x)
})

## Inline GMM with cond_threshold control
run_gmm_cond <- function(design_mat, h_x, r_vec, cond_threshold, trust_radius = 2.0) {
  n <- length(r_vec)
  k <- ncol(design_mat)
  h_dim <- ncol(h_x)
  init <- rep(0, k)
  tol <- 1e-6
  max_outer <- 200L
  stall_limit <- 20L

  g_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    (r_vec / pi_v - 1) * h_x
  }
  Gamma_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    cf <- r_vec * (1 - pi_v) / pi_v
    crossprod(h_x * cf, design_mat) / n
  }

  forward_cond <- function(alpha, W_hat) {
    Gamma <- Gamma_fn(alpha)
    H <- crossprod(Gamma, W_hat %*% Gamma)
    eig <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    max(eig) / max(min(eig), 1e-15)
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
          # Degeneracy guard
          Gamma_c <- Gamma_fn(cand)
          H_c <- crossprod(Gamma_c, W_hat %*% Gamma_c)
          eig_c <- eigen(H_c, symmetric = TRUE, only.values = TRUE)$values
          cc <- max(eig_c) / max(min(eig_c), 1e-15)
          if (cc < cond_threshold) { accepted <- TRUE; break }
        }
        step <- step * 0.5
      }
      if (accepted) alpha <- cand else break
    }
    alpha
  }

  # Step 1: W = I
  estimates <- nr_inner(init, diag(h_dim))

  # Step 2: Iterative GMM
  best_grad_norm <- Inf
  best_estimates <- estimates
  no_improve_count <- 0L
  n_restarts <- 0L
  max_restarts <- 5L

  for (t in 1:max_outer) {
    g_mat <- g_fn(estimates)
    G_vec <- matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
    W_hat <- tryCatch(W_fn(g_mat), error = function(e) diag(h_dim))
    Gamma <- Gamma_fn(estimates)
    cg <- 2 * as.vector(crossprod(Gamma, W_hat %*% G_vec))
    grad_norm <- max(abs(cg))
    if (grad_norm < tol) break

    if (grad_norm < best_grad_norm) {
      best_grad_norm <- grad_norm
      best_estimates <- estimates
      no_improve_count <- 0L
    } else {
      no_improve_count <- no_improve_count + 1L
      if (no_improve_count >= stall_limit) {
        if (n_restarts < max_restarts) {
          n_restarts <- n_restarts + 1L
          if (n_restarts <= 2L) {
            scale <- 0.1 * n_restarts
            perturb <- rnorm(k) * scale * (abs(best_estimates) + 0.1)
            restart_init <- best_estimates + perturb
          } else if (n_restarts == 3L) {
            g_b <- g_fn(best_estimates)
            G_b <- matrix(.colMeans(g_b, n, h_dim), h_dim, 1)
            W_b <- tryCatch(W_fn(g_b), error = function(e) diag(h_dim))
            grad_b <- 2 * as.vector(crossprod(Gamma_fn(best_estimates), W_b %*% G_b))
            restart_init <- best_estimates - trust_radius * grad_b / max(sqrt(sum(grad_b^2)), 1e-10)
          } else {
            restart_init <- rnorm(k) * 0.5
          }
          estimates <- nr_inner(restart_init, W_hat)
          g_rc <- g_fn(estimates)
          G_rc <- matrix(.colMeans(g_rc, n, h_dim), h_dim, 1)
          W_rc <- tryCatch(W_fn(g_rc), error = function(e) diag(h_dim))
          rg <- max(abs(2 * as.vector(crossprod(Gamma_fn(estimates), W_rc %*% G_rc))))
          if (rg < best_grad_norm) {
            best_grad_norm <- rg
            best_estimates <- estimates
          }
          no_improve_count <- 0L
          next
        } else {
          estimates <- best_estimates
          break
        }
      }
    }

    estimates <- nr_inner(estimates, W_hat)
  }

  # Final
  pi_final <- compute_pi_fn(as.vector(design_mat %*% estimates))
  g_mat_f <- g_fn(estimates)
  W_f <- tryCatch(W_fn(g_mat_f), error = function(e) diag(h_dim))
  cond_final <- forward_cond(estimates, W_f)
  mu_ipw <- mean(r_vec * dat$teacher_report / pi_final)

  list(alpha = estimates, pi = pi_final, cond = cond_final, mu = mu_ipw)
}

## Test different cond_threshold values
thresholds <- c(1e8, 1e6, 1e4, 1e3, 1e2, 50, 20)
r_vec <- dat$r

cat("=== Testing cond_threshold (inline GMM) ===\n\n")

for (thresh in thresholds) {
  cat(sprintf("--- cond_threshold = %.0e ---\n", thresh))
  set.seed(42)
  for (i in 1:3) {
    di <- model_info[[i]]
    res <- tryCatch(
      run_gmm_cond(di$design_mat, di$h_x, r_vec, cond_threshold = thresh),
      error = function(e) list(alpha = rep(NA, ncol(di$design_mat)),
                               pi = rep(NA, n), cond = NA, mu = NA)
    )
    ps <- res$pi
    cat(sprintf("  M%d: alpha=(%s)\n", i, paste(round(res$alpha, 3), collapse=", ")))
    cat(sprintf("      PS=[%.4f,%.4f], >0.99:%d, <0.01:%d, cond=%.1e, mu=%.4f\n",
        min(ps, na.rm=TRUE), max(ps, na.rm=TRUE),
        sum(ps > 0.99, na.rm=TRUE), sum(ps < 0.01, na.rm=TRUE),
        res$cond, res$mu))
  }
  cat("\n")
}
