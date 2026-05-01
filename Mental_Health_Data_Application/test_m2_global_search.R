## Global search for M2 optimal solution
## Try many starting points + multiple optimizers to find the true global minimum
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Mental_Health_Data_Application")
library(nloptr); library(Matrix)
source("MHD_functions.R")

original_data <- read.csv("data_application.csv")
percent <- original_data$Percentage
n <- 2486
class_n <- round(n * percent / 100)
dat <- gen_data(original_data, class_n, n)
dat$y <- dat$teacher_report

compute_pi <- function(eta) plogis(-pmin(pmax(eta, -20), 20))
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

# M2: r ~ teacher_report + health + parent_report + teacher_report:health + teacher_report:parent_report
dm <- cbind(1, dat$y, dat$health, dat$parent_report, dat$y*dat$health, dat$y*dat$parent_report)
hx <- cbind(1, dat$health, dat$father, dat$parent_report, dat$fp, dat$fh, dat$hp)
r_vec <- dat$r
k <- ncol(dm); h_dim <- ncol(hx)

g_fn <- function(a) (r_vec / compute_pi(as.vector(dm %*% a)) - 1) * hx
G_fn <- function(a) { g <- g_fn(a); matrix(.colMeans(g, n, h_dim), h_dim, 1) }
Gamma_fn <- function(a) {
  pi_v <- compute_pi(as.vector(dm %*% a))
  crossprod(hx * (r_vec * (1-pi_v)/pi_v), dm) / n
}

## Compute CUE objective: G(a)' W(a) G(a), W evaluated at a
cue_obj <- function(a) {
  g_mat <- g_fn(a)
  G <- matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
  W <- tryCatch(solve(crossprod(g_mat)/n), error = function(e) diag(h_dim))
  as.numeric(crossprod(G, W %*% G))
}

## Fixed-W objective and gradient
obj_fixed <- function(a, W) { G <- G_fn(a); as.numeric(crossprod(G, W %*% G)) }
grad_fixed <- function(a, W) { 2 * as.vector(crossprod(Gamma_fn(a), W %*% G_fn(a))) }

## Iterative GMM with given optimizer and starting point
run_iterative_gmm <- function(init, algo = "NLOPT_LD_LBFGS", max_outer = 200) {
  W_hat <- diag(h_dim)
  est <- tryCatch({
    res <- nloptr(x0 = init, eval_f = function(x) obj_fixed(x, W_hat),
                  eval_grad_f = function(x) grad_fixed(x, W_hat),
                  opts = list(algorithm = algo, maxeval = 5000, xtol_rel = 1e-12))
    res$solution
  }, error = function(e) init)

  for (t in 1:max_outer) {
    g_mat <- g_fn(est)
    W_hat <- tryCatch(W_fn(g_mat), error = function(e) diag(h_dim))
    gn <- max(abs(grad_fixed(est, W_hat)))
    if (gn < 1e-7) break
    est <- tryCatch({
      res <- nloptr(x0 = est, eval_f = function(x) obj_fixed(x, W_hat),
                    eval_grad_f = function(x) grad_fixed(x, W_hat),
                    opts = list(algorithm = algo, maxeval = 5000, xtol_rel = 1e-12))
      res$solution
    }, error = function(e) est)
  }
  est
}

cat("=== Global search for M2 optimal solution ===\n\n")

# Generate many starting points
set.seed(42)
starts <- list(rep(0, k))  # zero
for (i in 1:50) starts[[i+1]] <- rnorm(k) * 2  # random
for (i in 1:20) starts[[i+51]] <- rnorm(k) * 5  # wider random
for (i in 1:10) starts[[i+71]] <- rnorm(k) * 10  # very wide
# Also try grid on alpha[1] and alpha[2] (intercept and teacher_report coef)
for (a1 in seq(-3, 3, by = 1)) {
  for (a2 in seq(-3, 3, by = 1)) {
    starts[[length(starts)+1]] <- c(a1, a2, 0, 0, 0, 0)
  }
}

cat(sprintf("Total starting points: %d\n\n", length(starts)))

# Run from each start with multiple algorithms
results <- list()
algos <- c("NLOPT_LD_LBFGS", "NLOPT_LD_MMA")

for (algo in algos) {
  cat(sprintf("--- %s ---\n", algo))
  best_obj <- Inf
  best_alpha <- NULL

  for (s in seq_along(starts)) {
    tryCatch({
      alpha <- run_iterative_gmm(starts[[s]], algo = algo)
      obj <- cue_obj(alpha)
      pi_hat <- compute_pi(as.vector(dm %*% alpha))
      mu <- mean(r_vec * dat$teacher_report / pi_hat)

      if (obj < best_obj) {
        best_obj <- obj
        best_alpha <- alpha
        best_mu <- mu
        best_ps <- pi_hat
        cat(sprintf("  Start %3d: NEW BEST obj=%.6e, mu=%.4f, PS=[%.4f,%.4f], >0.99:%d\n",
            s, obj, mu, min(pi_hat), max(pi_hat), sum(pi_hat > 0.99)))
      }
    }, error = function(e) {})
  }

  cat(sprintf("\n  BEST: obj=%.6e, mu=%.4f, alpha=(%s)\n",
      best_obj, best_mu, paste(round(best_alpha, 4), collapse=", ")))
  cat(sprintf("  PS=[%.4f,%.4f], >0.99:%d, <0.01:%d\n\n",
      min(best_ps), max(best_ps), sum(best_ps > 0.99), sum(best_ps < 0.01)))

  results[[algo]] <- list(alpha = best_alpha, obj = best_obj, mu = best_mu, ps = best_ps)
}

# Also try R's optim with multiple methods
cat("--- R optim methods ---\n")
for (method in c("L-BFGS-B", "Nelder-Mead", "BFGS", "CG")) {
  best_obj <- Inf; best_alpha <- NULL
  for (s in seq_along(starts)) {
    tryCatch({
      # W=I first
      W_hat <- diag(h_dim)
      res <- optim(starts[[s]], function(x) obj_fixed(x, W_hat),
                   function(x) grad_fixed(x, W_hat), method = method,
                   control = list(maxit = 5000))
      est <- res$par
      # Iterate
      for (t in 1:100) {
        g_mat <- g_fn(est)
        W_hat <- tryCatch(W_fn(g_mat), error = function(e) diag(h_dim))
        gn <- max(abs(grad_fixed(est, W_hat)))
        if (gn < 1e-7) break
        res <- optim(est, function(x) obj_fixed(x, W_hat),
                     function(x) grad_fixed(x, W_hat), method = method,
                     control = list(maxit = 5000))
        est <- res$par
      }
      obj <- cue_obj(est)
      if (obj < best_obj) {
        best_obj <- obj
        best_alpha <- est
        pi_hat <- compute_pi(as.vector(dm %*% est))
        best_mu <- mean(r_vec * dat$teacher_report / pi_hat)
        best_ps <- pi_hat
      }
    }, error = function(e) {})
  }
  cat(sprintf("  %12s: obj=%.6e, mu=%.4f, PS=[%.4f,%.4f], >0.99:%d\n",
      method, best_obj, best_mu, min(best_ps), max(best_ps), sum(best_ps > 0.99)))
}
