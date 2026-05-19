## Test: stricter convergence (1e-8) on the bad reps
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Mental_Health_Data_Application")
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)
library(Matrix); library(numDeriv); library(nloptr)
source("MHD_functions.R")

original_data <- read.csv("data_application.csv")
percent <- original_data$Percentage
n <- 2486
class_n <- round(n * percent / 100)
dat <- gen_data(original_data, class_n, n)
dat$y <- dat$teacher_report

W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
compute_pi <- function(eta) plogis(-pmin(pmax(eta, -20), 20))

h_alpha_names <- c("health", "father", "parent_report", "fp", "fh", "hp", "fhp")

# Build design matrix and h_x manually for a given model
build_dm <- function(d, m) {
  if (m == 1) cbind(1, d$y, d$health, d$father, d$y*d$health, d$y*d$father, d$father*d$health)
  else if (m == 2) cbind(1, d$y, d$parent_report, d$health, d$y*d$parent_report, d$y*d$health, d$parent_report*d$health)
  else cbind(1, d$y, d$parent_report, d$father, d$y*d$parent_report, d$y*d$father, d$parent_report*d$father)
}

build_hx <- function(d) {
  h1 <- d[h_alpha_names]
  for (j in 1:ncol(h1)) h1[, j] <- as.factor(h1[, j])
  cbind(model.matrix(lm(rep(1, nrow(d)) ~ ., data = h1)))
}

## Strict GMM: iterate until grad < tol
run_strict_gmm <- function(dm, hx, r_vec, nn, tol = 1e-8, max_outer = 1000) {
  k <- ncol(dm); h_dim <- ncol(hx)
  g_fn <- function(a) (r_vec / compute_pi(as.vector(dm %*% a)) - 1) * hx
  G_fn <- function(a) { g <- g_fn(a); matrix(.colMeans(g, nn, h_dim), h_dim, 1) }
  Gamma_fn <- function(a) {
    pi_v <- compute_pi(as.vector(dm %*% a))
    crossprod(hx * (r_vec * (1-pi_v)/pi_v), dm) / nn
  }
  obj_f <- function(a, W) { G <- G_fn(a); as.numeric(crossprod(G, W %*% G)) }
  grad_f <- function(a, W) { 2 * as.vector(crossprod(Gamma_fn(a), W %*% G_fn(a))) }

  # W=I step
  W_hat <- diag(h_dim)
  res <- nloptr(x0 = rep(0, k), eval_f = function(x) obj_f(x, W_hat),
                eval_grad_f = function(x) grad_f(x, W_hat),
                opts = list(algorithm = "NLOPT_LD_LBFGS", maxeval = 10000, xtol_rel = 1e-14))
  est <- res$solution

  # Iterate with strict tolerance
  for (t in 1:max_outer) {
    g_mat <- g_fn(est)
    W_hat <- tryCatch(W_fn(g_mat), error = function(e) diag(h_dim))
    gn <- max(abs(grad_f(est, W_hat)))
    if (gn < tol) break
    res <- nloptr(x0 = est, eval_f = function(x) obj_f(x, W_hat),
                  eval_grad_f = function(x) grad_f(x, W_hat),
                  opts = list(algorithm = "NLOPT_LD_LBFGS", maxeval = 10000, xtol_rel = 1e-14))
    est <- res$solution
  }
  list(alpha = est, pi = compute_pi(as.vector(dm %*% est)),
       grad = max(abs(grad_f(est, W_hat))), obj = obj_f(est, W_hat), iter = t)
}

## Strict nu GMM
run_strict_nu <- function(ps_mat, hx_nu, r_vec, nn, tol = 1e-8, max_outer = 1000) {
  J <- ncol(ps_mat); h_dim <- ncol(hx_nu)
  g_fn <- function(nu) { ps_nu <- as.vector(ps_mat %*% nu); as.vector(r_vec/ps_nu - 1) * hx_nu }
  G_fn <- function(nu) { g <- g_fn(nu); matrix(.colMeans(g, nn, h_dim), h_dim, 1) }
  Gamma_fn <- function(nu) {
    ps_nu <- as.vector(ps_mat %*% nu)
    crossprod(hx_nu * (-r_vec/(ps_nu^2)), ps_mat) / nn
  }
  obj_f <- function(nu, W) { G <- G_fn(nu); as.numeric(crossprod(G, W %*% G)) }
  grad_f <- function(nu, W) { 2 * as.vector(crossprod(Gamma_fn(nu), W %*% G_fn(nu))) }

  W_hat <- diag(h_dim)
  est <- rep(1/J, J)
  res <- nloptr(x0 = est, eval_f = function(x) obj_f(x, W_hat),
                eval_grad_f = function(x) grad_f(x, W_hat),
                opts = list(algorithm = "NLOPT_LD_LBFGS", maxeval = 10000, xtol_rel = 1e-14))
  est <- res$solution

  for (t in 1:max_outer) {
    g_mat <- g_fn(est)
    W_hat <- tryCatch(W_fn(g_mat), error = function(e) diag(h_dim))
    gn <- max(abs(grad_f(est, W_hat)))
    if (gn < tol) break
    res <- nloptr(x0 = est, eval_f = function(x) obj_f(x, W_hat),
                  eval_grad_f = function(x) grad_f(x, W_hat),
                  opts = list(algorithm = "NLOPT_LD_LBFGS", maxeval = 10000, xtol_rel = 1e-14))
    est <- res$solution
  }
  w <- est^2 / sum(est^2)
  list(nu = est, w = w, grad = max(abs(grad_f(est, W_hat))),
       obj = obj_f(est, W_hat), iter = t)
}

# Bad reps from previous scan
bad_reps <- c(1, 2, 6, 7, 8)
good_reps <- c(3, 4, 5)

cat("=== Strict convergence (tol=1e-8) on bad and good reps ===\n\n")

for (b in c(bad_reps, good_reps)) {
  set.seed(12345 + b)
  idx <- sample(1:n, n, replace = TRUE)
  dat_b <- dat[idx, ]
  n_b <- nrow(dat_b)
  hx_b <- build_hx(dat_b)
  h_nu_mat <- cbind(1, dat_b$health, dat_b$father, dat_b$parent_report,
                     dat_b$fp, dat_b$fh, dat_b$hp, dat_b$fhp)

  cat(sprintf("--- Rep %d ---\n", b))

  # Fit each model with strict convergence
  ps_list <- list()
  for (m in 1:3) {
    dm_b <- build_dm(dat_b, m)
    res_m <- run_strict_gmm(dm_b, hx_b, dat_b$r, n_b, tol = 1e-8)
    ps_list[[m]] <- res_m$pi
    cat(sprintf("  M%d: obj=%.2e, grad=%.2e, iter=%d, PS=[%.4f,%.4f], >0.99:%d\n",
        m, res_m$obj, res_m$grad, res_m$iter, min(res_m$pi), max(res_m$pi), sum(res_m$pi > 0.99)))
  }

  # Strict nu
  ps_mat <- do.call(cbind, ps_list)
  nu_res <- run_strict_nu(ps_mat, h_nu_mat, dat_b$r, n_b, tol = 1e-8)
  ps_nu <- as.vector(ps_mat %*% nu_res$nu)
  mu <- mean(dat_b$r * dat_b$teacher_report / ps_nu)
  cat(sprintf("  Nu: obj=%.2e, grad=%.2e, iter=%d, w=(%s), mu=%.4f\n\n",
      nu_res$obj, nu_res$grad, nu_res$iter,
      paste(round(nu_res$w, 4), collapse=", "), mu))
}
