## Test: Replace inner-loop optimizer with nloptr (SLSQP, augmented Lagrangian)
## to see if it produces more stable M2/M3 alpha estimates
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Mental_Health_Data_Application")
library(nloptr)
library(Matrix)
source("MHD_functions.R")

original_data <- read.csv("data_application.csv")
percent <- original_data$Percentage
n <- 2486
class_n <- round(n * percent / 100)
dat <- gen_data(original_data, class_n, n)
dat$y <- dat$teacher_report

compute_pi <- function(eta) plogis(-pmin(pmax(eta, -20), 20))
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

build_h_x <- function(d) {
  cbind(1, d$health, d$father, d$parent_report,
        d$father * d$parent_report, d$father * d$health,
        d$health * d$parent_report)
}

build_dm <- function(d, model) {
  if (model == 1) cbind(1, d$y, d$health, d$father, d$y*d$health, d$y*d$father)
  else if (model == 2) cbind(1, d$y, d$health, d$parent_report, d$y*d$health, d$y*d$parent_report)
  else cbind(1, d$y, d$parent_report, d$father, d$y*d$parent_report, d$y*d$father)
}

## Iterative GMM using nloptr for the inner loop
run_gmm_nloptr <- function(dm, hx, r_vec, n, algorithm = "NLOPT_LD_SLSQP", max_outer = 200) {
  k <- ncol(dm); h_dim <- ncol(hx)

  g_fn <- function(alpha) {
    pi_v <- compute_pi(as.vector(dm %*% alpha))
    (r_vec / pi_v - 1) * hx
  }
  G_fn <- function(alpha) {
    g_mat <- g_fn(alpha)
    matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
  }
  Gamma_fn <- function(alpha) {
    pi_v <- compute_pi(as.vector(dm %*% alpha))
    cf <- r_vec * (1 - pi_v) / pi_v
    crossprod(hx * cf, dm) / n
  }

  # Inner objective: Q(alpha) = G(alpha)' W G(alpha)
  obj_fn <- function(alpha, W_hat) {
    G <- G_fn(alpha)
    as.numeric(crossprod(G, W_hat %*% G))
  }
  grad_fn <- function(alpha, W_hat) {
    G <- G_fn(alpha)
    Gamma <- Gamma_fn(alpha)
    2 * as.vector(crossprod(Gamma, W_hat %*% G))
  }

  # Step 1: W = I
  W_hat <- diag(h_dim)
  res <- nloptr(x0 = rep(0, k),
                eval_f = function(x) obj_fn(x, W_hat),
                eval_grad_f = function(x) grad_fn(x, W_hat),
                opts = list(algorithm = algorithm, maxeval = 5000, xtol_rel = 1e-12))
  estimates <- res$solution

  # Step 2: Iterative W updates
  best_grad <- Inf; best_est <- estimates
  for (t in 1:max_outer) {
    g_mat <- g_fn(estimates)
    W_hat <- tryCatch(W_fn(g_mat), error = function(e) diag(h_dim))
    cg <- grad_fn(estimates, W_hat)
    gn <- max(abs(cg))
    if (gn < 1e-6) break
    if (gn < best_grad) { best_grad <- gn; best_est <- estimates }

    res <- nloptr(x0 = estimates,
                  eval_f = function(x) obj_fn(x, W_hat),
                  eval_grad_f = function(x) grad_fn(x, W_hat),
                  opts = list(algorithm = algorithm, maxeval = 5000, xtol_rel = 1e-12))
    estimates <- res$solution
  }

  pi_final <- compute_pi(as.vector(dm %*% estimates))
  g_f <- g_fn(estimates)
  W_f <- tryCatch(W_fn(g_f), error = function(e) diag(h_dim))
  final_grad <- max(abs(grad_fn(estimates, W_f)))
  obj_final <- obj_fn(estimates, W_f)

  list(alpha = estimates, pi = pi_final, grad = final_grad, obj = obj_final)
}

cat("=== Testing nloptr optimizers on M1, M2, M3 ===\n\n")

# Available gradient-based algorithms in nloptr
algos <- c("NLOPT_LD_SLSQP", "NLOPT_LD_MMA", "NLOPT_LD_LBFGS",
           "NLOPT_LD_TNEWTON_PRECOND_RESTART", "NLOPT_LD_VAR2")

cat("--- Point estimates ---\n\n")
for (m in 1:3) {
  dm <- build_dm(dat, m)
  hx <- build_h_x(dat)
  cat(sprintf("Model %d:\n", m))

  for (algo in algos) {
    tryCatch({
      res <- run_gmm_nloptr(dm, hx, dat$r, n, algorithm = algo)
      cat(sprintf("  %s: mu=%.4f, PS=[%.4f,%.4f], >0.99:%d, <0.01:%d, grad=%.1e, obj=%.1e\n",
          algo, mean(dat$r * dat$teacher_report / res$pi),
          min(res$pi), max(res$pi), sum(res$pi > 0.99), sum(res$pi < 0.01),
          res$grad, res$obj))
    }, error = function(e) {
      cat(sprintf("  %s: ERROR: %s\n", algo, e$message))
    })
  }

  # Also run Fast4 L-BFGS-B for comparison
  tryCatch({
    devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)
    ps_spec <- list(
      formula.list = list(switch(m,
        r ~ teacher_report + health + father + teacher_report:health + teacher_report:father,
        r ~ teacher_report + health + parent_report + teacher_report:health + teacher_report:parent_report,
        r ~ teacher_report + parent_report + father + teacher_report:parent_report + teacher_report:father)),
      h_alpha.list = list(c("health", "father", "parent_report", "fp", "fh", "hp")),
      outcome = "teacher_report",
      inv_link = function(eta) 1 / (1 + exp(eta)))
    W <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
    ebmr <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat, W)
    pf <- ebmr$ps_fit.list[[1]]
    ps <- pf$fitted.values
    cat(sprintf("  Fast4 L-BFGS-B: mu=%.4f, PS=[%.4f,%.4f], >0.99:%d, <0.01:%d, grad=%.1e, obj=%.1e\n",
        mean(dat$r * dat$teacher_report / ps), min(ps), max(ps),
        sum(ps > 0.99), sum(ps < 0.01), pf$gmm_fit$opt$final_grad_norm, pf$gmm_fit$opt$objective))
  }, error = function(e) {
    cat(sprintf("  Fast4: ERROR: %s\n", e$message))
  })
  cat("\n")
}

# Bootstrap stability for the best nloptr algo on M2 and M3
cat("--- Bootstrap stability (B=200) for M2 and M3 ---\n\n")
B <- 200

for (m in 2:3) {
  cat(sprintf("Model %d:\n", m))

  for (algo in c("NLOPT_LD_SLSQP", "NLOPT_LD_LBFGS")) {
    mus <- rep(NA, B)
    alpha_mat <- matrix(NA, B, 6)

    for (b in 1:B) {
      set.seed(12345 + b)
      idx <- sample(1:n, n, replace = TRUE)
      dat_b <- dat[idx, ]
      dm_b <- build_dm(dat_b, m)
      hx_b <- build_h_x(dat_b)

      tryCatch({
        res <- run_gmm_nloptr(dm_b, hx_b, dat_b$r, nrow(dat_b), algorithm = algo)
        mus[b] <- mean(dat_b$r * dat_b$teacher_report / res$pi)
        alpha_mat[b, ] <- res$alpha
      }, error = function(e) {})
    }

    valid <- !is.na(mus)
    cat(sprintf("  %s: Valid=%d/%d, mu_sd=%.4f, alpha_sd=(%s)\n",
        algo, sum(valid), B, sd(mus, na.rm=TRUE),
        paste(round(apply(alpha_mat, 2, sd, na.rm=TRUE), 2), collapse=", ")))
  }

  # Fast4 comparison
  mus_f4 <- rep(NA, B)
  alpha_f4 <- matrix(NA, B, 6)
  for (b in 1:B) {
    set.seed(12345 + b)
    idx <- sample(1:n, n, replace = TRUE)
    dat_b <- dat[idx, ]
    tryCatch({
      ps_spec <- list(
        formula.list = list(switch(m,
          NULL,
          r ~ teacher_report + health + parent_report + teacher_report:health + teacher_report:parent_report,
          r ~ teacher_report + parent_report + father + teacher_report:parent_report + teacher_report:father)),
        h_alpha.list = list(c("health", "father", "parent_report", "fp", "fh", "hp")),
        outcome = "teacher_report",
        inv_link = function(eta) 1 / (1 + exp(eta)))
      ebmr_b <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat_b, W)
      pf <- ebmr_b$ps_fit.list[[1]]
      mus_f4[b] <- mean(dat_b$r * dat_b$teacher_report / pf$fitted.values)
      alpha_f4[b, ] <- pf$coefficients
    }, error = function(e) {})
  }
  valid_f4 <- !is.na(mus_f4)
  cat(sprintf("  Fast4 L-BFGS-B: Valid=%d/%d, mu_sd=%.4f, alpha_sd=(%s)\n",
      sum(valid_f4), B, sd(mus_f4, na.rm=TRUE),
      paste(round(apply(alpha_f4, 2, sd, na.rm=TRUE), 2), collapse=", ")))
  cat("\n")
}
