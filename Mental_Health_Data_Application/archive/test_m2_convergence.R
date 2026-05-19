## Investigate M2 convergence: track alpha, objective, PS at each outer iteration
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Mental_Health_Data_Application")
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

obj_fixed <- function(a, W) { G <- G_fn(a); as.numeric(crossprod(G, W %*% G)) }
grad_fixed <- function(a, W) { 2 * as.vector(crossprod(Gamma_fn(a), W %*% G_fn(a))) }
cue_obj <- function(a) {
  g_mat <- g_fn(a)
  G <- matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
  W <- tryCatch(solve(crossprod(g_mat)/n), error = function(e) diag(h_dim))
  as.numeric(crossprod(G, W %*% G))
}

## Inner NR solver (same as Fast4 constrained_nr)
nr_inner <- function(start, W_hat, max_iter = 500, trust = 2.0) {
  alpha <- start
  for (iter in 1:max_iter) {
    g_mat <- g_fn(alpha)
    G_vec <- matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
    Gamma <- Gamma_fn(alpha)
    gv <- 2 * as.vector(crossprod(Gamma, W_hat %*% G_vec))
    if (max(abs(gv)) < 1e-8) break
    H_gn <- crossprod(Gamma, W_hat %*% Gamma)
    direction <- tryCatch(solve(H_gn, -gv), error = function(e) -gv)
    sn <- sqrt(sum(direction^2))
    if (sn > trust) direction <- direction * (trust / sn)
    obj_cur <- as.numeric(crossprod(G_vec, W_hat %*% G_vec))
    step <- 1.0; accepted <- FALSE
    for (ls in 1:20) {
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

cat("=== M2 convergence trace ===\n\n")

# Step 1: W = I
alpha <- nr_inner(rep(0, k), diag(h_dim))
pi_v <- compute_pi(as.vector(dm %*% alpha))
mu <- mean(r_vec * dat$teacher_report / pi_v)

cat(sprintf("%-5s %12s %12s %12s %8s %8s %40s\n",
    "Iter", "CUE_obj", "Fixed_obj", "Grad_norm", ">0.99", "<0.01", "alpha"))
cat(paste(rep("-", 130), collapse=""), "\n")

cat(sprintf("W=I  %12.4e %12.4e %12.4e %8d %8d %40s  mu=%.4f\n",
    cue_obj(alpha), obj_fixed(alpha, diag(h_dim)), max(abs(grad_fixed(alpha, diag(h_dim)))),
    sum(pi_v > 0.99), sum(pi_v < 0.01),
    paste(round(alpha, 3), collapse=", "), mu))

# Step 2: Iterate
for (t in 1:50) {
  g_mat <- g_fn(alpha)
  W_hat <- tryCatch(W_fn(g_mat), error = function(e) diag(h_dim))

  cg <- grad_fixed(alpha, W_hat)
  grad_norm <- max(abs(cg))

  alpha_new <- nr_inner(alpha, W_hat)
  pi_v <- compute_pi(as.vector(dm %*% alpha_new))
  mu <- mean(r_vec * dat$teacher_report / pi_v)

  cue_val <- cue_obj(alpha_new)
  fixed_val <- obj_fixed(alpha_new, W_hat)

  cat(sprintf("%4d %12.4e %12.4e %12.4e %8d %8d %40s  mu=%.4f\n",
      t, cue_val, fixed_val, grad_norm,
      sum(pi_v > 0.99), sum(pi_v < 0.01),
      paste(round(alpha_new, 3), collapse=", "), mu))

  if (max(abs(alpha_new - alpha)) < 1e-8) {
    cat("Converged (alpha change < 1e-8)\n")
    break
  }
  if (grad_norm < 1e-6) {
    cat("Converged (grad < 1e-6)\n")
    break
  }

  alpha <- alpha_new
}
