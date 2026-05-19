## Test: nloptr L-BFGS for all models + ensemble, compute SE2 and bootstrap SE
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Mental_Health_Data_Application")
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)
library(nloptr); library(Matrix); library(numDeriv)
library(parallel); library(foreach); library(doSNOW)
source("MHD_functions.R")

original_data <- read.csv("data_application.csv")
percent <- original_data$Percentage
n <- 2486
class_n <- round(n * percent / 100)
dat <- gen_data(original_data, class_n, n)
dat$y <- dat$teacher_report

W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
compute_pi <- function(eta) plogis(-pmin(pmax(eta, -20), 20))

h_alpha <- c("health", "father", "parent_report", "fp", "fh", "hp")
h_nu <- function(data) cbind(health=data$health, father=data$father, parent_report=data$parent_report,
                              fp=data$fp, fh=data$fh, hp=data$hp, fhp=data$fhp)

ps_specifications <- list(
  formula.list = list(
    r ~ teacher_report + health + father + teacher_report:health + teacher_report:father,
    r ~ teacher_report + health + parent_report + teacher_report:health + teacher_report:parent_report,
    r ~ teacher_report + parent_report + father + teacher_report:parent_report + teacher_report:father
  ),
  h_alpha.list = list(h_alpha, h_alpha, h_alpha),
  outcome = "teacher_report",
  inv_link = function(eta) 1 / (1 + exp(eta))
)

## Custom GMM inner loop using nloptr L-BFGS
run_gmm_nloptr <- function(dm, hx, r_vec, nn) {
  k <- ncol(dm); h_dim <- ncol(hx)
  g_fn <- function(a) (r_vec / compute_pi(as.vector(dm %*% a)) - 1) * hx
  G_fn <- function(a) { g <- g_fn(a); matrix(.colMeans(g, nn, h_dim), h_dim, 1) }
  Gamma_fn <- function(a) {
    pi_v <- compute_pi(as.vector(dm %*% a))
    crossprod(hx * (r_vec * (1-pi_v)/pi_v), dm) / nn
  }
  obj <- function(a, W) { G <- G_fn(a); as.numeric(crossprod(G, W %*% G)) }
  grad <- function(a, W) { 2 * as.vector(crossprod(Gamma_fn(a), W %*% G_fn(a))) }

  W_hat <- diag(h_dim)
  res <- nloptr(x0 = rep(0, k),
                eval_f = function(x) obj(x, W_hat),
                eval_grad_f = function(x) grad(x, W_hat),
                opts = list(algorithm = "NLOPT_LD_LBFGS", maxeval = 5000, xtol_rel = 1e-12))
  est <- res$solution

  for (t in 1:200) {
    g_mat <- g_fn(est)
    W_hat <- tryCatch(W_fn(g_mat), error = function(e) diag(h_dim))
    gn <- max(abs(grad(est, W_hat)))
    if (gn < 1e-6) break
    res <- nloptr(x0 = est,
                  eval_f = function(x) obj(x, W_hat),
                  eval_grad_f = function(x) grad(x, W_hat),
                  opts = list(algorithm = "NLOPT_LD_LBFGS", maxeval = 5000, xtol_rel = 1e-12))
    est <- res$solution
  }
  list(alpha = est, pi = compute_pi(as.vector(dm %*% est)))
}

## Build design matrix manually
build_dm <- function(d, m) {
  if (m == 1) cbind(1, d$y, d$health, d$father, d$y*d$health, d$y*d$father)
  else if (m == 2) cbind(1, d$y, d$health, d$parent_report, d$y*d$health, d$y*d$parent_report)
  else cbind(1, d$y, d$parent_report, d$father, d$y*d$parent_report, d$y*d$father)
}
build_hx <- function(d) {
  cbind(1, d$health, d$father, d$parent_report, d$fp, d$fh, d$hp)
}

## SE2 computation (CUE)
compute_se2 <- function(dm, hx, r_vec, y_vec, alpha, nn) {
  pi_hat <- compute_pi(as.vector(dm %*% alpha))
  mu_hat <- mean(r_vec * y_vec / pi_hat)
  g_mat <- (r_vec / pi_hat - 1) * hx
  W_hat <- tryCatch(solve(crossprod(g_mat)/nn), error = function(e) NULL)
  if (is.null(W_hat)) return(list(mu = mu_hat, se = NA))
  cf <- r_vec * (1 - pi_hat) / pi_hat
  Gamma <- crossprod(hx * cf, dm) / nn
  GtWG <- crossprod(Gamma, W_hat %*% Gamma)
  k <- length(alpha); h <- ncol(hx)
  GtW <- crossprod(Gamma, W_hat)
  G_vec <- colMeans(g_mat)
  R_mat <- matrix(0, h*k, k)
  for (j in 1:k) { block <- crossprod(cf*dm, dm[,j]*hx)/nn; for (l in 1:h) R_mat[(l-1)*k+j,] <- block[,l] }
  eta_s <- matrix(G_vec, h, 1)
  WG <- as.vector(W_hat %*% eta_s); GW_row <- as.vector(t(eta_s) %*% W_hat)
  S_mat <- matrix(0, h^2, k)
  for (j in 1:k) { dg_j <- cf*dm[,j]*hx; cross <- (crossprod(dg_j,g_mat)+crossprod(g_mat,dg_j))/nn; S_mat[,j] <- as.vector(cross) }
  H_mat <- GtWG + kronecker(matrix(GW_row,1,h), diag(k)) %*% R_mat - kronecker(matrix(GW_row,1,h), GtW) %*% S_mat
  H_cond <- tryCatch({ ev <- eigen(H_mat, symmetric=FALSE, only.values=TRUE)$values; max(Mod(ev))/max(min(Mod(ev)),1e-15) }, error=function(e) Inf)
  if (!is.finite(H_cond) || H_cond > 1e15) return(list(mu = mu_hat, se = NA))
  H_inv <- tryCatch(solve(H_mat), error = function(e) NULL)
  if (is.null(H_inv)) return(list(mu = mu_hat, se = NA))
  GtW_g <- GtW %*% t(g_mat); g_WG <- as.vector(g_mat %*% WG)
  GammaT_WG <- t(dm * (cf * as.vector(hx %*% WG)))
  Q_mat <- GtW_g + GammaT_WG - sweep(GtW_g, 2, g_WG, `*`)
  psi <- -(H_inv %*% Q_mat)
  dot_pi <- -dm * (pi_hat * (1 - pi_hat))
  ry_ps_inv2 <- as.vector(r_vec * y_vec * (pi_hat^(-2)))
  H_alpha <- colMeans(dot_pi * ry_ps_inv2)
  mu_iid <- as.vector(t(r_vec/pi_hat*y_vec) - t(H_alpha) %*% psi)
  se <- sqrt(var(mu_iid)/nn)
  list(mu = mu_hat, se = se, psi = psi)
}

cat("=== nloptr L-BFGS: Point estimates with SE2 ===\n\n")

# Point estimates for each model
point_results <- list()
for (m in 1:3) {
  dm <- build_dm(dat, m); hx <- build_hx(dat)
  res <- run_gmm_nloptr(dm, hx, dat$r, n)
  se_res <- compute_se2(dm, hx, dat$r, dat$teacher_report, res$alpha, n)
  point_results[[m]] <- list(alpha = res$alpha, pi = res$pi, mu = se_res$mu, se = se_res$se)
  cat(sprintf("M%d: mu=%.4f, se=%.4f, PS=[%.4f,%.4f], >0.99:%d\n",
      m, se_res$mu, se_res$se, min(res$pi), max(res$pi), sum(res$pi > 0.99)))
}

# Ensemble using package (with nloptr-estimated alphas injected)
# For now, just use the package directly with its own optimizer for comparison
cat("\n--- Package results (Fast4 L-BFGS-B) for comparison ---\n\n")
ebmr <- EBMRAlgorithmFast4$new("teacher_report", ps_specifications, dat, W_fn)
all_sets <- list(c(1),c(2),c(3),c(1,2),c(1,3),c(2,3),c(1,2,3))
labels <- c("100","010","001","110","101","011","111")
cat(sprintf("%-8s %10s %10s\n", "Label", "Estimate", "SE2"))
cat(paste(rep("-", 30), collapse=""), "\n")
for (j in 1:7) {
  res <- ebmr$EBMR_IPW(h_nu = h_nu, model_indices = all_sets[[j]], true_ps = NULL)
  cat(sprintf("%-8s %10.4f %10.4f\n", labels[j], res$mu_ipw, res$se_ipw))
}

# Bootstrap with nloptr for single models M1, M2, M3
cat("\n=== Bootstrap (B=500): nloptr L-BFGS vs Fast4 ===\n\n")
B <- 500

for (m in 1:3) {
  mus_nl <- mus_f4 <- rep(NA, B)
  for (b in 1:B) {
    set.seed(12345 + b)
    idx <- sample(1:n, n, replace = TRUE)
    dat_b <- dat[idx, ]
    dm_b <- build_dm(dat_b, m); hx_b <- build_hx(dat_b)

    # nloptr
    tryCatch({
      res_nl <- run_gmm_nloptr(dm_b, hx_b, dat_b$r, nrow(dat_b))
      mus_nl[b] <- mean(dat_b$r * dat_b$teacher_report / res_nl$pi)
    }, error = function(e) {})

    # Fast4
    tryCatch({
      ps_spec_b <- list(formula.list = list(ps_specifications$formula.list[[m]]),
                        h_alpha.list = list(h_alpha), outcome = "teacher_report",
                        inv_link = function(eta) 1/(1+exp(eta)))
      ebmr_b <- EBMRAlgorithmFast4$new("teacher_report", ps_spec_b, dat_b, W_fn)
      mus_f4[b] <- mean(dat_b$r * dat_b$teacher_report / ebmr_b$ps_fit.list[[1]]$fitted.values)
    }, error = function(e) {})
  }

  se_analytical <- point_results[[m]]$se
  cat(sprintf("M%d: SE2=%.4f | nloptr boot_sd=%.4f (ratio=%.2f) | Fast4 boot_sd=%.4f (ratio=%.2f)\n",
      m, se_analytical,
      sd(mus_nl, na.rm=TRUE), se_analytical/sd(mus_nl, na.rm=TRUE),
      sd(mus_f4, na.rm=TRUE), se_analytical/sd(mus_f4, na.rm=TRUE)))
}
