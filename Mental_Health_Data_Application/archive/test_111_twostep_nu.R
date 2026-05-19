## Test: Two-step nu GMM (fix W after 1 update) vs iterative
## Also test: W=I only (no W update)
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
h_alpha <- c("health", "father", "parent_report", "fp", "fh", "hp", "fhp")
h_nu_fn <- function(data) cbind(health=data$health, father=data$father, parent_report=data$parent_report,
                                  fp=data$fp, fh=data$fh, hp=data$hp, fhp=data$fhp)

ps_spec <- list(
  formula.list = list(
    r ~ teacher_report + health + father + teacher_report:health + teacher_report:father,
    r ~ teacher_report + health + parent_report + teacher_report:health + teacher_report:parent_report,
    r ~ teacher_report + parent_report + father + teacher_report:parent_report + teacher_report:father
  ),
  h_alpha.list = list(h_alpha, h_alpha, h_alpha),
  outcome = "teacher_report", inv_link = function(eta) 1/(1+exp(eta)), optimizer = "L-BFGS-B"
)

## Custom nu GMM with controlled W iterations
run_nu_gmm <- function(ps_mat, h_x, r_vec, nn, max_w_updates = 200) {
  J <- ncol(ps_mat); h_dim <- ncol(h_x)
  g_fn <- function(nu) as.vector(r_vec / as.vector(ps_mat %*% nu) - 1) * h_x
  G_fn <- function(nu) { g <- g_fn(nu); matrix(.colMeans(g, nn, h_dim), h_dim, 1) }
  Gamma_fn <- function(nu) {
    ps_nu <- as.vector(ps_mat %*% nu)
    crossprod(h_x * (-r_vec/(ps_nu^2)), ps_mat) / nn
  }
  obj_f <- function(nu, W) { G <- G_fn(nu); as.numeric(crossprod(G, W %*% G)) }
  grad_f <- function(nu, W) { 2 * as.vector(crossprod(Gamma_fn(nu), W %*% G_fn(nu))) }

  # Step 1: W=I
  W_hat <- diag(h_dim)
  res <- nloptr(x0 = rep(1/J, J), eval_f = function(x) obj_f(x, W_hat),
                eval_grad_f = function(x) grad_f(x, W_hat),
                opts = list(algorithm = "NLOPT_LD_LBFGS", maxeval = 5000, xtol_rel = 1e-12))
  est <- res$solution

  # Step 2+: iterate W
  for (t in 1:max_w_updates) {
    g_mat <- g_fn(est)
    W_hat <- tryCatch(W_fn(g_mat), error = function(e) diag(h_dim))
    gn <- max(abs(grad_f(est, W_hat)))
    if (gn < 1e-7) break
    res <- nloptr(x0 = est, eval_f = function(x) obj_f(x, W_hat),
                  eval_grad_f = function(x) grad_f(x, W_hat),
                  opts = list(algorithm = "NLOPT_LD_LBFGS", maxeval = 5000, xtol_rel = 1e-12))
    est <- res$solution
  }
  w <- est^2 / sum(est^2)
  list(nu = est, w = w)
}

B <- 300
cat("=== Testing nu W-update strategies (B=300) ===\n\n")

strategies <- list(
  list(name = "W=I only (0 updates)", max_w = 0),
  list(name = "Two-step (1 update)", max_w = 1),
  list(name = "Three-step (2 updates)", max_w = 2),
  list(name = "Five-step (4 updates)", max_w = 4),
  list(name = "Full iterative (200)", max_w = 200)
)

# Also get the package default (for analytical SE)
ebmr_pt <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat, W_fn)
res_pt <- ebmr_pt$EBMR_IPW(h_nu = h_nu_fn, true_ps = NULL)
analytical_se <- res_pt$se_ipw

for (strat in strategies) {
  cat(sprintf("--- %s ---\n", strat$name))

  mus <- rep(NA, B)
  w1s <- rep(NA, B)

  for (b in 1:B) {
    set.seed(12345 + b)
    idx <- sample(1:n, n, replace = TRUE)
    dat_b <- dat[idx, ]

    tryCatch({
      ebmr_b <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat_b, W_fn)
      ps_mat <- do.call(cbind, lapply(ebmr_b$ps_fit.list, function(pf) pf$fitted.values))
      h_x <- cbind(1, h_nu_fn(dat_b))
      n_b <- nrow(dat_b)

      nu_res <- run_nu_gmm(ps_mat, h_x, dat_b$r, n_b, max_w_updates = strat$max_w)
      ps_nu <- as.vector(ps_mat %*% nu_res$nu)
      mus[b] <- mean(dat_b$r * dat_b$teacher_report / ps_nu)
      w1s[b] <- nu_res$w[1]
    }, error = function(e) {})
  }

  valid <- !is.na(mus)
  boot_se <- sd(mus, na.rm = TRUE)
  cat(sprintf("  boot_se=%.4f, ratio=%.3f, w1>0.9: %d/%d (%.1f%%)\n\n",
      boot_se, analytical_se / boot_se,
      sum(w1s > 0.9, na.rm = TRUE), sum(valid), 100*mean(w1s > 0.9, na.rm = TRUE)))
}
