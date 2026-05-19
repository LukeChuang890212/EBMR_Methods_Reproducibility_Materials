## Test: Constrain nu >= -0.001 (prevent negative nu)
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
h_alpha <- c("health", "father", "parent_report", "fp", "fh", "hp", "fhp")
h_nu <- function(data) cbind(health=data$health, father=data$father, parent_report=data$parent_report,
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

## Bounded nu solver: nu >= -0.001
solve_nu_bounded <- function(ps_mat, h_x, r_vec, nn) {
  J <- ncol(ps_mat); h_dim <- ncol(h_x)

  obj_f <- function(nu) {
    ps_nu <- as.vector(ps_mat %*% nu)
    g_mat <- as.vector(r_vec / ps_nu - 1) * h_x
    G <- colMeans(g_mat)
    W_hat <- tryCatch(solve(crossprod(g_mat) / nn), error = function(e) diag(h_dim))
    as.numeric(t(G) %*% W_hat %*% G)
  }
  grad_f <- function(nu) {
    ps_nu <- as.vector(ps_mat %*% nu)
    g_mat <- as.vector(r_vec / ps_nu - 1) * h_x
    G <- matrix(colMeans(g_mat), h_dim, 1)
    W_hat <- tryCatch(solve(crossprod(g_mat) / nn), error = function(e) diag(h_dim))
    Gamma <- crossprod(h_x * (-r_vec/(ps_nu^2)), ps_mat) / nn
    2 * as.vector(crossprod(Gamma, W_hat %*% G))
  }

  # W=I first
  W_I <- diag(h_dim)
  obj_I <- function(nu) {
    ps_nu <- as.vector(ps_mat %*% nu)
    g_mat <- as.vector(r_vec / ps_nu - 1) * h_x
    G <- colMeans(g_mat)
    sum(G^2)
  }
  grad_I <- function(nu) {
    ps_nu <- as.vector(ps_mat %*% nu)
    g_mat <- as.vector(r_vec / ps_nu - 1) * h_x
    G <- matrix(colMeans(g_mat), h_dim, 1)
    Gamma <- crossprod(h_x * (-r_vec/(ps_nu^2)), ps_mat) / nn
    2 * as.vector(crossprod(Gamma, G))
  }

  lb <- rep(-0.001, J)

  # Step 1: W=I with bounds
  res <- nloptr(x0 = rep(1/J, J), eval_f = obj_I, eval_grad_f = grad_I,
                lb = lb, opts = list(algorithm = "NLOPT_LD_LBFGS", maxeval = 5000, xtol_rel = 1e-12))
  est <- res$solution

  # Step 2: Iterate W with bounds
  for (t in 1:200) {
    ps_nu <- as.vector(ps_mat %*% est)
    g_mat <- as.vector(r_vec / ps_nu - 1) * h_x
    W_hat <- tryCatch(solve(crossprod(g_mat) / nn), error = function(e) diag(h_dim))

    obj_W <- function(nu) {
      ps_nu <- as.vector(ps_mat %*% nu)
      g_mat <- as.vector(r_vec / ps_nu - 1) * h_x
      G <- colMeans(g_mat)
      as.numeric(t(G) %*% W_hat %*% G)
    }
    grad_W <- function(nu) {
      ps_nu <- as.vector(ps_mat %*% nu)
      g_mat <- as.vector(r_vec / ps_nu - 1) * h_x
      G <- matrix(colMeans(g_mat), h_dim, 1)
      Gamma <- crossprod(h_x * (-r_vec/(ps_nu^2)), ps_mat) / nn
      2 * as.vector(crossprod(Gamma, W_hat %*% G))
    }

    gn <- max(abs(grad_W(est)))
    if (gn < 1e-7) break

    res <- nloptr(x0 = est, eval_f = obj_W, eval_grad_f = grad_W,
                  lb = lb, opts = list(algorithm = "NLOPT_LD_LBFGS", maxeval = 5000, xtol_rel = 1e-12))
    est <- res$solution
  }

  w <- est^2 / sum(est^2)
  list(nu = est, w = w)
}

# Analytical SE
ebmr_pt <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat, W_fn)
res_pt <- ebmr_pt$EBMR_IPW(h_nu = h_nu, true_ps = NULL)
anal_se <- res_pt$se_ipw

B <- 300
cat(sprintf("=== Bounded nu >= -0.001 (B=%d) ===\n\n", B))

mus_bounded <- mus_default <- rep(NA, B)

for (b in 1:B) {
  if (b %% 100 == 0) cat(sprintf("  rep %d...\n", b))
  set.seed(12345 + b)
  idx <- sample(1:n, n, replace = TRUE)
  dat_b <- dat[idx, ]

  tryCatch({
    ebmr_b <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat_b, W_fn)

    # Default
    res_d <- ebmr_b$EBMR_IPW(h_nu = h_nu, true_ps = NULL, se.fit = FALSE)
    mus_default[b] <- res_d$mu_ipw

    # Bounded
    ps_mat <- do.call(cbind, lapply(ebmr_b$ps_fit.list, function(pf) pf$fitted.values))
    h_x <- cbind(1, h_nu(dat_b))
    n_b <- nrow(dat_b)
    nu_res <- solve_nu_bounded(ps_mat, h_x, dat_b$r, n_b)
    ps_nu <- as.vector(ps_mat %*% nu_res$nu)
    mus_bounded[b] <- mean(dat_b$r * dat_b$teacher_report / ps_nu)
  }, error = function(e) {})
}

cat(sprintf("\nDefault:  boot_se=%.4f, ratio=%.3f\n",
    sd(mus_default, na.rm=TRUE), anal_se / sd(mus_default, na.rm=TRUE)))
cat(sprintf("Bounded:  boot_se=%.4f, ratio=%.3f\n",
    sd(mus_bounded, na.rm=TRUE), anal_se / sd(mus_bounded, na.rm=TRUE)))
