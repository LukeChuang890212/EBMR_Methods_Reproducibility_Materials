## Test: Use R's gmm package for alpha estimation
## Compare stability with Fast4's optimizer on 111 ensemble
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Mental_Health_Data_Application")
library(gmm)
library(Matrix)
source("MHD_functions.R")

original_data <- read.csv("data_application.csv")
percent <- original_data$Percentage
n <- 2486
class_n <- round(n * percent / 100)
dat <- gen_data(original_data, class_n, n)
dat$y <- dat$teacher_report

compute_pi <- function(eta) 1 / (1 + exp(eta))  # logistic complement

## Define moment functions for gmm package
## M1: r ~ teacher_report + health + father + teacher_report:health + teacher_report:father
## design_mat = (1, y, health, father, y:health, y:father)
## h_x = (1, health, father, parent_report, fp, fh, hp)
## Moment: g_i(alpha) = (r_i/pi_i - 1) * h_x_i, where pi_i = 1/(1+exp(design_mat_i %*% alpha))

build_design_mat <- function(d) {
  cbind(1, d$y, d$health, d$father, d$y * d$health, d$y * d$father)
}

build_h_x <- function(d) {
  cbind(1, d$health, d$father, d$parent_report,
        d$father * d$parent_report, d$father * d$health,
        d$health * d$parent_report)
}

## Moment function for gmm: g(theta, x)
## x is the data matrix, theta is the parameter vector
## Returns n x h_dim matrix of moment conditions
make_moment_fn <- function(design_mat, h_x, r_vec) {
  function(theta, x) {
    pi_v <- compute_pi(as.vector(design_mat %*% theta))
    as.matrix((r_vec / pi_v - 1) * h_x)
  }
}

## Gradient of moment function for gmm
make_grad_fn <- function(design_mat, h_x, r_vec) {
  n <- length(r_vec)
  k <- ncol(design_mat)
  h_dim <- ncol(h_x)
  function(theta, x) {
    pi_v <- compute_pi(as.vector(design_mat %*% theta))
    cf <- r_vec * (1 - pi_v) / pi_v
    # Return Gamma = E[dg/dtheta], h_dim x k matrix
    crossprod(h_x * cf, design_mat) / n
  }
}

cat("=== Testing gmm package on MHD data ===\n\n")

# Build matrices for M1
dm1 <- build_design_mat(dat)
hx <- build_h_x(dat)
r_vec <- dat$r

cat(sprintf("Design matrix: %d x %d\n", nrow(dm1), ncol(dm1)))
cat(sprintf("h_x: %d x %d\n", nrow(hx), ncol(hx)))
cat(sprintf("overid: %d\n\n", ncol(hx) - ncol(dm1)))

g_fn <- make_moment_fn(dm1, hx, r_vec)
dg_fn <- make_grad_fn(dm1, hx, r_vec)

# Test different gmm optimizers
init <- rep(0, ncol(dm1))
optfcts <- c("optim", "nlminb")
gmm_types <- c("iterative", "cue")

cat("--- Point estimates ---\n\n")

for (gtype in gmm_types) {
  for (optf in optfcts) {
    cat(sprintf("gmm type=%s, optfct=%s: ", gtype, optf))
    tryCatch({
      fit <- gmm(g_fn, x = dat, t0 = init, type = gtype, optfct = optf,
                 gradv = dg_fn,
                 control = list(maxit = 5000, reltol = 1e-10))
      alpha <- coef(fit)
      pi_hat <- compute_pi(as.vector(dm1 %*% alpha))
      mu <- mean(r_vec * dat$teacher_report / pi_hat)
      cat(sprintf("mu=%.4f, alpha=(%s)\n",
          mu, paste(round(alpha, 4), collapse=", ")))
      cat(sprintf("  PS=[%.4f, %.4f], conv=%d, obj=%.2e\n",
          min(pi_hat), max(pi_hat), fit$algoInfo$convergence, fit$objective))
    }, error = function(e) {
      cat(sprintf("ERROR: %s\n", e$message))
    })
  }
}

# Bootstrap comparison: gmm CUE with nlminb vs Fast4
cat("\n--- Bootstrap stability (B=200): gmm CUE nlminb vs optim ---\n\n")

B <- 200

for (optf in c("nlminb", "optim")) {
  cat(sprintf("gmm CUE %s:\n", optf))
  mus <- alphas <- rep(NA, B)
  alpha_mat <- matrix(NA, B, ncol(dm1))

  for (b in 1:B) {
    set.seed(12345 + b)
    idx <- sample(1:n, n, replace = TRUE)
    dat_b <- dat[idx, ]

    dm_b <- build_design_mat(dat_b)
    hx_b <- build_h_x(dat_b)
    r_b <- dat_b$r

    g_b <- make_moment_fn(dm_b, hx_b, r_b)
    dg_b <- make_grad_fn(dm_b, hx_b, r_b)

    tryCatch({
      fit_b <- gmm(g_b, x = dat_b, t0 = init, type = "cue", optfct = optf,
                    gradv = dg_b,
                    control = list(maxit = 5000, reltol = 1e-10))
      alpha_b <- coef(fit_b)
      pi_b <- compute_pi(as.vector(dm_b %*% alpha_b))
      mus[b] <- mean(r_b * dat_b$teacher_report / pi_b)
      alpha_mat[b, ] <- alpha_b
    }, error = function(e) {})
  }

  valid <- !is.na(mus)
  cat(sprintf("  Valid: %d/%d\n", sum(valid), B))
  cat(sprintf("  mu: mean=%.4f, sd(bootstrap SE)=%.4f\n",
      mean(mus, na.rm=TRUE), sd(mus, na.rm=TRUE)))
  cat(sprintf("  Alpha sd: %s\n",
      paste(round(apply(alpha_mat, 2, sd, na.rm=TRUE), 2), collapse=", ")))
  cat("\n")
}

# Also run Fast4 for direct comparison on same bootstrap samples
cat("--- Fast4 L-BFGS-B (same seeds, M1 only) ---\n")
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

W <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
ps_spec_m1 <- list(
  formula.list = list(r ~ teacher_report + health + father + teacher_report:health + teacher_report:father),
  h_alpha.list = list(c("health", "father", "parent_report", "fp", "fh", "hp")),
  outcome = "teacher_report",
  inv_link = function(eta) 1 / (1 + exp(eta))
)
h_nu_fn <- function(data) cbind(health=data$health, father=data$father, parent_report=data$parent_report,
                                 fp=data$fp, fh=data$fh, hp=data$hp, fhp=data$fhp)

mus_f4 <- rep(NA, B)
alpha_mat_f4 <- matrix(NA, B, 6)

for (b in 1:B) {
  set.seed(12345 + b)
  idx <- sample(1:n, n, replace = TRUE)
  dat_b <- dat[idx, ]
  tryCatch({
    ebmr_b <- EBMRAlgorithmFast4$new("teacher_report", ps_spec_m1, dat_b, W)
    res_b <- ebmr_b$EBMR_IPW(h_nu = h_nu_fn, true_ps = NULL, se.fit = FALSE)
    mus_f4[b] <- res_b$mu_ipw
    alpha_mat_f4[b, ] <- ebmr_b$ps_fit.list[[1]]$coefficients
  }, error = function(e) {})
}

valid_f4 <- !is.na(mus_f4)
cat(sprintf("  Valid: %d/%d\n", sum(valid_f4), B))
cat(sprintf("  mu: mean=%.4f, sd(bootstrap SE)=%.4f\n",
    mean(mus_f4, na.rm=TRUE), sd(mus_f4, na.rm=TRUE)))
cat(sprintf("  Alpha sd: %s\n",
    paste(round(apply(alpha_mat_f4, 2, sd, na.rm=TRUE), 2), collapse=", ")))
