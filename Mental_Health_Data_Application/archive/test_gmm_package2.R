## Test: gmm package on M2 and M3 (the degenerate models)
## Compare alpha stability with Fast4
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

compute_pi <- function(eta) 1 / (1 + exp(eta))

build_h_x <- function(d) {
  cbind(1, d$health, d$father, d$parent_report,
        d$father * d$parent_report, d$father * d$health,
        d$health * d$parent_report)
}

# Design matrices for each model
build_dm <- function(d, model) {
  if (model == 1) {
    cbind(1, d$y, d$health, d$father, d$y * d$health, d$y * d$father)
  } else if (model == 2) {
    cbind(1, d$y, d$health, d$parent_report, d$y * d$health, d$y * d$parent_report)
  } else {
    cbind(1, d$y, d$parent_report, d$father, d$y * d$parent_report, d$y * d$father)
  }
}

make_moment_fn <- function(design_mat, h_x, r_vec) {
  function(theta, x) {
    pi_v <- compute_pi(as.vector(design_mat %*% theta))
    as.matrix((r_vec / pi_v - 1) * h_x)
  }
}

make_grad_fn <- function(design_mat, h_x, r_vec) {
  n <- length(r_vec)
  function(theta, x) {
    pi_v <- compute_pi(as.vector(design_mat %*% theta))
    cf <- r_vec * (1 - pi_v) / pi_v
    crossprod(h_x * cf, design_mat) / n
  }
}

cat("=== Testing gmm package on M1, M2, M3 ===\n\n")

# Point estimates for all 3 models
init <- rep(0, 6)

cat("--- Point estimates ---\n\n")
for (m in 1:3) {
  dm <- build_dm(dat, m)
  hx <- build_h_x(dat)
  g_fn <- make_moment_fn(dm, hx, dat$r)
  dg_fn <- make_grad_fn(dm, hx, dat$r)

  for (optf in c("nlminb", "optim")) {
    tryCatch({
      fit <- gmm(g_fn, x = dat, t0 = init, type = "cue", optfct = optf,
                 gradv = dg_fn,
                 control = list(maxit = 5000, reltol = 1e-12))
      alpha <- coef(fit)
      pi_hat <- compute_pi(as.vector(dm %*% alpha))
      mu <- mean(dat$r * dat$teacher_report / pi_hat)
      cat(sprintf("M%d %s: mu=%.4f, alpha=(%s)\n", m, optf, mu,
          paste(round(alpha, 3), collapse=", ")))
      cat(sprintf("  PS=[%.4f,%.4f], >0.99:%d, <0.01:%d, conv=%d, obj=%.2e\n",
          min(pi_hat), max(pi_hat), sum(pi_hat>0.99), sum(pi_hat<0.01),
          fit$algoInfo$convergence, fit$objective))
    }, error = function(e) {
      cat(sprintf("M%d %s: ERROR: %s\n", m, optf, e$message))
    })
  }
  cat("\n")
}

# Bootstrap comparison: gmm CUE nlminb vs Fast4 on all 3 models
cat("--- Bootstrap stability (B=200) ---\n\n")
B <- 200

for (optf in c("nlminb", "optim")) {
  cat(sprintf("=== gmm CUE %s ===\n", optf))

  for (m in 1:3) {
    mus <- rep(NA, B)
    alpha_mat <- matrix(NA, B, 6)

    for (b in 1:B) {
      set.seed(12345 + b)
      idx <- sample(1:n, n, replace = TRUE)
      dat_b <- dat[idx, ]
      dm_b <- build_dm(dat_b, m)
      hx_b <- build_h_x(dat_b)
      g_b <- make_moment_fn(dm_b, hx_b, dat_b$r)
      dg_b <- make_grad_fn(dm_b, hx_b, dat_b$r)

      tryCatch({
        fit_b <- gmm(g_b, x = dat_b, t0 = init, type = "cue", optfct = optf,
                      gradv = dg_b,
                      control = list(maxit = 5000, reltol = 1e-12))
        alpha_b <- coef(fit_b)
        pi_b <- compute_pi(as.vector(dm_b %*% alpha_b))
        mus[b] <- mean(dat_b$r * dat_b$teacher_report / pi_b)
        alpha_mat[b, ] <- alpha_b
      }, error = function(e) {})
    }

    valid <- !is.na(mus)
    cat(sprintf("  M%d: Valid=%d/%d, mu_sd=%.4f, alpha_sd=(%s)\n",
        m, sum(valid), B, sd(mus, na.rm=TRUE),
        paste(round(apply(alpha_mat, 2, sd, na.rm=TRUE), 2), collapse=", ")))
  }
  cat("\n")
}

# Fast4 comparison (same seeds)
cat("=== Fast4 L-BFGS-B ===\n")
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)
W <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

formulas <- list(
  r ~ teacher_report + health + father + teacher_report:health + teacher_report:father,
  r ~ teacher_report + health + parent_report + teacher_report:health + teacher_report:parent_report,
  r ~ teacher_report + parent_report + father + teacher_report:parent_report + teacher_report:father
)
h_alpha <- c("health", "father", "parent_report", "fp", "fh", "hp")

for (m in 1:3) {
  mus <- rep(NA, B)
  alpha_mat <- matrix(NA, B, 6)

  for (b in 1:B) {
    set.seed(12345 + b)
    idx <- sample(1:n, n, replace = TRUE)
    dat_b <- dat[idx, ]
    ps_spec <- list(formula.list = list(formulas[[m]]),
                    h_alpha.list = list(h_alpha),
                    outcome = "teacher_report",
                    inv_link = function(eta) 1 / (1 + exp(eta)))
    tryCatch({
      ebmr_b <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat_b, W)
      pf <- ebmr_b$ps_fit.list[[1]]
      pi_b <- pf$fitted.values
      mus[b] <- mean(dat_b$r * dat_b$teacher_report / pi_b)
      alpha_mat[b, ] <- pf$coefficients
    }, error = function(e) {})
  }

  valid <- !is.na(mus)
  cat(sprintf("  M%d: Valid=%d/%d, mu_sd=%.4f, alpha_sd=(%s)\n",
      m, sum(valid), B, sd(mus, na.rm=TRUE),
      paste(round(apply(alpha_mat, 2, sd, na.rm=TRUE), 2), collapse=", ")))
}
