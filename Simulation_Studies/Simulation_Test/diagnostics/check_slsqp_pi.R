## Check if SLSQP constrained GMM addresses pi->1 for y=1 obs in S4 M3
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EPS", quiet = TRUE)
library(nloptr)

n_val <- 2000
ps_spec <- get_ps_spec("9-alt1")
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
compute_pi_fn <- function(eta) plogis(-eta)
tol <- 1e-6; max_outer <- 100L; cond_limit <- 1e8

run_constrained_gmm <- function(dat, design_mat, h_x, link_type) {
  r_vec <- dat[["r"]]
  n <- nrow(dat); k <- ncol(design_mat); h_dim <- ncol(h_x)
  init <- rep(0, k)
  g_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    (r_vec / pi_v - 1) * h_x
  }
  G_bar_fn <- function(alpha) {
    g_mat <- g_fn(alpha)
    matrix(.colMeans(g_mat, n, h_dim), h_dim, 1)
  }
  Gamma_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    crossprod(h_x * (r_vec * (1 - pi_v) / pi_v), design_mat) / n
  }
  constraint_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    mean(r_vec / pi_v - 1)
  }
  constraint_jac <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    as.vector(colMeans(r_vec * (1 - pi_v) / pi_v * design_mat))
  }
  solve_constrained <- function(start, W_hat) {
    obj_fn <- function(alpha) {
      G_v <- G_bar_fn(alpha)
      as.numeric(crossprod(G_v, W_hat %*% G_v))
    }
    grad_fn <- function(alpha) {
      G_v <- G_bar_fn(alpha)
      Gamma <- Gamma_fn(alpha)
      2 * as.vector(crossprod(Gamma, W_hat %*% G_v))
    }
    res <- nloptr::slsqp(x0 = start, fn = obj_fn, gr = grad_fn,
      heq = constraint_fn, heqjac = constraint_jac,
      control = list(maxeval = 2000, xtol_rel = 1e-10, ftol_rel = 1e-12))
    res$par
  }
  estimates <- solve_constrained(init, diag(h_dim))
  for (t in 1:max_outer) {
    g_mat <- g_fn(estimates)
    G_vec <- G_bar_fn(estimates)
    W_hat <- tryCatch(W_fn(g_mat), error = function(e) diag(h_dim))
    Gamma <- Gamma_fn(estimates)
    cg <- 2 * as.vector(crossprod(Gamma, W_hat %*% G_vec))
    if (max(abs(cg)) < tol) break
    estimates <- solve_constrained(estimates, W_hat)
  }
  pi_final <- compute_pi_fn(as.vector(design_mat %*% estimates))
  g_mat_f <- g_fn(estimates)
  W_f <- tryCatch(solve(crossprod(g_mat_f)/n), error = function(e) diag(h_dim))
  fc <- tryCatch({
    Gamma <- Gamma_fn(estimates)
    H <- crossprod(Gamma, W_f %*% Gamma)
    eig <- eigen(H, symmetric=TRUE, only.values=TRUE)$values
    max(eig)/max(min(eig),1e-15)
  }, error = function(e) Inf)
  list(estimates = estimates, pi_hat = pi_final, cond = fc)
}

data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)

cat("Checking pi distribution for y=1 obs in S4 M3 (SLSQP constrained)\n\n")

n_check <- 100
n_degen <- 0; n_nondegen <- 0
degen_pi_y1 <- c(); nondegen_pi_y1 <- c()
degen_pi_y0 <- c(); nondegen_pi_y0 <- c()

for (rep_i in 1:n_check) {
  if (rep_i %% 20 == 0) cat(sprintf("  rep %d...\n", rep_i))
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[3]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[3]]),
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL), optimizer = "L-BFGS-B"
  )
  tryCatch({
    ebmr <- EPS$new("y", single_ps, dat, W_fn)
    di <- list(design_mat = ebmr$ps_fit.list[[1]]$design_matrix,
               h_x = ebmr$ps_fit.list[[1]]$h_x,
               link_type = ebmr$ps_fit.list[[1]]$link_type)
    res <- run_constrained_gmm(dat, di$design_mat, di$h_x, di$link_type)
    y <- dat$y; r <- dat$r
    pi_hat <- res$pi_hat
    is_degen <- res$cond >= cond_limit
    if (is_degen) {
      n_degen <- n_degen + 1
      degen_pi_y1 <- c(degen_pi_y1, pi_hat[y == 1])
      degen_pi_y0 <- c(degen_pi_y0, pi_hat[y == 0])
    } else {
      n_nondegen <- n_nondegen + 1
      nondegen_pi_y1 <- c(nondegen_pi_y1, pi_hat[y == 1])
      nondegen_pi_y0 <- c(nondegen_pi_y0, pi_hat[y == 0])
    }
  }, error = function(e) cat(sprintf("  Rep %d ERROR: %s\n", rep_i, e$message)))
}

cat(sprintf("\nOut of %d reps: %d degenerate, %d non-degenerate\n\n", n_check, n_degen, n_nondegen))

if (n_degen > 0) {
  cat("=== DEGENERATE reps: pi for y=1 obs ===\n")
  cat(sprintf("  n obs: %d\n", length(degen_pi_y1)))
  cat(sprintf("  min=%.6f, Q1=%.6f, median=%.6f, Q3=%.6f, max=%.6f\n",
      min(degen_pi_y1), quantile(degen_pi_y1, 0.25), median(degen_pi_y1),
      quantile(degen_pi_y1, 0.75), max(degen_pi_y1)))
  cat(sprintf("  frac pi > 0.99: %.4f\n", mean(degen_pi_y1 > 0.99)))
  cat(sprintf("  frac pi > 0.999: %.4f\n", mean(degen_pi_y1 > 0.999)))
  cat(sprintf("  frac pi > 0.9999: %.4f\n", mean(degen_pi_y1 > 0.9999)))
  cat("\n=== DEGENERATE reps: pi for y=0 obs ===\n")
  cat(sprintf("  n obs: %d\n", length(degen_pi_y0)))
  cat(sprintf("  min=%.6f, Q1=%.6f, median=%.6f, Q3=%.6f, max=%.6f\n",
      min(degen_pi_y0), quantile(degen_pi_y0, 0.25), median(degen_pi_y0),
      quantile(degen_pi_y0, 0.75), max(degen_pi_y0)))
  cat(sprintf("  E[r/pi] for y=1 in degen: %.4f\n",
      mean(degen_pi_y1)))  # approximate
}

if (n_nondegen > 0) {
  cat("\n=== NON-DEGENERATE reps: pi for y=1 obs ===\n")
  cat(sprintf("  n obs: %d\n", length(nondegen_pi_y1)))
  cat(sprintf("  min=%.6f, Q1=%.6f, median=%.6f, Q3=%.6f, max=%.6f\n",
      min(nondegen_pi_y1), quantile(nondegen_pi_y1, 0.25), median(nondegen_pi_y1),
      quantile(nondegen_pi_y1, 0.75), max(nondegen_pi_y1)))
  cat("\n=== NON-DEGENERATE reps: pi for y=0 obs ===\n")
  cat(sprintf("  n obs: %d\n", length(nondegen_pi_y0)))
  cat(sprintf("  min=%.6f, Q1=%.6f, median=%.6f, Q3=%.6f, max=%.6f\n",
      min(nondegen_pi_y0), quantile(nondegen_pi_y0, 0.25), median(nondegen_pi_y0),
      quantile(nondegen_pi_y0, 0.75), max(nondegen_pi_y0)))
}

# Also check: among degenerate reps, what does r/pi look like for y=1?
if (n_degen > 0) {
  cat("\n=== DEGENERATE reps: constraint check ===\n")
  # Re-run one degenerate rep to get full details
  for (rep_i in 1:n_check) {
    dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
    single_ps <- list(
      formula.list = list(ps_spec[["formula.list"]][[3]]),
      h_alpha.list = list(ps_spec[["h_alpha.list"]][[3]]),
      inv_link = ps_spec[["inv_link"]],
      outcome = ps_spec[["outcome"]],
      alpha_init.list = list(NULL), optimizer = "L-BFGS-B"
    )
    tryCatch({
      ebmr <- EPS$new("y", single_ps, dat, W_fn)
      di <- list(design_mat = ebmr$ps_fit.list[[1]]$design_matrix,
                 h_x = ebmr$ps_fit.list[[1]]$h_x,
                 link_type = ebmr$ps_fit.list[[1]]$link_type)
      res <- run_constrained_gmm(dat, di$design_mat, di$h_x, di$link_type)
      if (res$cond >= cond_limit) {
        y <- dat$y; r <- dat$r; pi_hat <- res$pi_hat
        cat(sprintf("\n  Rep %d (degenerate, cond=%.2e):\n", rep_i, res$cond))
        cat(sprintf("    alpha = %s\n", paste(round(res$estimates, 3), collapse=", ")))
        cat(sprintf("    E[r/pi-1] = %.2e\n", mean(r/pi_hat - 1)))
        cat(sprintf("    E[r/pi-1 | y=1] = %.4f  (n_y1=%d, n_resp_y1=%d)\n",
            mean(r[y==1]/pi_hat[y==1] - 1), sum(y==1), sum(y==1 & r==1)))
        cat(sprintf("    E[r/pi-1 | y=0] = %.4f  (n_y0=%d, n_resp_y0=%d)\n",
            mean(r[y==0]/pi_hat[y==0] - 1), sum(y==0), sum(y==0 & r==0)))
        cat(sprintf("    pi(y=1): min=%.6f, mean=%.6f, max=%.6f\n",
            min(pi_hat[y==1]), mean(pi_hat[y==1]), max(pi_hat[y==1])))
        cat(sprintf("    pi(y=0): min=%.6f, mean=%.6f, max=%.6f\n",
            min(pi_hat[y==0]), mean(pi_hat[y==0]), max(pi_hat[y==0])))
        cat(sprintf("    r/pi(y=1,r=1): min=%.4f, mean=%.4f, max=%.4f\n",
            min(r[y==1&r==1]/pi_hat[y==1&r==1]), mean(r[y==1&r==1]/pi_hat[y==1&r==1]),
            max(r[y==1&r==1]/pi_hat[y==1&r==1])))
        cat(sprintf("    r/pi(y=0,r=1): min=%.4f, mean=%.4f, max=%.4f\n",
            min(r[y==0&r==1]/pi_hat[y==0&r==1]), mean(r[y==0&r==1]/pi_hat[y==0&r==1]),
            max(r[y==0&r==1]/pi_hat[y==0&r==1])))
        break  # just show first degenerate rep
      }
    }, error = function(e) NULL)
  }
}
