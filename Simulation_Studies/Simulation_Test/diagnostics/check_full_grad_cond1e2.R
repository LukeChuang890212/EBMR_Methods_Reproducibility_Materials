## Check full dQ/dα at NR cond=1e2 solution
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EPS", quiet = TRUE)
library(numDeriv)

ps_spec_base <- get_ps_spec("7")
h_alpha_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2)

data_file <- correct_model_all_data_file.list[["setting3"]][["miss50"]]
all_data <- NULL
for (fi in seq_along(data_file)) {
  if (file.exists(data_file[[fi]])) {
    test <- readRDS(data_file[[fi]])
    if (nrow(test) / 1000 == 2000) { all_data <- test; break }
  }
}

nn <- 2000
design_formula <- ps_spec_base$formula.list[[2]]
compute_pi <- function(eta) 1/(1+exp(eta))

g_fn_maker <- function(dat) {
  design_mat <- model.matrix(design_formula, dat)
  r_vec <- as.vector(dat$r)
  h_x <- cbind(1, h_alpha_fn(dat))
  function(alpha) {
    eta <- pmin(pmax(as.vector(design_mat %*% alpha), -20), 20)
    pi_vec <- compute_pi(eta)
    as.vector(r_vec/pi_vec - 1) * h_x
  }
}

# Full CUE objective
Q_cue_maker <- function(dat) {
  g_fn <- g_fn_maker(dat)
  n <- nrow(dat)
  h_dim <- ncol(cbind(1, h_alpha_fn(dat)))
  function(alpha) {
    g_mat <- g_fn(alpha)
    G <- colMeans(g_mat)
    W <- tryCatch(solve(t(g_mat) %*% g_mat / n), error = function(e) diag(h_dim))
    as.numeric(t(G) %*% W %*% G)
  }
}

# Get NR cond=1e2 solution (reuse code from earlier)
Gamma_fn_maker <- function(dat) {
  design_mat <- model.matrix(design_formula, dat)
  r_vec <- as.vector(dat$r)
  h_x <- cbind(1, h_alpha_fn(dat))
  h_dim <- ncol(h_x)
  alpha_dim <- ncol(design_mat)
  function(alpha) {
    eta <- pmin(pmax(as.vector(design_mat %*% alpha), -20), 20)
    pi_vec <- compute_pi(eta)
    common <- r_vec * (1 - pi_vec) / pi_vec
    Gamma_mat <- matrix(0, h_dim, alpha_dim)
    for (j in 1:alpha_dim) Gamma_mat[, j] <- colMeans(common * design_mat[, j] * h_x)
    Gamma_mat
  }
}

run_nr_cond <- function(dat, cond_thresh) {
  g_fn <- g_fn_maker(dat)
  G_fn <- function(alpha) colMeans(g_fn(alpha))
  Gamma_fn <- Gamma_fn_maker(dat)
  n <- nrow(dat)
  h_dim <- ncol(cbind(1, h_alpha_fn(dat)))
  alpha_dim <- ncol(model.matrix(design_formula, dat))

  nr_inner <- function(start, W_hat) {
    alpha <- start
    for (k in 1:500) {
      Gamma_hat <- Gamma_fn(alpha)
      G_vec <- matrix(G_fn(alpha), ncol=1)
      g_vec <- 2 * as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))
      if (max(abs(g_vec)) < 1e-8) break
      H_gn <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
      direction <- tryCatch(solve(H_gn, -g_vec), error = function(e) -g_vec)
      sn <- sqrt(sum(direction^2))
      if (sn > 2.0) direction <- direction * (2.0 / sn)
      obj_cur <- as.numeric(crossprod(G_vec, W_hat %*% G_vec))
      step <- 1.0; accepted <- FALSE
      for (ls in 1:30) {
        cand <- alpha + step * direction
        G_cand <- matrix(G_fn(cand), ncol=1)
        obj_new <- as.numeric(crossprod(G_cand, W_hat %*% G_cand))
        if (is.finite(obj_new) && obj_new < obj_cur - 1e-4 * step * sum(g_vec * direction)) {
          Gamma_c <- Gamma_fn(cand)
          H_c <- crossprod(Gamma_c, W_hat %*% Gamma_c)
          if (rcond(H_c) > 1/cond_thresh) { accepted <- TRUE; break }
        }
        step <- step * 0.5
      }
      if (accepted) { alpha <- cand } else { break }
    }
    alpha
  }

  alpha <- nr_inner(rep(0, alpha_dim), diag(h_dim))
  best_grad <- Inf; best_alpha <- alpha
  for (t in 1:200) {
    g_mat <- g_fn(alpha)
    W_hat <- tryCatch(solve(t(g_mat) %*% g_mat / n), error = function(e) diag(h_dim))
    G_vec <- G_fn(alpha)
    grad_norm <- max(abs(2 * as.vector(t(Gamma_fn(alpha)) %*% W_hat %*% G_vec)))
    if (grad_norm < 1e-6) break
    if (grad_norm < best_grad) { best_grad <- grad_norm; best_alpha <- alpha }
    alpha <- nr_inner(alpha, W_hat)
  }
  best_alpha
}

cat("=== Full dQ/dα comparison: L-BFGS-B vs NR cond=1e8 vs NR cond=1e2 ===\n\n")

W_sm <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
ps_spec_lb <- list(formula.list = list(design_formula), h_alpha.list = list(h_alpha_fn),
                   inv_link = ps_spec_base$inv_link, outcome = ps_spec_base$outcome)

for (rep_i in c(1, 7, 13)) {
  dat <- all_data[((rep_i-1)*nn + 1):(rep_i*nn), ]
  Q_cue <- Q_cue_maker(dat)
  g_fn <- g_fn_maker(dat)
  Gamma_fn <- Gamma_fn_maker(dat)
  n <- nn; h_dim <- ncol(cbind(1, h_alpha_fn(dat)))

  # L-BFGS-B
  ebmr <- EPS$new("y", ps_spec_lb, dat, W_sm)
  alpha_lb <- ebmr$ps_fit.list[[1]]$coefficients

  # NR cond=1e8
  alpha_nr8 <- run_nr_cond(dat, 1e8)

  # NR cond=1e2
  alpha_nr2 <- run_nr_cond(dat, 1e2)

  cat(sprintf("--- Rep %d ---\n", rep_i))
  for (label in c("L-BFGS-B", "NR 1e8", "NR 1e2")) {
    alpha <- switch(label, "L-BFGS-B" = alpha_lb, "NR 1e8" = alpha_nr8, "NR 1e2" = alpha_nr2)

    g_mat <- g_fn(alpha)
    G_vec <- colMeans(g_mat)
    W_hat <- solve(t(g_mat) %*% g_mat / n)
    iter_grad <- 2 * as.vector(t(Gamma_fn(alpha)) %*% W_hat %*% G_vec)
    full_grad <- numDeriv::grad(Q_cue, alpha)

    cat(sprintf("  %-8s: max|2Γ'WG|=%.2e  max|full dQ/dα|=%.2e  Q=%.4e  mu=%.2f\n",
        label, max(abs(iter_grad)), max(abs(full_grad)), Q_cue(alpha),
        mean(dat$r * dat$y / compute_pi(pmin(pmax(as.vector(model.matrix(design_formula, dat) %*% alpha), -20), 20)))))
  }
  cat("\n")
}
