setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

n_val <- 2000
ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
mu_true <- get_mu_true("setting4")
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

model_idx <- 3

# First, collect typical normal-rep solutions as reference
cat("=== Collecting normal rep solutions as reference ===\n")
normal_alphas <- list()
for (rep_i in c(1, 2, 3, 5, 10, 20, 50)) {
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL),
    optimizer = "constrained_nr"
  )
  ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
  gmm <- ebmr$ps_fit.list[[1]]$gmm_fit
  if (gmm$opt$solution_cond < 1e8) {
    normal_alphas[[length(normal_alphas)+1]] <- gmm$estimates
    cat(sprintf("  Rep %d: alpha=%s, cond=%.1e\n", rep_i,
        paste(round(gmm$estimates, 3), collapse=", "), gmm$opt$solution_cond))
  }
}

# Known degenerate reps
degen_reps <- c(35, 38, 68, 109, 133, 155, 161, 166)

compute_pi_fn <- function(eta) plogis(-eta)

for (rep_i in degen_reps) {
  cat(sprintf("\n========== Rep %d ==========\n", rep_i))
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
  r_vec <- dat[["r"]]; n <- n_val

  # Get design_mat and h_x from package
  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL),
    optimizer = "L-BFGS-B"
  )
  ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
  design_mat <- ebmr$ps_fit.list[[1]]$design_matrix
  h_x <- ebmr$ps_fit.list[[1]]$h_x
  k <- ncol(design_mat)
  h_dim <- ncol(h_x)

  g_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    (r_vec / pi_v - 1) * h_x
  }
  Gamma_fn <- function(alpha) {
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    cf <- r_vec * (1 - pi_v) / pi_v
    crossprod(h_x * cf, design_mat) / n
  }
  obj_fn <- function(alpha) {
    g_mat <- g_fn(alpha)
    G <- colMeans(g_mat)
    W <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(h_dim))
    as.numeric(t(G) %*% W %*% G)
  }
  grad_fn <- function(alpha) {
    g_mat <- g_fn(alpha)
    G <- matrix(colMeans(g_mat), h_dim, 1)
    W <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(h_dim))
    Gamma <- Gamma_fn(alpha)
    max(abs(2 * as.vector(crossprod(Gamma, W %*% G))))
  }
  cond_fn <- function(alpha) {
    g_mat <- g_fn(alpha)
    W <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(h_dim))
    Gamma <- Gamma_fn(alpha)
    H <- crossprod(Gamma, W %*% Gamma)
    eig <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    max(eig) / max(min(eig), 1e-15)
  }

  # L-BFGS-B iterative GMM from a starting point (no cond constraint)
  run_lbfgsb <- function(start) {
    st <- new.env(parent = emptyenv())
    st$W <- NULL; st$cp <- NULL; st$cg <- NULL; st$cG <- NULL
    cache <- function(p) {
      if (!is.null(st$cp) && identical(p, st$cp)) return()
      st$cp <- p; st$cg <- g_fn(p); st$cG <- matrix(colMeans(st$cg), h_dim, 1)
    }
    of <- function(p) { cache(p); v <- as.numeric(crossprod(st$cG, st$W %*% st$cG)); if(!is.finite(v)) 1e8 else v }
    gf <- function(p) { cache(p); 2 * as.vector(crossprod(Gamma_fn(p), st$W %*% st$cG)) }

    est <- start
    st$W <- diag(h_dim)
    opt0 <- tryCatch(optim(est, of, gr = gf, method = "L-BFGS-B", control = list(maxit = 500)),
                     error = function(e) list(par = est))
    est <- opt0$par
    for (t in 1:50) {
      cache(est)
      st$W <- tryCatch(solve(crossprod(st$cg) / n), error = function(e) diag(h_dim))
      fg <- max(abs(gf(est)))
      if (fg < 1e-6) break
      st$cp <- NULL
      opt_t <- tryCatch(optim(est, of, gr = gf, method = "L-BFGS-B", control = list(maxit = 500)),
                        error = function(e) list(par = est))
      est <- opt_t$par
    }
    g_f <- g_fn(est); G_f <- matrix(colMeans(g_f), h_dim, 1)
    W_f <- tryCatch(solve(crossprod(g_f) / n), error = function(e) diag(h_dim))
    gn <- max(abs(2 * as.vector(crossprod(Gamma_fn(est), W_f %*% G_f))))
    list(est = est, grad = gn)
  }

  # ========== Thorough search ==========
  results <- data.frame(method=character(), grad=numeric(), cond=numeric(),
                        obj=numeric(), alpha2=numeric(), mu=numeric(),
                        stringsAsFactors=FALSE)

  add_result <- function(method, alpha) {
    gn <- tryCatch(grad_fn(alpha), error = function(e) NA)
    cn <- tryCatch(cond_fn(alpha), error = function(e) NA)
    ob <- tryCatch(obj_fn(alpha), error = function(e) NA)
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    mu <- mean(r_vec * dat[["y"]] / pi_v)
    results[nrow(results)+1, ] <<- list(method, gn, cn, ob, alpha[2], mu)
  }

  # 1. Start from each normal rep's solution
  for (i in seq_along(normal_alphas)) {
    res <- run_lbfgsb(normal_alphas[[i]])
    add_result(sprintf("normal_init_%d", i), res$est)
  }

  # 2. Grid over alpha[2] from -1 to -8, fix others at normal-rep average
  avg_alpha <- Reduce(`+`, normal_alphas) / length(normal_alphas)
  for (a2 in seq(-1, -8, by = -1)) {
    init <- avg_alpha; init[2] <- a2
    res <- run_lbfgsb(init)
    add_result(sprintf("grid_a2=%.0f", a2), res$est)
  }

  # 3. 100 random starts from N(0, 1)
  set.seed(rep_i)
  n_found_good <- 0
  for (trial in 1:100) {
    init <- rnorm(k) * 1.0
    res <- tryCatch(run_lbfgsb(init), error = function(e) list(est = init, grad = NA))
    if (!is.na(res$grad)) {
      cn <- tryCatch(cond_fn(res$est), error = function(e) NA)
      if (!is.na(cn) && cn < 1e4 && res$grad < 0.01) {
        n_found_good <- n_found_good + 1
        add_result(sprintf("random_%d", trial), res$est)
      }
    }
  }
  # Also record total random that converged with low cond
  cat(sprintf("  Random starts with cond<1e4 & grad<0.01: %d/100\n", n_found_good))

  # 4. Two-step GMM: use W from identity, then one W update
  for (init_alpha in c(list(rep(0, k)), normal_alphas[1:2])) {
    g0 <- g_fn(init_alpha)
    W_I <- diag(h_dim)
    # Step 1: optimize with W = I
    st <- new.env(parent = emptyenv())
    st$W <- W_I; st$cp <- NULL; st$cg <- NULL; st$cG <- NULL
    cache2 <- function(p) {
      if (!is.null(st$cp) && identical(p, st$cp)) return()
      st$cp <- p; st$cg <- g_fn(p); st$cG <- matrix(colMeans(st$cg), h_dim, 1)
    }
    of2 <- function(p) { cache2(p); v <- as.numeric(crossprod(st$cG, st$W %*% st$cG)); if(!is.finite(v)) 1e8 else v }
    gf2 <- function(p) { cache2(p); 2 * as.vector(crossprod(Gamma_fn(p), st$W %*% st$cG)) }
    opt1 <- tryCatch(optim(init_alpha, of2, gr = gf2, method = "L-BFGS-B", control = list(maxit = 1000)),
                     error = function(e) list(par = init_alpha))
    est1 <- opt1$par
    # Step 2: one W update
    cache2(est1)
    st$W <- tryCatch(solve(crossprod(st$cg) / n), error = function(e) diag(h_dim))
    st$cp <- NULL
    opt2 <- tryCatch(optim(est1, of2, gr = gf2, method = "L-BFGS-B", control = list(maxit = 1000)),
                     error = function(e) list(par = est1))
    add_result("twostep", opt2$par)
  }

  # Report best results
  results <- results[order(results$cond), ]
  cat(sprintf("  %20s  %10s  %10s  %10s  %10s  %10s\n",
      "method", "grad", "cond", "obj", "alpha2", "mu"))
  n_show <- min(nrow(results), 15)
  for (i in 1:n_show) {
    r <- results[i, ]
    cat(sprintf("  %20s  %10.2e  %10.2e  %10.6f  %10.4f  %10.4f\n",
        r$method, r$grad, r$cond, r$obj, r$alpha2, r$mu))
  }
}
