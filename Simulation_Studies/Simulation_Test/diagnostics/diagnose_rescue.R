setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
  library(EBMRalgorithmFast4)
})

W_func  <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

ps_spec <- get_ps_spec("9-alt1")
all_data <- readRDS("Simulation_Data/Setting3.B1_n2000_replicate1000.RDS")
n_val <- 2000

at_boundary <- function(alpha, lo, hi, tol = 1e-4) {
  any(alpha <= lo + tol | alpha >= hi - tol)
}

# Run model 2 with verbose rescue diagnostics
# Reimplement the rescue logic inline so we can print what's happening
run_with_rescue_diag <- function(i) {
  set.seed(12345 + i)
  dat  <- all_data[((i-1)*n_val + 1):(i*n_val), ]

  # Build model 2 components manually
  ps_sub <- list(
    formula.list = ps_spec[["formula.list"]][2],
    h_alpha.list = ps_spec[["h_alpha.list"]][2],
    inv_link     = ps_spec[["inv_link"]],
    outcome      = ps_spec[["outcome"]]
  )

  ebmr <- EBMRAlgorithmFast4$new("y", ps_sub, dat, W_func)
  # Access the internals via the environment trick
  # We need to call WangShaoKim2014 directly and trace the rescue
  # Instead: run normally and report what happened
  h_nu_fn <- function(d) cbind(u1=d$u1, u2=d$u2, z1=d$z1, z2=d$z2, u1_u2=d$u1*d$u2)
  res <- ebmr$EBMR_IPW(h_nu=h_nu_fn, type="HT", se.fit=TRUE)
  fit <- ebmr$ps_fit.list[[1]]$gmm_fit

  hpd <- isTRUE(fit$opt$hessian_pd)
  alpha <- fit$estimates
  obj <- fit$opt$objective

  list(hpd=hpd, alpha=alpha, obj=obj, mu=res$mu_ipw, se=res$se_ipw)
}

# Find non-PD cases
cat("Scanning for non-PD cases...\n")
results <- list()
for (i in 1:50) {
  r <- run_with_rescue_diag(i)
  results[[i]] <- r
  if (!r$hpd) cat(sprintf("rep %d: hess_pd=FALSE  alpha=[%s]  obj=%.5f  mu=%.3f\n",
      i, paste(round(r$alpha, 3), collapse=", "), r$obj, r$mu))
}

# Now re-run non-PD cases with manual rescue diagnostics
cat("\n=== Manual rescue diagnostics for non-PD cases ===\n")
non_pd_reps <- which(sapply(results, function(r) !r$hpd))
cat("Non-PD reps:", non_pd_reps, "\n\n")

for (i in non_pd_reps[1:min(3, length(non_pd_reps))]) {
  cat(sprintf("--- Rep %d ---\n", i))
  set.seed(12345 + i)
  dat <- all_data[((i-1)*n_val + 1):(i*n_val), ]

  # Rebuild the MNAR model 2 components
  formula_2 <- ps_spec[["formula.list"]][[2]]
  h_alpha_2 <- ps_spec[["h_alpha.list"]][[2]]
  inv_link   <- ps_spec[["inv_link"]]

  # Build design matrix
  r_col <- dat$r
  parsed_r <- "r"
  dm <- model.matrix(formula_2, data = dat)
  design_matrix <- dm
  alpha_dim <- ncol(design_matrix)

  h_x_raw <- if (is.character(h_alpha_2)) dat[h_alpha_2] else h_alpha_2(dat)
  h_x <- cbind(1, as.matrix(h_x_raw))
  h_dim <- ncol(h_x)

  r_vec <- as.vector(r_col)
  eta_max <- 20
  compute_pi <- function(eta) {
    eta_c <- pmin(pmax(eta, -eta_max), eta_max)
    1 / (1 + exp(eta_c))  # logistic_complement
  }

  Phi_alpha <- function(param) {
    eta <- as.vector(design_matrix %*% param)
    pi_vec <- compute_pi(eta)
    (r_vec / pi_vec - 1) * h_x
  }

  W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

  Q_full <- function(alpha) {
    gm <- Phi_alpha(alpha)
    G  <- matrix(colMeans(gm), h_dim, 1)
    W  <- tryCatch(W_fn(gm), error = function(e) diag(h_dim))
    as.numeric(crossprod(G, W %*% G))
  }

  n <- nrow(dat)
  init <- rep(0, alpha_dim)

  # Run unconstrained first
  state <- new.env(parent=emptyenv())
  state$W.hat <- diag(h_dim)
  state$cache_g <- NULL; state$cache_G <- NULL; state$cache_param <- NULL

  cache_fn <- function(param) {
    if (!is.null(state$cache_param) && identical(param, state$cache_param)) return()
    state$cache_param <- param
    state$cache_g <- Phi_alpha(param)
    state$cache_G <- matrix(colMeans(state$cache_g), h_dim, 1)
  }
  obj_fn <- function(param) {
    cache_fn(param)
    v <- as.numeric(crossprod(state$cache_G, state$W.hat %*% state$cache_G))
    if (!is.finite(v)) 1e8 else v
  }

  # Step 1: identity W
  opt0 <- optim(init, obj_fn, method="L-BFGS-B", control=list(maxit=1000))
  alpha0 <- opt0$par
  cat(sprintf("  Unconstrained init opt: alpha=[%s] obj=%.5f\n",
      paste(round(alpha0, 3), collapse=", "), Q_full(alpha0)))

  # Iterate
  for (t in 1:50) {
    cache_fn(alpha0)
    state$W.hat <- tryCatch(W_fn(state$cache_g), error=function(e) diag(h_dim))
    state$cache_param <- NULL
    opt_t <- optim(alpha0, obj_fn, method="L-BFGS-B", control=list(maxit=1000))
    alpha0 <- opt_t$par
  }
  cat(sprintf("  Unconstrained final:    alpha=[%s] obj=%.5f\n",
      paste(round(alpha0, 3), collapse=", "), Q_full(alpha0)))

  # Now try rescue bounds
  for (bound_val in c(10, 5, 3, 2)) {
    lo <- rep(-bound_val, alpha_dim)
    hi <- rep( bound_val, alpha_dim)
    init_c <- pmin(pmax(init, lo), hi)

    state$W.hat <- diag(h_dim); state$cache_param <- NULL
    opt0b <- optim(init_c, obj_fn, method="L-BFGS-B", lower=lo, upper=hi, control=list(maxit=1000))
    alpha_b <- opt0b$par

    for (t in 1:50) {
      cache_fn(alpha_b)
      state$W.hat <- tryCatch(W_fn(state$cache_g), error=function(e) diag(h_dim))
      state$cache_param <- NULL
      opt_t <- optim(alpha_b, obj_fn, method="L-BFGS-B", lower=lo, upper=hi, control=list(maxit=1000))
      alpha_b <- opt_t$par
    }

    is_bound <- at_boundary(alpha_b, lo, hi)
    cat(sprintf("  Bound ±%d: alpha=[%s] obj=%.5f  at_boundary=%s\n",
        bound_val, paste(round(alpha_b, 3), collapse=", "),
        Q_full(alpha_b), is_bound))
  }
  cat("\n")
}
