## Test: Does multistart for nu change the ensemble results?
## Check if the GMM objective for nu has multiple local minima
## and if the current single-start (0.5, 0.5) is missing the global min
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

n_val <- 2000
ps_spec <- get_ps_spec("9-alt1")
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

h_nu_fn <- function(dat) cbind(
  u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2,
  u1u2 = dat$u1*dat$u2, u1z1 = dat$u1*dat$z1, u1z2 = dat$u1*dat$z2,
  u2z1 = dat$u2*dat$z1, u2z2 = dat$u2*dat$z2, z1z2 = dat$z1*dat$z2,
  u1u2z1 = dat$u1*dat$u2*dat$z1, u1u2z2 = dat$u1*dat$u2*dat$z2,
  u1z1z2 = dat$u1*dat$z1*dat$z2, u2z1z2 = dat$u2*dat$z1*dat$z2
)

data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
mu_true <- get_mu_true("setting4")

cat("=== Diagnosing ensemble nu optimization across reps ===\n\n")

# Check 20 reps in detail
for (rep_i in 1:20) {
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]

  ebmr <- EBMRAlgorithmFast4$new("y", ps_spec, dat, W_fn)

  # Get PS matrices for M2 and M3
  ps2 <- ebmr$ps_fit.list[[2]]$fitted.values
  ps3 <- ebmr$ps_fit.list[[3]]$fitted.values
  ps_mat <- cbind(ps2, ps3)
  r_vec <- dat$r

  # h_nu with auto-intercept (as ensemble does)
  h_x2 <- h_nu_fn(dat)
  h_x <- cbind(1, h_x2)
  n <- nrow(dat)
  h_dim <- ncol(h_x)

  # GMM objective for nu
  gmm_obj <- function(nu) {
    ps_nu <- as.vector(ps_mat %*% nu)
    g_mat <- as.vector(r_vec / ps_nu - 1) * h_x
    G_vec <- colMeans(g_mat)
    W_hat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(h_dim))
    as.numeric(t(G_vec) %*% W_hat %*% G_vec)
  }

  # Try multiple starting points
  starts <- list(
    c(0.5, 0.5),
    c(0.9, 0.1),
    c(0.1, 0.9),
    c(1.0, 0.0001),
    c(0.0001, 1.0),
    c(0.3, 0.7),
    c(0.7, 0.3),
    c(1, 1),
    c(2, 1),
    c(1, 2)
  )

  best_obj <- Inf
  best_nu <- NULL
  all_results <- list()

  for (s in seq_along(starts)) {
    nu_init <- starts[[s]]
    tryCatch({
      # Run EBMR_IPW with this nu_init
      result <- ebmr$EBMR_IPW(
        h_nu = h_nu_fn,
        model_indices = c(2, 3),
        nu_init = nu_init,
        se.fit = FALSE,
        type = "HT"
      )
      w <- result$w.hat
      obj <- gmm_obj(result$nu.hat)
      all_results[[s]] <- list(nu = result$nu.hat, w = w, obj = obj, mu = result$mu_ipw)
      if (obj < best_obj) {
        best_obj <- obj
        best_nu <- result$nu.hat
      }
    }, error = function(e) {
      all_results[[s]] <<- list(nu = NA, w = NA, obj = Inf, mu = NA)
    })
  }

  # Summarize
  objs <- sapply(all_results, function(r) r$obj)
  mus <- sapply(all_results, function(r) r$mu)
  w1s <- sapply(all_results, function(r) if(length(r$w) >= 1) r$w[1] else NA)

  n_distinct <- length(unique(round(mus[is.finite(mus)], 4)))
  obj_range <- diff(range(objs[is.finite(objs)]))

  cat(sprintf("Rep %2d: %d distinct solutions, obj range=%.2e\n", rep_i, n_distinct, obj_range))
  cat(sprintf("  w1 values: %s\n", paste(round(w1s[is.finite(w1s)], 3), collapse=", ")))
  cat(sprintf("  mu values: %s\n", paste(round(mus[is.finite(mus)], 4), collapse=", ")))
  cat(sprintf("  obj values: %s\n", paste(format(objs[is.finite(objs)], digits=4, scientific=TRUE), collapse=", ")))
  cat("\n")
}
