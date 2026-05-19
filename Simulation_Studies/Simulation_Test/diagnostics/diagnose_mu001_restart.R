## Test M3 with multiple random restarts to find a better solution
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

ps_spec_base <- get_ps_spec("9-alt1")
h_alpha_fn <- function(dat) cbind(
  u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2,
  u1u2 = dat$u1*dat$u2, z1z2 = dat$z1*dat$z2,
  u1z1 = dat$u1*dat$z1, z1u2 = dat$z1*dat$u2,
  u1z2 = dat$u1*dat$z2, u2z2 = dat$u2*dat$z2)

data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]]
all_data <- NULL
for (fi in seq_along(data_file)) {
  if (file.exists(data_file[[fi]])) {
    test <- readRDS(data_file[[fi]])
    if (nrow(test) / 1000 == 2000) { all_data <- test; break }
  }
}

dat <- all_data[1:2000, ]

# Setup
design_mat <- model.matrix(ps_spec_base$formula.list[[3]], dat)
r_vec <- as.vector(dat$r)
n <- 2000
h_x <- cbind(1, h_alpha_fn(dat))
h_dim <- ncol(h_x)
alpha_dim <- ncol(design_mat)

compute_pi <- function(eta) 1/(1+exp(eta))

g_fn <- function(alpha) {
  eta <- as.vector(design_mat %*% alpha)
  eta <- pmin(pmax(eta, -20), 20)
  pi_vec <- compute_pi(eta)
  as.vector(r_vec/pi_vec - 1) * h_x
}
G_fn <- function(alpha) colMeans(g_fn(alpha))

# Run iterative CUE from a given starting point
run_cue <- function(init, max_iter = 50) {
  alpha <- init
  best_grad <- Inf
  best_alpha <- alpha
  for (t in 1:max_iter) {
    g_mat <- g_fn(alpha)
    W_hat <- tryCatch(solve(var(g_mat)), error = function(e) diag(h_dim))
    opt <- optim(alpha, function(a) {
      G <- matrix(G_fn(a), ncol=1)
      as.numeric(t(G) %*% W_hat %*% G)
    }, method = "L-BFGS-B", control = list(maxit = 1000))
    alpha <- opt$par
    G_new <- G_fn(alpha)
    grad_norm <- max(abs(G_new))
    if (grad_norm < best_grad) {
      best_grad <- grad_norm
      best_alpha <- alpha
    }
    if (grad_norm < 1e-6) break
  }
  list(alpha = best_alpha, grad_norm = best_grad,
       G = G_fn(best_alpha), obj = opt$value)
}

cat("=== Testing multiple restarts for M3, W=solve(var(g)) ===\n\n")
cat(sprintf("%4s | %12s | %12s | %s\n", "#", "grad_norm", "obj", "alpha"))
cat(strrep("-", 80), "\n")

# Try from zeros
res0 <- run_cue(rep(0, alpha_dim))
cat(sprintf("%4s | %12.4e | %12.4e | %s\n", "zero",
    res0$grad_norm, res0$obj, paste(round(res0$alpha, 3), collapse=", ")))

# Try random restarts
set.seed(42)
for (r in 1:20) {
  init <- rnorm(alpha_dim, sd = 2)
  res <- run_cue(init)
  cat(sprintf("%4d | %12.4e | %12.4e | %s\n", r,
      res$grad_norm, res$obj, paste(round(res$alpha, 3), collapse=", ")))
}
