## Test different optimizers for M3 CUE with W=solve(var(g))
## The issue: iterative CUE oscillates in a 2-cycle, never converges
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EPS", quiet = TRUE)

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

# CUE objective: Q(alpha) = G(alpha)' W(alpha) G(alpha)
# where W(alpha) = solve(var(g(alpha)))
cue_obj <- function(alpha) {
  g_mat <- g_fn(alpha)
  G <- colMeans(g_mat)
  W <- tryCatch(solve(var(g_mat)), error = function(e) diag(h_dim))
  as.numeric(t(G) %*% W %*% G)
}

# CUE gradient (numerical)
cue_grad_num <- function(alpha) {
  numDeriv::grad(cue_obj, alpha)
}

cat("=== Testing optimizers for M3 CUE with W=solve(var(g)), Rep 1 ===\n\n")

# 1. Direct CUE optimization with optim (Nelder-Mead)
cat("--- 1. Nelder-Mead on CUE objective ---\n")
res_nm <- optim(rep(0, alpha_dim), cue_obj, method = "Nelder-Mead",
                control = list(maxit = 10000))
cat(sprintf("   obj=%.6e, convergence=%d, alpha=(%s)\n",
    res_nm$value, res_nm$convergence, paste(round(res_nm$par, 4), collapse=", ")))
G_nm <- G_fn(res_nm$par)
cat(sprintf("   max|G|=%.4e\n\n", max(abs(G_nm))))

# 2. BFGS on CUE objective (not L-BFGS-B, no bounds)
cat("--- 2. BFGS on CUE objective ---\n")
res_bfgs <- optim(rep(0, alpha_dim), cue_obj, method = "BFGS",
                  control = list(maxit = 5000))
cat(sprintf("   obj=%.6e, convergence=%d, alpha=(%s)\n",
    res_bfgs$value, res_bfgs$convergence, paste(round(res_bfgs$par, 4), collapse=", ")))
G_bfgs <- G_fn(res_bfgs$par)
cat(sprintf("   max|G|=%.4e\n\n", max(abs(G_bfgs))))

# 3. L-BFGS-B on CUE objective directly (not iterative)
cat("--- 3. L-BFGS-B on CUE objective directly ---\n")
res_lbfgsb <- optim(rep(0, alpha_dim), cue_obj, method = "L-BFGS-B",
                    control = list(maxit = 5000))
cat(sprintf("   obj=%.6e, convergence=%d, alpha=(%s)\n",
    res_lbfgsb$value, res_lbfgsb$convergence, paste(round(res_lbfgsb$par, 4), collapse=", ")))
G_lbfgsb <- G_fn(res_lbfgsb$par)
cat(sprintf("   max|G|=%.4e\n\n", max(abs(G_lbfgsb))))

# 4. nlminb
cat("--- 4. nlminb on CUE objective ---\n")
res_nlminb <- nlminb(rep(0, alpha_dim), cue_obj,
                     control = list(iter.max = 5000, eval.max = 10000))
cat(sprintf("   obj=%.6e, convergence=%d, message=%s\n",
    res_nlminb$objective, res_nlminb$convergence, res_nlminb$message))
cat(sprintf("   alpha=(%s)\n", paste(round(res_nlminb$par, 4), collapse=", ")))
G_nlminb <- G_fn(res_nlminb$par)
cat(sprintf("   max|G|=%.4e\n\n", max(abs(G_nlminb))))

# 5. Nelder-Mead with numerical gradient (CUE)
cat("--- 5. CG (conjugate gradient) on CUE objective ---\n")
res_cg <- optim(rep(0, alpha_dim), cue_obj, cue_grad_num, method = "CG",
                control = list(maxit = 5000))
cat(sprintf("   obj=%.6e, convergence=%d, alpha=(%s)\n",
    res_cg$value, res_cg$convergence, paste(round(res_cg$par, 4), collapse=", ")))
G_cg <- G_fn(res_cg$par)
cat(sprintf("   max|G|=%.4e\n\n", max(abs(G_cg))))

# 6. Try solving G(alpha) = 0 directly with nleqslv (if available)
cat("--- 6. nleqslv (solving G(alpha)=0 directly) ---\n")
if (requireNamespace("nleqslv", quietly = TRUE)) {
  # Need h_dim == alpha_dim for exactly identified; we're overidentified
  # Use min ||G||^2 approach instead
  cat("   Overidentified (h_dim=%d > alpha_dim=%d), skipping nleqslv\n\n", h_dim, alpha_dim)
} else {
  cat("   nleqslv not available\n\n")
}

# 7. Try DEoptim (differential evolution) - global optimizer
cat("--- 7. DEoptim (global optimizer) ---\n")
if (requireNamespace("DEoptim", quietly = TRUE)) {
  library(DEoptim)
  res_de <- DEoptim(cue_obj, lower = rep(-10, alpha_dim), upper = rep(10, alpha_dim),
                    control = DEoptim.control(itermax = 500, trace = FALSE))
  cat(sprintf("   obj=%.6e, alpha=(%s)\n",
      res_de$optim$bestval, paste(round(res_de$optim$bestmem, 4), collapse=", ")))
  G_de <- G_fn(res_de$optim$bestmem)
  cat(sprintf("   max|G|=%.4e\n\n", max(abs(G_de))))
} else {
  cat("   DEoptim not available\n\n")
}

# 8. Try GenSA (generalized simulated annealing)
cat("--- 8. GenSA (simulated annealing) ---\n")
if (requireNamespace("GenSA", quietly = TRUE)) {
  library(GenSA)
  res_sa <- GenSA(rep(0, alpha_dim), cue_obj,
                  lower = rep(-10, alpha_dim), upper = rep(10, alpha_dim),
                  control = list(maxit = 5000, verbose = FALSE))
  cat(sprintf("   obj=%.6e, alpha=(%s)\n",
      res_sa$value, paste(round(res_sa$par, 4), collapse=", ")))
  G_sa <- G_fn(res_sa$par)
  cat(sprintf("   max|G|=%.4e\n\n", max(abs(G_sa))))
} else {
  cat("   GenSA not available\n\n")
}

# 9. Damped iterative CUE: average W across iterations
cat("--- 9. Damped iterative CUE (W averaging) ---\n")
alpha <- rep(0, alpha_dim)
W_prev <- diag(h_dim)
for (t in 1:100) {
  g_mat <- g_fn(alpha)
  W_new <- tryCatch(solve(var(g_mat)), error = function(e) diag(h_dim))
  # Average with previous W (damping)
  W_hat <- 0.5 * W_prev + 0.5 * W_new
  W_prev <- W_hat
  opt <- optim(alpha, function(a) {
    G <- matrix(G_fn(a), ncol=1)
    as.numeric(t(G) %*% W_hat %*% G)
  }, method = "L-BFGS-B", control = list(maxit = 1000))
  alpha <- opt$par
}
G_damp <- G_fn(alpha)
g_damp <- g_fn(alpha)
obj_damp <- as.numeric(t(G_damp) %*% solve(var(g_damp)) %*% G_damp)
cat(sprintf("   obj=%.6e, max|G|=%.4e, alpha=(%s)\n\n",
    obj_damp, max(abs(G_damp)), paste(round(alpha, 4), collapse=", ")))

# 10. Fixed W from large sample (use all 2000 obs to form stable W, then optimize once)
cat("--- 10. One-shot: W from initial alpha=0, then optimize once ---\n")
g_init <- g_fn(rep(0, alpha_dim))
W_fixed <- solve(var(g_init))
opt_fixed <- optim(rep(0, alpha_dim), function(a) {
  G <- matrix(G_fn(a), ncol=1)
  as.numeric(t(G) %*% W_fixed %*% G)
}, method = "L-BFGS-B", control = list(maxit = 5000))
G_fixed <- G_fn(opt_fixed$par)
g_fixed <- g_fn(opt_fixed$par)
obj_fixed <- as.numeric(t(G_fixed) %*% solve(var(g_fixed)) %*% G_fixed)
cat(sprintf("   obj(CUE)=%.6e, max|G|=%.4e, alpha=(%s)\n\n",
    obj_fixed, max(abs(G_fixed)), paste(round(opt_fixed$par, 4), collapse=", ")))
