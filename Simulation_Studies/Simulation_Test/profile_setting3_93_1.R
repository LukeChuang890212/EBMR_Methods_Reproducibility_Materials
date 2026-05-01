setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
source("Basic_setup.r"); source("Data_Generation.r")
source("config/scenarios.R"); source("Simulation.r")
library(EBMRalgorithmFast4)

# Setting: setting3-miss50-scenario9-3_1  (model 1 only, misspecified)
ps_spec <- get_ps_spec("9-alt1")
ps_spec_1 <- list(
  formula.list = ps_spec[["formula.list"]][1],
  h_alpha.list = ps_spec[["h_alpha.list"]][1],
  inv_link     = ps_spec[["inv_link"]],
  outcome      = ps_spec[["outcome"]]
)
W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

data_file <- "Simulation_Data/Setting3.B1_n2000_replicate1000.RDS"
if (!file.exists(data_file)) {
  cat("Data file not found:", data_file, "\n")
  cat("Available data files:\n")
  print(list.files("Simulation_Data", pattern="Setting3"))
  stop("Missing data file")
}
all_data <- readRDS(data_file)
n_val <- 2000

cat("=== Profiling setting3-miss50-scenario9-3_1, n=2000 ===\n\n")

# ─── 1. Single-replicate timing breakdown ────────────────────────────────────
n_probe <- 20
cat("--- 1. Per-replicate timing (", n_probe, "reps) ---\n")

times_total   <- numeric(n_probe)
hess_pd_vals  <- logical(n_probe)
iter_counts   <- integer(n_probe)
grad_norms    <- numeric(n_probe)
alpha_norms   <- numeric(n_probe)
had_rescue    <- logical(n_probe)

for (i in 1:n_probe) {
  set.seed(12345 + i)
  dat <- all_data[((i-1)*n_val + 1):(i*n_val), ]

  t0 <- proc.time()["elapsed"]
  ebmr   <- EBMRAlgorithmFast4$new("y", ps_spec_1, dat, W_func)
  result <- ebmr$EBMR_IPW(
    h_nu = function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1_u2=dat$u1*dat$u2),
    type = "HT", se.fit = TRUE
  )
  times_total[i] <- proc.time()["elapsed"] - t0

  ps_fit  <- ebmr$ps_fit.list[[1]]
  gmm_opt <- ps_fit$gmm_fit$opt
  hess_pd_vals[i] <- isTRUE(gmm_opt$hessian_pd)
  iter_counts[i]  <- gmm_opt$iterations
  grad_norms[i]   <- ifelse(is.na(gmm_opt$final_grad_norm), Inf, gmm_opt$final_grad_norm)
  alpha_norms[i]  <- sqrt(sum(ps_fit$coefficients^2))
}

cat(sprintf("  Time (sec): mean=%.2f  min=%.2f  max=%.2f  total=%.1f\n",
    mean(times_total), min(times_total), max(times_total), sum(times_total)))
cat(sprintf("  hessian_pd=TRUE: %d/%d\n", sum(hess_pd_vals), n_probe))
cat(sprintf("  outer iterations: mean=%.0f  max=%d\n",
    mean(iter_counts), max(iter_counts)))
cat(sprintf("  final_grad_norm: mean=%.2e  max=%.2e\n",
    mean(grad_norms), max(grad_norms)))
cat(sprintf("  alpha_norm: mean=%.2f  max=%.2f\n",
    mean(alpha_norms), max(alpha_norms)))

cat("\n  Per-rep details:\n")
cat(sprintf("  %4s  %6s  %8s  %5s  %10s  %10s\n",
    "rep", "time", "alpha_nm", "iters", "grad_norm", "hess_pd"))
for (i in 1:n_probe) {
  cat(sprintf("  %4d  %6.2f  %8.3f  %5d  %10.2e  %10s\n",
      i, times_total[i], alpha_norms[i], iter_counts[i],
      grad_norms[i], as.character(hess_pd_vals[i])))
}

# ─── 2. Isolate component costs ──────────────────────────────────────────────
cat("\n--- 2. Isolated component timing (rep 1) ---\n")
set.seed(12346)
dat <- all_data[1:n_val, ]
ebmr <- EBMRAlgorithmFast4$new("y", ps_spec_1, dat, W_func)

# Time just the hessian check portion
# We re-create the g/W functions manually to measure
alpha_dim <- 4L  # r ~ y + u1 + u2 → intercept + 3 = 4
h_dim <- 6L      # full_sq1: u1,u2,z1,z2,u2_sq + intercept d = 6

# Estimate once to get a representative alpha
res <- ebmr$EBMR_IPW(
  h_nu = function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1_u2=dat$u1*dat$u2),
  type = "HT", se.fit = FALSE
)
alpha_hat <- ebmr$ps_fit.list[[1]]$coefficients
cat(sprintf("  Representative alpha: [%s], norm=%.3f, hess_pd=%s\n",
    paste(round(alpha_hat, 3), collapse=", "),
    sqrt(sum(alpha_hat^2)),
    as.character(ebmr$ps_fit.list[[1]]$gmm_fit$opt$hessian_pd)))

# Time hessian() call alone
t_hess <- system.time({
  for (k in 1:10) {
    ebmr2 <- EBMRAlgorithmFast4$new("y", ps_spec_1, dat, W_func)
    ebmr2$EBMR_IPW(
      h_nu = function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1_u2=dat$u1*dat$u2),
      type = "HT", se.fit = FALSE
    )
  }
})["elapsed"]
cat(sprintf("  10× full EBMR_IPW (no SE): %.2f sec  → %.3f sec/call\n",
    t_hess, t_hess/10))

t_se <- system.time({
  for (k in 1:10) {
    ebmr2 <- EBMRAlgorithmFast4$new("y", ps_spec_1, dat, W_func)
    ebmr2$EBMR_IPW(
      h_nu = function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1_u2=dat$u1*dat$u2),
      type = "HT", se.fit = TRUE
    )
  }
})["elapsed"]
cat(sprintf("  10× full EBMR_IPW (with SE): %.2f sec  → %.3f sec/call\n",
    t_se, t_se/10))
cat(sprintf("  SE computation overhead: %.3f sec/call\n", (t_se - t_hess)/10))

# ─── 3. Estimate total simulation time ───────────────────────────────────────
cat("\n--- 3. Projected total simulation time ---\n")
n_reps <- 1000
n_cores <- max(1, parallel::detectCores() - 2)
time_per_rep <- mean(times_total)
projected_wall <- n_reps * time_per_rep / n_cores
cat(sprintf("  Average time/rep: %.2f sec\n", time_per_rep))
cat(sprintf("  Parallel cores: %d\n", n_cores))
cat(sprintf("  Projected wall time: %.1f min\n", projected_wall/60))

# ─── 4. Profile hessian check contribution ───────────────────────────────────
cat("\n--- 4. Profile: does hessian check dominate? ---\n")
cat("  The hessian check calls numDeriv::hessian(Q_full, alpha)\n")
cat("  For alpha_dim=4, this requires ~4*(4+1)/2 + 4 ≈ 12-20 Q_full evaluations\n")
cat("  Each Q_full computes g() + W = O(n * h_dim * alpha_dim)\n")
cat("  To measure directly, run with Rprof:\n")
cat("    Rprof('profile.out')\n")
cat("    <run EBMR_IPW 50 times>\n")
cat("    Rprof(NULL)\n")
cat("    summaryRprof('profile.out')\n\n")

# Approximate hessian cost
t_hess_only <- system.time({
  library(numDeriv)
  # Simulate what Q_full does
  ebmr3 <- EBMRAlgorithmFast4$new("y", ps_spec_1, dat, W_func)
  # Get a g function proxy
  invisible(ebmr3$EBMR_IPW(
    h_nu = function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1_u2=dat$u1*dat$u2),
    type="HT", se.fit=FALSE
  ))
})["elapsed"]

cat(sprintf("  Note: one rep includes hessian check.\n"))
cat(sprintf("  Run Rprof to get exact breakdown.\n"))

# ─── 5. Rprof profiling ───────────────────────────────────────────────────────
cat("\n--- 5. Rprof profiling (30 reps) ---\n")
Rprof("profile_93_1.out", interval = 0.01)
for (i in 1:30) {
  dat <- all_data[((i-1)*n_val + 1):(i*n_val), ]
  ebmr2 <- EBMRAlgorithmFast4$new("y", ps_spec_1, dat, W_func)
  ebmr2$EBMR_IPW(
    h_nu = function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1_u2=dat$u1*dat$u2),
    type = "HT", se.fit = TRUE
  )
}
Rprof(NULL)
prof_summary <- summaryRprof("profile_93_1.out")
cat("\nTop 15 functions by self time:\n")
print(head(prof_summary$by.self, 15))
cat("\nTop 15 functions by total time:\n")
print(head(prof_summary$by.total, 15))
