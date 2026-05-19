setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EPS", quiet = TRUE)

n_val <- 2000
ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

model_idx <- 3
dat <- all_data[1:n_val, ]
single_ps <- list(
  formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
  h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]],
  alpha_init.list = list(NULL),
  optimizer = "constrained_nr"
)
ebmr <- EPS$new("y", single_ps, dat, W_fn)
gmm_fit <- ebmr$ps_fit.list[[1]]$gmm_fit
design_mat <- ebmr$ps_fit.list[[1]]$design_matrix
alpha_hat <- gmm_fit$estimates
pi_hat <- ebmr$ps_fit.list[[1]]$fitted.values
n <- n_val; k <- length(alpha_hat)
r_vec <- dat[["r"]]; y_vec <- dat[["y"]]

h_alpha_vars <- ps_spec[["h_alpha.list"]][[model_idx]]
h_x <- cbind(1, as.matrix(dat[, h_alpha_vars, drop = FALSE]))
h <- ncol(h_x)

compute_pi_fn <- function(eta) plogis(-eta)

g_fn <- function(alpha) {
  eta <- as.vector(design_mat %*% alpha)
  pi_v <- compute_pi_fn(eta)
  (r_vec / pi_v - 1) * h_x
}

# Verify: does g_fn(alpha_hat) match gmm_fit$g.matrix?
g_mine <- g_fn(alpha_hat)
g_pkg <- gmm_fit$g.matrix
cat(sprintf("g_fn vs package g.matrix: max diff = %.2e\n", max(abs(g_mine - g_pkg))))

# Step 1: Verify dg/dalpha_j numerically for j=1
j <- 1
eps <- 1e-6
dg_num <- (g_fn(alpha_hat + eps * (1:k == j)) - g_fn(alpha_hat - eps * (1:k == j))) / (2*eps)

# Analytical dg for j=1
cf <- r_vec * (1 - pi_hat) / pi_hat
dg_ana <- cf * design_mat[, j] * h_x  # n x h

cat(sprintf("\ndg/dalpha_%d: max diff = %.2e, rel = %.2e\n", j,
    max(abs(dg_num - dg_ana)), max(abs(dg_num - dg_ana)) / max(abs(dg_ana))))

# Step 2: Compute S_mat column j analytically and numerically
# Analytical
cross_ana <- (crossprod(dg_ana, g_pkg) + crossprod(g_pkg, dg_ana)) / n
S_col_ana <- as.vector(cross_ana)

# Numerical
g_plus <- g_fn(alpha_hat + eps * (1:k == j))
g_minus <- g_fn(alpha_hat - eps * (1:k == j))
Winv_plus <- crossprod(g_plus) / n
Winv_minus <- crossprod(g_minus) / n
S_col_num <- as.vector((Winv_plus - Winv_minus) / (2*eps))

cat(sprintf("\nS_mat col %d:\n", j))
cat(sprintf("  max abs diff = %.6e\n", max(abs(S_col_ana - S_col_num))))
cat(sprintf("  rel diff     = %.6e\n", max(abs(S_col_ana - S_col_num)) / max(abs(S_col_num))))

# Print detailed comparison
cat("\n  S_analytical vs S_numerical:\n")
for (idx in 1:length(S_col_ana)) {
  if (abs(S_col_ana[idx] - S_col_num[idx]) > 0.001 * max(abs(S_col_num))) {
    l1 <- ((idx-1) %% h) + 1
    l2 <- ((idx-1) %/% h) + 1
    cat(sprintf("  [%d,%d]: ana=%.6f, num=%.6f, diff=%.6f\n", l1, l2,
        S_col_ana[idx], S_col_num[idx], S_col_ana[idx] - S_col_num[idx]))
  }
}

# Step 3: Decompose the issue
# The numerical S captures: d(g'g/n)/dalpha_j
# = (1/n) [dg'*g + g'*dg] where dg includes BOTH the explicit alpha derivative AND
# the effect of alpha on g through the estimating equations.
# Wait - dg_ana already captures the full derivative of g_i(alpha) w.r.t. alpha_j.
# So the analytical and numerical should match.

# Let me check: is there an issue with g_fn not matching g.matrix?
cat(sprintf("\nSanity check: g_fn(alpha_hat) max diff from g.matrix = %.2e\n",
    max(abs(g_fn(alpha_hat) - g_pkg))))

# Step 4: Maybe the issue is with the package's t_Gamma_arr?
# Let me compare the package's S_mat with my analytical S_mat
# To access the package's S_mat, I need to look at what gmm_fit stores
# The package computes S_mat internally and uses it in H2_mat
# Let me reconstruct it from the package's data

# Actually, let me just check all 4 columns
cat("\n=== Full S_mat comparison ===\n")
S_ana_full <- matrix(0, h^2, k)
S_num_full <- matrix(0, h^2, k)
for (jj in 1:k) {
  dg_j_ana <- cf * design_mat[, jj] * h_x
  cross_j <- (crossprod(dg_j_ana, g_pkg) + crossprod(g_pkg, dg_j_ana)) / n
  S_ana_full[, jj] <- as.vector(cross_j)

  g_p <- g_fn(alpha_hat + eps * (1:k == jj))
  g_m <- g_fn(alpha_hat - eps * (1:k == jj))
  S_num_full[, jj] <- as.vector((crossprod(g_p) - crossprod(g_m)) / (2*n*eps))
}
cat(sprintf("||S_ana - S_num||: %.6e\n", norm(S_ana_full - S_num_full, "F")))
cat(sprintf("||S_num||: %.6f\n", norm(S_num_full, "F")))
cat(sprintf("Relative: %.6e\n", norm(S_ana_full - S_num_full, "F") / norm(S_num_full, "F")))

# Step 5: The analytical S captures (1/n)(dg'g + g'dg)
# The numerical S captures d((1/n)g'g)/dalpha
# By product rule: d(g'g)/dalpha_j = dg'g + g'dg
# So they should match exactly. Unless dg_ana != dg_num!

cat("\n=== dg verification for all columns ===\n")
for (jj in 1:k) {
  dg_j_num <- (g_fn(alpha_hat + eps * (1:k == jj)) - g_fn(alpha_hat - eps * (1:k == jj))) / (2*eps)
  dg_j_ana <- cf * design_mat[, jj] * h_x
  rel <- norm(dg_j_num - dg_j_ana, "F") / norm(dg_j_num, "F")
  cat(sprintf("  j=%d: ||dg_num - dg_ana|| / ||dg_num|| = %.6e\n", jj, rel))
}

# Step 6: Maybe g used in S differs? Use dg_num in S computation
cat("\n=== S_mat using numerical dg ===\n")
S_dg_num_full <- matrix(0, h^2, k)
for (jj in 1:k) {
  dg_j_num <- (g_fn(alpha_hat + eps * (1:k == jj)) - g_fn(alpha_hat - eps * (1:k == jj))) / (2*eps)
  cross_j <- (crossprod(dg_j_num, g_pkg) + crossprod(g_pkg, dg_j_num)) / n
  S_dg_num_full[, jj] <- as.vector(cross_j)
}
cat(sprintf("||S_dg_num - S_num||: %.6e (rel=%.6e)\n",
    norm(S_dg_num_full - S_num_full, "F"),
    norm(S_dg_num_full - S_num_full, "F") / norm(S_num_full, "F")))
