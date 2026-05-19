## Diagnose the ~8% that flip to M3 in two-step GMM
## Questions:
## 1. At which step do they flip? Step 1 or Step 2?
## 2. What's different about these datasets?
## 3. Is the M3 solution a lower CUE objective?
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

ps_spec_base <- get_ps_spec("9")
h_alpha_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2)
W_sm <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_nu_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1u2=dat$u1*dat$u2)

data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]]
all_data <- NULL
for (fi in seq_along(data_file)) {
  if (file.exists(data_file[[fi]])) {
    test <- readRDS(data_file[[fi]])
    if (nrow(test) / 1000 >= 2) { all_data <- test; break }
  }
}

nn <- 2000
ps_spec_m23 <- list(
  formula.list = ps_spec_base$formula.list[2:3],
  h_alpha.list = list(h_alpha_fn, h_alpha_fn),
  inv_link = ps_spec_base$inv_link,
  outcome = ps_spec_base$outcome
)

cat("=== Diagnose: which reps flip to M3, at which step, and why ===\n\n")

n_check <- 100
flip_at_s1 <- 0; flip_at_s2 <- 0; stay_m2 <- 0

cat(sprintf("%-4s | %-8s | %-8s | %-8s | %-12s | %-12s | %-12s | %s\n",
    "Rep", "w1_init", "w1_s1", "w1_s2", "obj_s1(I)", "obj_s2(W)", "CUE_obj", "Status"))
cat(strrep("-", 110), "\n")

for (i in 1:n_check) {
  dat <- all_data[((i-1)*nn + 1):(i*nn), ]
  ebmr <- EBMRAlgorithmFast4$new("y", ps_spec_m23, dat, W_sm)
  ps_mat <- do.call(cbind, lapply(ebmr$ps_fit.list, function(pf) pf$fitted.values))
  r_vec <- as.vector(dat$r)
  h_x <- cbind(1, h_nu_fn(dat))
  h_dim <- ncol(h_x); n <- nn

  g_nu <- function(nu) {
    ps_nu <- as.vector(ps_mat %*% nu)
    as.vector(r_vec / ps_nu - 1) * h_x
  }
  G_nu <- function(nu) colMeans(g_nu(nu))

  # CUE objective
  cue_obj <- function(nu) {
    g_mat <- g_nu(nu)
    G <- colMeans(g_mat)
    W <- tryCatch(solve(t(g_mat) %*% g_mat / n), error = function(e) diag(h_dim))
    as.numeric(t(G) %*% W %*% G)
  }

  nu_init <- c(0.95, 0.05)
  w_init <- nu_init^2 / sum(nu_init^2)

  # Step 1: W = I
  obj_I <- function(nu) { G <- G_nu(nu); sum(G^2) }
  opt1 <- optim(nu_init, obj_I, method = "L-BFGS-B", control = list(maxit = 1000))
  nu_s1 <- opt1$par
  w_s1 <- nu_s1^2 / sum(nu_s1^2)
  obj_s1 <- opt1$value

  # Step 2: one W update
  g1 <- g_nu(nu_s1)
  W2 <- tryCatch(solve(t(g1) %*% g1 / n), error = function(e) diag(h_dim))
  obj2_fn <- function(nu) { G <- matrix(G_nu(nu), ncol=1); as.numeric(t(G) %*% W2 %*% G) }
  opt2 <- optim(nu_s1, obj2_fn, method = "L-BFGS-B", control = list(maxit = 1000))
  nu_s2 <- opt2$par
  w_s2 <- nu_s2^2 / sum(nu_s2^2)
  obj_s2 <- opt2$value

  cue_val <- cue_obj(nu_s2)

  # Classify
  if (w_s1[1] < 0.5) {
    status <- "FLIP@S1"
    flip_at_s1 <- flip_at_s1 + 1
  } else if (w_s2[1] < 0.5) {
    status <- "FLIP@S2"
    flip_at_s2 <- flip_at_s2 + 1
  } else {
    status <- "M2"
    stay_m2 <- stay_m2 + 1
  }

  # Only print flips and a few M2 for comparison
  if (status != "M2" || i <= 5) {
    cat(sprintf("%4d | %.4f  | %.4f  | %.4f  | %.4e    | %.4e    | %.4e    | %s\n",
        i, w_init[1], w_s1[1], w_s2[1], obj_s1, obj_s2, cue_val, status))
  }
}

cat(sprintf("\n--- Summary (first %d reps) ---\n", n_check))
cat(sprintf("Stay M2:  %d (%.1f%%)\n", stay_m2, 100*stay_m2/n_check))
cat(sprintf("Flip@S1:  %d (%.1f%%)\n", flip_at_s1, 100*flip_at_s1/n_check))
cat(sprintf("Flip@S2:  %d (%.1f%%)\n", flip_at_s2, 100*flip_at_s2/n_check))

# For flipped reps, compare M2 vs M3 CUE objectives
cat("\n\n=== For flipped reps: compare CUE objective at M2 vs M3 solutions ===\n\n")
cat(sprintf("%-4s | %-12s | %-12s | %s\n", "Rep", "CUE@M2", "CUE@M3", "Lower"))
cat(strrep("-", 55), "\n")

for (i in 1:n_check) {
  dat <- all_data[((i-1)*nn + 1):(i*nn), ]
  ebmr <- EBMRAlgorithmFast4$new("y", ps_spec_m23, dat, W_sm)
  ps_mat <- do.call(cbind, lapply(ebmr$ps_fit.list, function(pf) pf$fitted.values))
  r_vec <- as.vector(dat$r)
  h_x <- cbind(1, h_nu_fn(dat))
  h_dim <- ncol(h_x); n <- nn

  g_nu <- function(nu) {
    ps_nu <- as.vector(ps_mat %*% nu)
    as.vector(r_vec / ps_nu - 1) * h_x
  }
  G_nu <- function(nu) colMeans(g_nu(nu))

  cue_obj <- function(nu) {
    g_mat <- g_nu(nu)
    G <- colMeans(g_mat)
    W <- tryCatch(solve(t(g_mat) %*% g_mat / n), error = function(e) diag(h_dim))
    as.numeric(t(G) %*% W %*% G)
  }

  # Run from M2-favoring init
  obj_I_m2 <- function(nu) { G <- G_nu(nu); sum(G^2) }
  opt_m2_s1 <- optim(c(0.99, 0.01), obj_I_m2, method = "L-BFGS-B", control = list(maxit = 1000))
  g_m2 <- g_nu(opt_m2_s1$par)
  W_m2 <- tryCatch(solve(t(g_m2) %*% g_m2 / n), error = function(e) diag(h_dim))
  obj_m2_fn <- function(nu) { G <- matrix(G_nu(nu), ncol=1); as.numeric(t(G) %*% W_m2 %*% G) }
  opt_m2_s2 <- optim(opt_m2_s1$par, obj_m2_fn, method = "L-BFGS-B", control = list(maxit = 1000))
  w_m2 <- opt_m2_s2$par^2 / sum(opt_m2_s2$par^2)

  # Run from M3-favoring init
  opt_m3_s1 <- optim(c(0.01, 0.99), obj_I_m2, method = "L-BFGS-B", control = list(maxit = 1000))
  g_m3 <- g_nu(opt_m3_s1$par)
  W_m3 <- tryCatch(solve(t(g_m3) %*% g_m3 / n), error = function(e) diag(h_dim))
  obj_m3_fn <- function(nu) { G <- matrix(G_nu(nu), ncol=1); as.numeric(t(G) %*% W_m3 %*% G) }
  opt_m3_s2 <- optim(opt_m3_s1$par, obj_m3_fn, method = "L-BFGS-B", control = list(maxit = 1000))
  w_m3 <- opt_m3_s2$par^2 / sum(opt_m3_s2$par^2)

  cue_m2 <- cue_obj(opt_m2_s2$par)
  cue_m3 <- cue_obj(opt_m3_s2$par)

  # Check if this rep flips with (0.95, 0.05) init
  opt_check <- optim(c(0.95, 0.05), obj_I_m2, method = "L-BFGS-B", control = list(maxit = 1000))
  g_ch <- g_nu(opt_check$par)
  W_ch <- tryCatch(solve(t(g_ch) %*% g_ch / n), error = function(e) diag(h_dim))
  obj_ch_fn <- function(nu) { G <- matrix(G_nu(nu), ncol=1); as.numeric(t(G) %*% W_ch %*% G) }
  opt_ch2 <- optim(opt_check$par, obj_ch_fn, method = "L-BFGS-B", control = list(maxit = 1000))
  w_ch <- opt_ch2$par^2 / sum(opt_ch2$par^2)

  if (w_ch[1] < 0.5) {
    lower <- if (abs(cue_m2 - cue_m3) < 1e-8) "TIE" else if (cue_m2 < cue_m3) "M2" else "M3"
    cat(sprintf("%4d | %.4e    | %.4e    | %s  (w_M2=%.3f, w_M3=%.3f)\n",
        i, cue_m2, cue_m3, lower, w_m2[1], w_m3[1]))
  }
}
