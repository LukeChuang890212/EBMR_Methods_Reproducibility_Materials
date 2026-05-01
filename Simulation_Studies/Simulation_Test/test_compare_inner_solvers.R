setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})

n_val <- 2000
ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)

formula_j <- ps_spec[["formula.list"]][[3]]
h_alpha_vars <- ps_spec[["h_alpha.list"]][[3]]

# Test a few reps: compare what L-BFGS-B vs NR find for the same fixed W
test_reps <- c(1, 3, 5, 7, 10, 15, 20, 34, 50)

for (rep_i in test_reps) {
  cat(sprintf("\n========== Rep %d ==========\n", rep_i))
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]

  r_vec <- dat[["r"]]
  design_mat <- model.matrix(formula_j, data=dat)
  h_alpha_mat <- as.matrix(dat[, h_alpha_vars, drop=FALSE])
  h_full <- cbind(1, h_alpha_mat)
  n <- n_val
  p <- ncol(design_mat)
  h_dim <- ncol(h_full)

  model_fn <- function(alpha) {
    eta <- as.vector(design_mat %*% alpha)
    1 / (1 + exp(-eta))
  }

  Phi_alpha <- function(alpha) {
    pi_hat <- model_fn(alpha)
    (r_vec / pi_hat - 1) * h_full
  }

  Gamma_fn <- function(alpha) {
    pi_hat <- model_fn(alpha)
    common <- -r_vec * (1 - pi_hat) / pi_hat
    crossprod(h_full * common, design_mat) / n
  }

  W_func <- function(g_mat) solve(crossprod(g_mat) / nrow(g_mat))

  # === Step 1: Get initial estimate with W=I using L-BFGS-B ===
  obj_I <- function(alpha) {
    g_mat <- Phi_alpha(alpha)
    G <- colMeans(g_mat)
    sum(G^2)
  }
  opt0 <- optim(rep(0, p), obj_I, method = "L-BFGS-B", control = list(maxit = 5000))
  alpha_step1 <- opt0$par

  # === Step 2: Compute W at step 1 estimate ===
  g_mat <- Phi_alpha(alpha_step1)
  W_hat <- tryCatch(W_func(g_mat), error = function(e) diag(h_dim))

  # Fixed-W objective
  obj_fixedW <- function(alpha) {
    g_mat <- Phi_alpha(alpha)
    G <- colMeans(g_mat)
    as.numeric(t(G) %*% W_hat %*% G)
  }

  # === Method A: L-BFGS-B from step1 estimate ===
  optA <- optim(alpha_step1, obj_fixedW, method = "L-BFGS-B", control = list(maxit = 5000))
  alphaA <- optA$par
  objA <- optA$value
  piA <- mean(model_fn(alphaA))

  # === Method B: L-BFGS-B from zero ===
  optB <- optim(rep(0, p), obj_fixedW, method = "L-BFGS-B", control = list(maxit = 5000))
  alphaB <- optB$par
  objB <- optB$value
  piB <- mean(model_fn(alphaB))

  # === Method C: L-BFGS-B with multiple random starts ===
  best_obj <- Inf
  best_alpha <- NULL
  for (s in 1:10) {
    init_s <- rnorm(p, sd = 0.5)
    opt_s <- optim(init_s, obj_fixedW, method = "L-BFGS-B", control = list(maxit = 5000))
    if (opt_s$value < best_obj) {
      best_obj <- opt_s$value
      best_alpha <- opt_s$par
    }
  }
  objC <- best_obj
  alphaC <- best_alpha
  piC <- mean(model_fn(alphaC))

  cat(sprintf("Step 1 alpha: [%.3f, %.3f, %.3f, %.3f], mean(pi)=%.4f\n",
      alpha_step1[1], alpha_step1[2], alpha_step1[3], alpha_step1[4],
      mean(model_fn(alpha_step1))))

  cat(sprintf("  A (LBFGSB from step1): obj=%.8f, alpha[2]=%.4f, mean(pi)=%.4f\n", objA, alphaA[2], piA))
  cat(sprintf("  B (LBFGSB from zero):  obj=%.8f, alpha[2]=%.4f, mean(pi)=%.4f\n", objB, alphaB[2], piB))
  cat(sprintf("  C (best of 10 random): obj=%.8f, alpha[2]=%.4f, mean(pi)=%.4f\n", objC, alphaC[2], piC))

  # Check: is the spurious solution actually lower Q?
  # Try extreme alpha[2]
  alpha_extreme <- alpha_step1
  alpha_extreme[2] <- 15
  # Re-optimize other params with alpha[2] fixed at 15
  obj_partial <- function(a_rest) {
    alpha_full <- c(a_rest[1], 15, a_rest[2], a_rest[3])
    obj_fixedW(alpha_full)
  }
  opt_ext <- optim(c(alpha_step1[1], alpha_step1[3], alpha_step1[4]), obj_partial,
                   method = "L-BFGS-B", control = list(maxit = 5000))
  alpha_ext <- c(opt_ext$par[1], 15, opt_ext$par[2], opt_ext$par[3])
  obj_ext <- obj_fixedW(alpha_ext)
  pi_ext <- mean(model_fn(alpha_ext))

  cat(sprintf("  D (alpha[2]=15 fixed): obj=%.8f, mean(pi)=%.4f\n", obj_ext, pi_ext))
}
