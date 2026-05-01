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

# Pick reps that converged to extreme alpha in previous run
test_reps <- c(3, 7, 15, 34, 35)

for (rep_i in test_reps) {
  cat(sprintf("\n========== Rep %d ==========\n", rep_i))
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]

  r_vec <- dat[["r"]]
  design_mat <- model.matrix(formula_j, data=dat)
  h_alpha_mat <- as.matrix(dat[, h_alpha_vars, drop=FALSE])
  h_full <- cbind(1, h_alpha_mat)
  n <- n_val
  p <- ncol(design_mat)

  model_fn <- function(alpha) {
    eta <- as.vector(design_mat %*% alpha)
    1 / (1 + exp(-eta))
  }

  Phi_alpha <- function(alpha) {
    pi_hat <- model_fn(alpha)
    (r_vec / pi_hat - 1) * h_full
  }

  # GLM init
  glm_fit <- glm(formula_j, data=dat, family=binomial(link="logit"))
  alpha_glm <- coef(glm_fit)
  cat(sprintf("GLM init: (%s)\n", paste(round(alpha_glm, 4), collapse=", ")))

  # Step 1: W0 from zero initials (matching current package code)
  alpha_zero <- rep(0, p)
  g_zero <- Phi_alpha(alpha_zero)
  W0 <- tryCatch(solve(crossprod(g_zero)/n), error=function(e) diag(ncol(h_full)))

  # Optimize from GLM init with W0
  obj_fn_w0 <- function(a) {
    g_m <- Phi_alpha(a)
    G_v <- colMeans(g_m)
    as.numeric(t(G_v) %*% W0 %*% G_v)
  }
  opt1 <- optim(alpha_glm, obj_fn_w0, method="L-BFGS-B", control=list(maxit=1000))
  cat(sprintf("Step 1 (W0 at zero): alpha=(%s)\n", paste(round(opt1$par, 4), collapse=", ")))
  cat(sprintf("  obj=%.6f, mean_pi=%.4f\n", opt1$value, mean(model_fn(opt1$par))))

  # Then trace iterative GMM from step 1 result
  alpha_t <- opt1$par
  for (t in 1:10) {
    g_t <- Phi_alpha(alpha_t)
    G_t <- colMeans(g_t)
    W_t <- tryCatch(solve(crossprod(g_t)/n), error=function(e) diag(ncol(h_full)))
    obj_t <- as.numeric(t(G_t) %*% W_t %*% G_t)

    obj_fn <- function(a) {
      g_m <- Phi_alpha(a)
      G_v <- colMeans(g_m)
      as.numeric(t(G_v) %*% W_t %*% G_v)
    }
    opt <- optim(alpha_t, obj_fn, method="L-BFGS-B", control=list(maxit=1000))
    alpha_new <- opt$par
    step_size <- sqrt(sum((alpha_new - alpha_t)^2))

    cat(sprintf("  t=%2d: alpha[2]=%8.3f, obj=%.6f, ||step||=%7.4f, mean_pi=%.4f\n",
        t, alpha_t[2], obj_t, step_size, mean(model_fn(alpha_t))))

    if (step_size < 1e-6 && t > 2) break
    alpha_t <- alpha_new
  }
}
