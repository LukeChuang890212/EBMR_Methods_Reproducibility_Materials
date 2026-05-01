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

# Reps that hit the -10 bound (from previous check, without bounds)
# Rep 15: alpha[2]=-18, Rep 34: alpha[2]=-21, Rep 35: alpha[2]=-18
test_reps <- c(15, 34, 35, 68, 161)

for (rep_i in test_reps) {
  cat(sprintf("\n========== Rep %d ==========\n", rep_i))
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]

  r_vec <- dat[["r"]]
  design_mat <- model.matrix(formula_j, data=dat)
  h_alpha_mat <- as.matrix(dat[, h_alpha_vars, drop=FALSE])
  h_full <- cbind(1, h_alpha_mat)
  n <- n_val

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

  # Trace iterative GMM without bounds
  alpha_t <- alpha_glm
  max_iter <- 50

  cat(sprintf("GLM init: alpha[2]=%.4f\n", alpha_glm[2]))

  for (t in 1:max_iter) {
    g_t <- Phi_alpha(alpha_t)
    G_t <- colMeans(g_t)
    W_t <- tryCatch(solve(crossprod(g_t)/n), error=function(e) diag(ncol(h_full)))
    obj_t <- as.numeric(t(G_t) %*% W_t %*% G_t)

    # Inner optim (no bounds)
    obj_fn <- function(a) {
      g_m <- Phi_alpha(a)
      G_v <- colMeans(g_m)
      as.numeric(t(G_v) %*% W_t %*% G_v)
    }
    opt <- optim(alpha_t, obj_fn, method="L-BFGS-B", control=list(maxit=1000))
    alpha_new <- opt$par

    step_size <- sqrt(sum((alpha_new - alpha_t)^2))
    pi_mean <- mean(model_fn(alpha_t))

    cat(sprintf("  t=%2d: a2=%8.3f, obj=%.6f, ||step||=%7.4f, mean_pi=%.4f\n",
        t, alpha_t[2], obj_t, step_size, pi_mean))

    # Stop if converged or alpha[2] goes extreme
    if (step_size < 1e-6 && t > 3) break
    if (abs(alpha_new[2]) > 40) {
      cat(sprintf("  --> alpha[2] exploded to %.2f, stopping\n", alpha_new[2]))
      break
    }

    alpha_t <- alpha_new
  }
  cat(sprintf("  Final: alpha=(%s)\n", paste(round(alpha_t, 3), collapse=", ")))
}
