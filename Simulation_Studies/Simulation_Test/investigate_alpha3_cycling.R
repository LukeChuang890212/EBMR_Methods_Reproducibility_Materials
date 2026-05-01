setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
  devtools::load_all("../EBMRalgorithmFast4")
})
W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
n_val <- 2000

ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)

formula_j <- ps_spec[["formula.list"]][[3]]
h_alpha_vars <- ps_spec[["h_alpha.list"]][[3]]

# Trace rep 3 (one of the cycling reps)
for (rep_i in c(3, 7)) {
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

  glm_fit <- glm(formula_j, data=dat, family=binomial(link="logit"))
  alpha_glm <- coef(glm_fit)
  cat(sprintf("GLM init: (%s)\n", paste(round(alpha_glm, 4), collapse=", ")))

  # Trace iterative GMM
  alpha_t <- alpha_glm
  sign_history <- c()

  for (t in 1:50) {
    g_t <- Phi_alpha(alpha_t)
    G_t <- colMeans(g_t)
    W_t <- tryCatch(solve(crossprod(g_t)/n), error=function(e) diag(ncol(h_full)))
    obj_t <- as.numeric(t(G_t) %*% W_t %*% G_t)

    # Convergence gradient
    pi_hat <- model_fn(alpha_t)
    neg_r_pi2 <- -r_vec / (pi_hat^2)
    dpi <- design_mat * (pi_hat * (1 - pi_hat))
    Gamma_mat <- crossprod(h_full * neg_r_pi2, dpi) / n
    conv_grad <- 2 * as.vector(crossprod(Gamma_mat, W_t %*% G_t))
    grad_norm <- max(abs(conv_grad))

    obj_fn <- function(a) {
      g_m <- Phi_alpha(a)
      G_v <- colMeans(g_m)
      as.numeric(t(G_v) %*% W_t %*% G_v)
    }
    opt <- optim(alpha_t, obj_fn, method="L-BFGS-B", control=list(maxit=1000))
    alpha_new <- opt$par

    if (t <= 20 || t %% 10 == 0) {
      cat(sprintf("  t=%2d: alpha[1]=% .4f, alpha[2]=% .4f, obj=%.6f, grad=%.2e, ||step||=%.4f\n",
          t, alpha_t[1], alpha_t[2], obj_t, grad_norm, sqrt(sum((alpha_new - alpha_t)^2))))
    }

    sign_history <- c(sign_history, sign(alpha_t[1]))
    alpha_t <- alpha_new
  }

  # Check if cycling between signs
  cat(sprintf("\nSign of alpha[1] over iterations: "))
  cat(paste(ifelse(sign_history > 0, "+", "-"), collapse=""))
  cat("\n")

  # Compare objectives at final alpha and its negative
  alpha_final <- alpha_t
  alpha_flipped <- -alpha_final
  g_f <- Phi_alpha(alpha_final); W_f <- solve(crossprod(g_f)/n)
  g_fl <- Phi_alpha(alpha_flipped); W_fl <- solve(crossprod(g_fl)/n)
  cat(sprintf("\nFinal alpha: (%s)\n", paste(round(alpha_final, 4), collapse=", ")))
  cat(sprintf("Obj(final, W(final)):     %.8f\n", as.numeric(t(colMeans(g_f)) %*% W_f %*% colMeans(g_f))))
  cat(sprintf("Obj(-final, W(-final)):   %.8f\n", as.numeric(t(colMeans(g_fl)) %*% W_fl %*% colMeans(g_fl))))

  # What's happening: check pi values
  pi_final <- model_fn(alpha_final)
  pi_flipped <- model_fn(alpha_flipped)
  cat(sprintf("\npi(final): mean=%.4f, range=[%.4f, %.4f]\n", mean(pi_final), min(pi_final), max(pi_final)))
  cat(sprintf("pi(-final): mean=%.4f, range=[%.4f, %.4f]\n", mean(pi_flipped), min(pi_flipped), max(pi_flipped)))
  cat(sprintf("Missing rate: %.4f\n", mean(r_vec)))
}
