setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})

n_val <- 2000
n_reps <- 200
ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
mu_true <- get_mu_true("setting4")

formula_j <- ps_spec[["formula.list"]][[3]]
h_alpha_vars <- ps_spec[["h_alpha.list"]][[3]]

mu_vals <- rep(NA_real_, n_reps)
alpha_mat <- matrix(NA, n_reps, 4)
conv_vals <- rep(NA, n_reps)
obj_vals <- rep(NA_real_, n_reps)

for (rep_i in 1:n_reps) {
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
  tryCatch({
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

    # Fixed W = I: minimize G(alpha)'G(alpha)
    obj_fn <- function(alpha) {
      g_mat <- Phi_alpha(alpha)
      G <- colMeans(g_mat)
      sum(G^2)
    }

    opt <- optim(rep(0, p), obj_fn, method = "L-BFGS-B", control = list(maxit = 5000))
    alpha_final <- opt$par

    pi_hat <- model_fn(alpha_final)
    mu_vals[rep_i] <- mean(r_vec * dat[["y"]] / pi_hat)
    alpha_mat[rep_i, ] <- alpha_final
    conv_vals[rep_i] <- opt$convergence == 0
    obj_vals[rep_i] <- opt$value
  }, error = function(e) NULL)
}

valid <- !is.na(mu_vals)
conv <- conv_vals[valid]

cat(sprintf("mu.true = %.4f\n", mu_true))
cat(sprintf("Total valid: %d, converged: %d\n\n", sum(valid), sum(conv)))

esd <- sd(mu_vals[valid])
cat(sprintf("=== Overall ===\n"))
cat(sprintf("Bias=%.4f, ESD=%.4f\n\n", mean(mu_vals[valid]) - mu_true, esd))

colnames(alpha_mat) <- c("intercept", "y", "u2", "z2")
cat(sprintf("=== alpha ===\n"))
for (j in 1:4) {
  cat(sprintf("alpha[%d] (%s): mean=%.4f, sd=%.4f, range=[%.4f, %.4f]\n",
      j, colnames(alpha_mat)[j],
      mean(alpha_mat[valid, j]), sd(alpha_mat[valid, j]),
      min(alpha_mat[valid, j]), max(alpha_mat[valid, j])))
}

cat(sprintf("\n=== Objective G'G ===\n"))
cat(sprintf("mean=%.6f, sd=%.6f, range=[%.6f, %.6f]\n",
    mean(obj_vals[valid]), sd(obj_vals[valid]),
    min(obj_vals[valid]), max(obj_vals[valid])))

# mean pi
pi_means <- rep(NA, n_reps)
for (i in which(valid)) {
  dat <- all_data[((i-1)*n_val + 1):(i*n_val), ]
  dm <- model.matrix(formula_j, data=dat)
  eta <- as.vector(dm %*% alpha_mat[i,])
  pi_means[i] <- mean(1/(1+exp(-eta)))
}
cat(sprintf("\n=== mean(pi) ===\n"))
cat(sprintf("mean=%.4f, sd=%.4f, range=[%.4f, %.4f]\n",
    mean(pi_means[valid]), sd(pi_means[valid]),
    min(pi_means[valid]), max(pi_means[valid])))
