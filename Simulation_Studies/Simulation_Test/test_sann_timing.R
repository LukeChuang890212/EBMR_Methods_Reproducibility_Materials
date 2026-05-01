setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
library(EBMRalgorithmFast4)

data_file <- misspecified_model_all_data_file.list$setting3$miss50[[1]]
all_data <- readRDS(data_file)
n_val <- 2000

ps_spec <- get_ps_spec("9-alt1")
subset_ps_spec <- list(
  formula.list = ps_spec$formula.list[2],
  h_alpha.list = ps_spec$h_alpha.list[2],
  inv_link = ps_spec$inv_link,
  outcome = ps_spec$outcome
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

test_reps <- c(929, 438, 1, 2, 3)

cat("=== Timing comparison: current L-BFGS-B vs SANN ===\n\n")

for (rep_i in test_reps) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]

  # Time L-BFGS-B (current)
  t1 <- system.time({
    ebmr1 <- EBMRAlgorithmFast4$new("y", subset_ps_spec, dat, W_func)
  })

  fit1 <- ebmr1$ps_fit.list[[1]]
  mu1 <- mean(as.numeric(dat$r) / fit1$fitted.values * dat$y)

  cat(sprintf("Rep %d  L-BFGS-B: %.3fs  obj=%.2e  mu=%.4f  alpha=[%s]\n",
              rep_i, t1["elapsed"], fit1$gmm_fit$opt$objective, mu1,
              paste(round(fit1$coefficients, 4), collapse=", ")))

  # Now manually run SANN on the same GMM objective
  # Need to reconstruct the objective function
  parsed <- EBMRalgorithmFast4:::parse_formula(subset_ps_spec$formula.list[[1]],
                                                subset_ps_spec$outcome)
  dm <- EBMRalgorithmFast4:::build_design_matrix(dat, parsed, include_intercept = TRUE)
  design_matrix <- dm$design
  r_vec <- as.numeric(dat$r)
  inv_link <- subset_ps_spec$inv_link

  # Detect link type
  eta_max <- 20
  compute_pi <- function(eta) {
    eta <- pmin(pmax(eta, -eta_max), eta_max)
    1 / (1 + exp(eta))
  }

  # Use the h_x from the already-fitted model
  h_x <- fit1$h_x
  n <- nrow(dat)

  g_func <- function(param) {
    eta <- design_matrix %*% param
    pi_vec <- compute_pi(eta)
    rw <- r_vec / as.vector(pi_vec)
    (rw - 1) * h_x
  }

  # Simple GMM objective with identity W first, then update
  obj_sann <- function(param) {
    g_mat <- g_func(param)
    G_bar <- colMeans(g_mat)
    W_mat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(ncol(h_x)))
    as.numeric(t(G_bar) %*% W_mat %*% G_bar)
  }

  # Use MAR init as starting point for SANN too
  mar_formula <- as.formula(paste(parsed$r_name, "~", paste(parsed$x_terms, collapse = " + ")))
  mar_fit <- glm(mar_formula, data = dat, family = binomial())
  init_sann <- rep(0, ncol(design_matrix))
  init_sann[1] <- coef(mar_fit)[1]
  init_sann[(1 + dm$n_y + 1):ncol(design_matrix)] <- coef(mar_fit)[-1]

  t2 <- system.time({
    opt_sann <- optim(init_sann, obj_sann, method = "SANN",
                      control = list(maxit = 5000, temp = 10, tmax = 500))
  })

  pi_sann <- as.vector(compute_pi(design_matrix %*% opt_sann$par))
  mu_sann <- mean(r_vec / pi_sann * dat$y)

  cat(sprintf("Rep %d  SANN:     %.3fs  obj=%.2e  mu=%.4f  alpha=[%s]\n\n",
              rep_i, t2["elapsed"], opt_sann$value, mu_sann,
              paste(round(opt_sann$par, 4), collapse=", ")))
}

cat("Done!\n")
