setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
  devtools::load_all("../EPS")
})
W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
n_val <- 2000

ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)

# Test reps that hit 200 iterations at default max_outer
test_reps <- c(3, 5, 7, 9, 33)

# First run at default (200) to confirm cycling
cat("=== Default max_outer=200 ===\n")
for (i in test_reps) {
  dat <- all_data[((i-1)*n_val + 1):(i*n_val), ]
  ps_sub3 <- list(
    formula.list = ps_spec[["formula.list"]][3],
    h_alpha.list = ps_spec[["h_alpha.list"]][3],
    inv_link     = ps_spec[["inv_link"]],
    outcome      = ps_spec[["outcome"]]
  )
  ebmr <- EPS[["new"]]("y", ps_sub3, dat, W_func)
  gf <- ebmr[["ps_fit.list"]][[1]][["gmm_fit"]]
  cat(sprintf("  Rep %d: conv=%s, iter=%d, grad=%.2e, obj=%.6f, alpha=(%s)\n",
      i, gf[["opt"]][["converged"]], gf[["opt"]][["iterations"]],
      gf[["opt"]][["final_grad_norm"]], gf[["opt"]][["objective"]],
      paste(round(gf[["estimates"]], 4), collapse=", ")))
}

# Now test with higher max_outer by passing max_outer_override
# Need to rebuild WangShaoKim2014 internals manually
cat("\n=== Testing with max_outer = 500, 1000, 2000 ===\n")

for (mo in c(500, 1000, 2000)) {
  cat(sprintf("\n--- max_outer = %d ---\n", mo))
  for (i in test_reps) {
    dat <- all_data[((i-1)*n_val + 1):(i*n_val), ]

    # Create EBMR object but we'll re-run the GMM for model 3 with higher max_outer
    ps_sub3 <- list(
      formula.list = ps_spec[["formula.list"]][3],
      h_alpha.list = ps_spec[["h_alpha.list"]][3],
      inv_link     = ps_spec[["inv_link"]],
      outcome      = ps_spec[["outcome"]]
    )

    # Build the PS model manually to access gmm with max_outer_override
    # Step 1: fit GLM for init
    formula_j <- ps_sub3$formula.list[[1]]
    glm_fit <- glm(formula_j, data=dat, family=binomial(link="logit"))
    init <- coef(glm_fit)

    # Step 2: build moment function
    r_vec <- dat[["r"]]
    y_vec <- dat[["y"]]
    n <- n_val
    design_mat <- model.matrix(formula_j, data=dat)
    h_alpha_vars <- ps_sub3$h_alpha.list[[1]]
    h_alpha_mat <- as.matrix(dat[, h_alpha_vars, drop=FALSE])

    model_fn <- function(alpha) {
      eta <- as.vector(design_mat %*% alpha)
      1 / (1 + exp(-eta))
    }

    Phi_alpha <- function(alpha) {
      pi_hat <- model_fn(alpha)
      (r_vec / pi_hat - 1) * cbind(1, h_alpha_mat)
    }

    h_dim <- 1 + ncol(h_alpha_mat)
    alpha_dim <- length(init)

    # Analytical gradient
    dg_alpha <- function(alpha) {
      pi_hat <- model_fn(alpha)
      neg_r_pi2 <- -r_vec / (pi_hat^2)
      dpi <- design_mat * (pi_hat * (1 - pi_hat))
      base_mat <- neg_r_pi2 * dpi
      h_full <- cbind(1, h_alpha_mat)
      Gamma_arr <- array(0, dim=c(n, alpha_dim, h_dim))
      for (l in 1:h_dim) Gamma_arr[,,l] <- base_mat * h_full[,l]
      Gamma_arr
    }

    Gamma_direct <- function(alpha) {
      pi_hat <- model_fn(alpha)
      neg_r_pi2 <- -r_vec / (pi_hat^2)
      dpi <- design_mat * (pi_hat * (1 - pi_hat))
      h_full <- cbind(1, h_alpha_mat)
      crossprod(h_full * neg_r_pi2, dpi) / n
    }

    # Use the package's gmm via an EBMR object's private method
    # Simpler: just create the object and access the gmm
    ebmr <- EPS[["new"]]("y", ps_sub3, dat, W_func)

    # Re-run gmm with higher max_outer
    # Access private$gmm through environment
    gmm_env <- environment(ebmr[["EBMR_IPW"]])
    gmm_fn <- gmm_env$private$gmm
    W_fn <- gmm_env$private$W

    gf <- gmm_fn(Phi_alpha, W_fn, n, h_dim, alpha_dim, init, se.fit=FALSE,
                  dg=dg_alpha, Gamma_direct=Gamma_direct, max_outer_override=mo)

    cat(sprintf("  Rep %d: conv=%s, iter=%d, grad=%.2e, obj=%.6f, alpha=(%s)\n",
        i, gf[["opt"]][["converged"]], gf[["opt"]][["iterations"]],
        gf[["opt"]][["final_grad_norm"]], gf[["opt"]][["objective"]],
        paste(round(gf[["estimates"]], 4), collapse=", ")))
  }
}
