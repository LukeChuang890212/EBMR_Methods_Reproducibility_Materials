setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
  library(EBMRalgorithmFast4)
})

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_nu_fn <- function(dat) cbind(u1=dat[["u1"]], u2=dat[["u2"]], z1=dat[["z1"]], z2=dat[["z2"]], u1_u2=dat[["u1"]]*dat[["u2"]])
n_val <- 2000

ps_spec <- get_ps_spec("9-alt1")
ps_sub <- list(
  formula.list = ps_spec[["formula.list"]][2:3],
  h_alpha.list = ps_spec[["h_alpha.list"]][2:3],
  inv_link     = ps_spec[["inv_link"]],
  outcome      = ps_spec[["outcome"]]
)

data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)

# For a few reps, manually run the ensemble GMM with many different starting values
# and check if they converge to different solutions
n_check <- 20

cat("=== Checking nu sensitivity to starting values ===\n\n")

for (rep_i in 1:n_check) {
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]

  ebmr <- EBMRAlgorithmFast4[["new"]]("y", ps_sub, dat, W_func)

  # Get ps.matrix from the fitted models
  ps_mat <- do.call(cbind, lapply(1:2, function(j) ebmr[["ps_fit.list"]][[j]][["fitted.values"]]))
  r_vec <- dat[["r"]]
  n <- n_val

  # Build ensemble moment conditions manually
  h_x2 <- h_nu_fn(dat)
  h_x <- cbind(1, h_x2)  # add intercept like ensemble() does
  h_dim <- ncol(h_x)
  J <- 2

  Phi_nu <- function(param) {
    ps_nu <- as.vector(ps_mat %*% param)
    (r_vec / ps_nu - 1) * h_x
  }

  # Try many starting values
  inits <- list(
    c(0.5, 0.5),
    c(1, 0),
    c(0, 1),
    c(0.1, 0.9),
    c(0.9, 0.1),
    c(0.3, 0.7),
    c(0.7, 0.3),
    c(1, 1),
    c(2, 1),
    c(1, 2),
    c(0.01, 0.99),
    c(0.99, 0.01),
    c(-0.5, 1.5),
    c(1.5, -0.5),
    c(0.5, -0.5),
    c(-0.5, 0.5)
  )

  results <- matrix(NA, length(inits), 5)
  colnames(results) <- c("nu1", "nu2", "w1", "obj", "grad_norm")

  for (k in seq_along(inits)) {
    tryCatch({
      # Use the package's gmm function via a workaround:
      # Create a temporary EBMR object and call ensemble directly
      # Instead, manually run the GMM
      fit <- tryCatch({
        # Access private gmm through the ensemble mechanism
        # Simpler: just compute the GMM objective at convergence
        g_mat <- Phi_nu(inits[[k]])
        G_vec <- colMeans(g_mat)
        W_hat <- tryCatch(solve(crossprod(g_mat)/n), error=function(e) diag(h_dim))
        obj_init <- as.numeric(t(G_vec) %*% W_hat %*% G_vec)

        # Run optim with iterative GMM
        state_W <- diag(h_dim)
        obj_fn <- function(param) {
          g_m <- Phi_nu(param)
          G_v <- matrix(colMeans(g_m), h_dim, 1)
          as.numeric(crossprod(G_v, state_W %*% G_v))
        }

        # Step 1: W=I
        opt1 <- optim(inits[[k]], obj_fn, method="L-BFGS-B", control=list(maxit=1000))
        est <- opt1$par

        # Step 2: iterate
        for (tt in 1:200) {
          g_m <- Phi_nu(est)
          state_W <- tryCatch(solve(crossprod(g_m)/n), error=function(e) diag(h_dim))
          opt_t <- optim(est, obj_fn, method="L-BFGS-B", control=list(maxit=1000))
          if (max(abs(opt_t$par - est)) < 1e-6) break
          est <- opt_t$par
        }

        g_final <- Phi_nu(est)
        G_final <- matrix(colMeans(g_final), h_dim, 1)
        W_final <- tryCatch(solve(crossprod(g_final)/n), error=function(e) diag(h_dim))
        obj_final <- as.numeric(crossprod(G_final, W_final %*% G_final))

        # Compute gradient norm
        dg_fn <- function(param) {
          ps_nu <- as.vector(ps_mat %*% param)
          neg_r_ps2 <- -r_vec / (ps_nu^2)
          base_mat <- neg_r_ps2 * ps_mat
          Gamma_arr <- array(0, dim = c(n, J, h_dim))
          for (l in 1:h_dim) Gamma_arr[,,l] <- base_mat * h_x[,l]
          Gamma_arr
        }
        dg_arr <- dg_fn(est)
        Gamma_mat <- matrix(0, h_dim, J)
        for (j in 1:J) Gamma_mat[,j] <- colMeans(dg_arr[,j,])
        grad_vec <- 2 * as.vector(crossprod(Gamma_mat, W_final %*% G_final))

        w_est <- est^2 / sum(est^2)
        list(nu=est, w=w_est, obj=obj_final, grad=max(abs(grad_vec)))
      }, error = function(e) NULL)

      if (!is.null(fit)) {
        results[k,] <- c(fit$nu, fit$w[1], fit$obj, fit$grad)
      }
    }, error = function(e) NULL)
  }

  # Report
  cat(sprintf("Rep %d:\n", rep_i))

  # Package result for comparison
  pkg_res <- ebmr[["EBMR_IPW"]](h_nu=h_nu_fn, type="HT", se.fit=FALSE)
  pkg_w <- pkg_res[["w.hat"]]
  pkg_nu <- pkg_res[["nu.hat"]]
  cat(sprintf("  Package: nu=(%.4f, %.4f)  w=(%.4f, %.4f)\n",
      pkg_nu[1], pkg_nu[2], pkg_w[1], pkg_w[2]))

  # Group by converged solution (cluster by w1)
  valid <- !is.na(results[,1])
  if (sum(valid) > 0) {
    w1_vals <- results[valid, "w1"]
    obj_vals <- results[valid, "obj"]
    nu1_vals <- results[valid, "nu1"]
    nu2_vals <- results[valid, "nu2"]

    # Find distinct solutions (w1 differs by > 0.05)
    ord <- order(w1_vals)
    w1_sorted <- w1_vals[ord]
    groups <- cumsum(c(TRUE, diff(w1_sorted) > 0.05))
    n_solutions <- max(groups)

    cat(sprintf("  %d distinct solutions found from %d starts:\n", n_solutions, sum(valid)))
    for (g in 1:n_solutions) {
      idx <- ord[groups == g]
      cat(sprintf("    Solution %d: w1=%.4f (n=%d), obj=%.6f, nu=(%.4f,%.4f)\n",
          g, mean(w1_vals[idx]), length(idx), mean(obj_vals[idx]),
          mean(nu1_vals[idx]), mean(nu2_vals[idx])))
    }
  }
  cat("\n")
}
