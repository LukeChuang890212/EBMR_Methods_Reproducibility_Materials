setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Data_Generation.r")
source("config/scenarios.R")
library(EPS)
library(numDeriv)
source("Basic_setup.r")

data_file <- misspecified_model_all_data_file.list$setting4$miss50[[1]]
all_data <- readRDS(data_file)
n_val <- 2000

ps_spec <- get_ps_spec("9-alt1")
subset_ps_spec <- list(
  formula.list = ps_spec$formula.list[c(1, 3)],
  h_alpha.list = ps_spec$h_alpha.list[c(1, 3)],
  inv_link = ps_spec$inv_link,
  outcome = ps_spec$outcome
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

# Get outlier indices
f <- "Simulation_Results/EBMR_IPW_setting4-miss50-scenario9-2_13_n2000_replicate1000_test57.RDS"
res <- readRDS(f)
mu <- res[1, ]; se <- res[3, ]
w1_stored <- res["w.hat1", ]; w2_stored <- res["w.hat2", ]

valid <- !is.na(mu) & !is.na(se)
mu_v <- mu[valid]; se_v <- se[valid]
q1 <- quantile(mu_v, 0.25); q3 <- quantile(mu_v, 0.75); iqr <- q3 - q1
is_out_mu <- mu_v < (q1 - 3*iqr) | mu_v > (q3 + 3*iqr)
q1s <- quantile(se_v, 0.25); q3s <- quantile(se_v, 0.75); iqrs <- q3s - q1s
is_out_se <- se_v < (q1s - 3*iqrs) | se_v > (q3s + 3*iqrs)
is_out <- is_out_mu | is_out_se
outlier_idx <- which(valid)[is_out]
normal_idx <- which(valid)[!is_out]

check_reps <- c(head(outlier_idx, 5), head(normal_idx, 3))
cat("Checking reps:", paste(check_reps, collapse=", "), "\n\n")

for (rep_i in check_reps) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]
  is_outlier_rep <- rep_i %in% outlier_idx
  label <- if (is_outlier_rep) "OUTLIER" else "NORMAL"

  cat(sprintf("\n=== Rep %d [%s] (stored w1=%.4f, w2=%.4f) ===\n",
              rep_i, label, w1_stored[rep_i], w2_stored[rep_i]))

  tryCatch({
    ebmr <- EPS$new("y", subset_ps_spec, dat, W_func)
    J <- 2
    ps.matrix <- do.call(cbind, lapply(ebmr$ps_fit.list, function(f) f$fitted.values))

    h_x2 <- cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2, u1_u2 = dat$u1*dat$u2)
    d <- ps.matrix[, -1, drop = FALSE] - ps.matrix[, 1]
    h_x <- cbind(d, h_x2)
    h_dim <- ncol(h_x)
    r_vec <- as.numeric(dat$r)
    n <- nrow(dat)

    Phi_nu <- function(nu) {
      ps_nu <- as.vector(ps.matrix %*% nu)
      (r_vec / ps_nu - 1) * h_x
    }

    run_gmm <- function(init_val) {
      estimates <- init_val
      for (t in 1:200) {
        g_mat <- Phi_nu(estimates)
        W_hat <- tryCatch(
          solve(t(g_mat) %*% g_mat / n),
          error = function(e) diag(h_dim)
        )

        obj_fn <- function(nu) {
          g_m <- Phi_nu(nu)
          G_v <- matrix(colMeans(g_m), h_dim, 1)
          v <- as.numeric(crossprod(G_v, W_hat %*% G_v))
          if (!is.finite(v)) 1e8 else v
        }

        opt <- optim(estimates, obj_fn, method = "L-BFGS-B",
                     lower = rep(-Inf, J), upper = rep(Inf, J),
                     control = list(maxit = 1000))

        conv_err <- max(abs(opt$par - estimates))
        estimates <- opt$par
        if (conv_err < 1e-6) break
      }

      g_final <- Phi_nu(estimates)
      W_final <- tryCatch(solve(t(g_final) %*% g_final / n), error = function(e) diag(h_dim))
      G_final <- matrix(colMeans(g_final), h_dim, 1)
      obj_final <- as.numeric(crossprod(G_final, W_final %*% G_final))
      w_final <- estimates^2 / sum(estimates^2)
      list(nu = estimates, w = w_final, obj = obj_final, iter = t)
    }

    inits <- list(
      "uniform"  = c(0.5, 0.5),
      "model1"   = c(1, 0.001),
      "model2"   = c(0.001, 1)
    )

    for (nm in names(inits)) {
      r <- tryCatch(
        run_gmm(inits[[nm]]),
        error = function(e) {
          cat(sprintf("  init=%-10s ERROR: %s\n", nm, conditionMessage(e)))
          list(nu=c(NA,NA), w=c(NA,NA), obj=NA, iter=NA)
        }
      )
      if (!is.na(r$obj)) {
        cat(sprintf("  init=%-10s -> nu=[%7.4f,%7.4f], w=[%.4f,%.4f], obj=%.8f, iter=%d\n",
                    nm, r$nu[1], r$nu[2], r$w[1], r$w[2], r$obj, r$iter))
      }
    }

  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", conditionMessage(e)))
  })
}

cat("\nDone!\n")
