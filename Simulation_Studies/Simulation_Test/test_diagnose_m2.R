setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})

n_val <- 2000
n_reps <- 50  # fewer reps for diagnostics
COND_THRESH <- 1e8
TRUST_RADIUS <- 2.0

diagnose <- function(setting, miss, ps_spec_id, model_idx) {
  ps_spec <- get_ps_spec(ps_spec_id)
  data_file <- misspecified_model_all_data_file.list[[setting]][[miss]][[1]]
  all_data <- readRDS(data_file)
  mu_true <- get_mu_true(setting)

  formula_j <- ps_spec[["formula.list"]][[model_idx]]
  h_alpha_vars <- ps_spec[["h_alpha.list"]][[model_idx]]

  cat(sprintf("\n========== %s, model %d ==========\n", setting, model_idx))
  cat(sprintf("Formula: %s\n", deparse(formula_j)))
  cat(sprintf("h_alpha: %s\n", paste(h_alpha_vars, collapse=", ")))
  cat(sprintf("param_dim: %d, esteq_dim: %d (overid=%d)\n",
      ncol(model.matrix(formula_j, data=all_data[1:10,])),
      1 + length(h_alpha_vars),
      1 + length(h_alpha_vars) - ncol(model.matrix(formula_j, data=all_data[1:10,]))))

  for (rep_i in 1:n_reps) {
    dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
    tryCatch({
      r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
      design_mat <- model.matrix(formula_j, data=dat)
      h_full <- cbind(1, as.matrix(dat[, h_alpha_vars, drop=FALSE]))
      n <- n_val; p <- ncol(design_mat); h_dim <- ncol(h_full)

      model_fn <- function(alpha) { eta <- as.vector(design_mat %*% alpha); 1/(1+exp(-eta)) }
      Phi_alpha <- function(alpha) { pi_hat <- model_fn(alpha); (r_vec/pi_hat - 1) * h_full }
      Gamma_fn <- function(alpha) {
        pi_hat <- model_fn(alpha)
        crossprod(h_full * (-r_vec*(1-pi_hat)/pi_hat), design_mat) / n
      }
      W_func <- function(g_mat) solve(crossprod(g_mat) / nrow(g_mat))

      # Run baseline L-BFGS-B
      est_bl <- rep(0, p)
      for (t in 1:200) {
        g_mat <- Phi_alpha(est_bl); G_vec <- colMeans(g_mat)
        if (t == 1) { W_hat <- diag(h_dim)
        } else { W_hat <- tryCatch(W_func(g_mat), error = function(e) diag(h_dim)) }
        Gamma_hat <- Gamma_fn(est_bl)
        conv_grad <- 2*as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))
        if (max(abs(conv_grad)) < 1e-6) break
        obj_fn <- function(alpha) {
          G <- colMeans(Phi_alpha(alpha))
          as.numeric(t(G) %*% W_hat %*% G)
        }
        opt <- optim(est_bl, obj_fn, method="L-BFGS-B", control=list(maxit=5000))
        est_bl <- opt$par
      }

      pi_bl <- model_fn(est_bl)
      Gamma_bl <- Gamma_fn(est_bl)
      g_mat_bl <- Phi_alpha(est_bl)
      W_bl <- tryCatch(W_func(g_mat_bl), error=function(e) diag(h_dim))
      H_bl <- crossprod(Gamma_bl, W_bl %*% Gamma_bl)
      eig_bl <- eigen(H_bl, symmetric=TRUE, only.values=TRUE)$values
      cond_bl <- max(eig_bl)/max(min(eig_bl), 1e-15)

      if (rep_i <= 10 || cond_bl > 1e6) {
        cat(sprintf("  Rep %3d: alpha=[%s], max|alpha|=%.2f, cond(H)=%.2e, min(pi)=%.4f, max(pi)=%.4f, mean(pi*(1-pi))=%.4f\n",
            rep_i,
            paste(sprintf("%.2f", est_bl), collapse=","),
            max(abs(est_bl)),
            cond_bl,
            min(pi_bl), max(pi_bl),
            mean(pi_bl*(1-pi_bl))))
      }
    }, error = function(e) cat(sprintf("  Rep %3d: ERROR - %s\n", rep_i, conditionMessage(e))))
  }
}

diagnose("setting4", "miss50", "9-alt1", 2)
diagnose("setting4", "miss50", "9-alt1", 3)
diagnose("setting3", "miss50", "9-alt1", 2)
