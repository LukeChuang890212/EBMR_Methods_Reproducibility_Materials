setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})

n_val <- 2000
n_reps <- 200
COND_THRESH <- 1e8
TRUST_RADIUS <- 2.0

run_test <- function(setting, miss, ps_spec_id, model_idx) {
  ps_spec <- get_ps_spec(ps_spec_id)
  data_file <- misspecified_model_all_data_file.list[[setting]][[miss]][[1]]
  all_data <- readRDS(data_file)
  mu_true <- get_mu_true(setting)

  formula_j <- ps_spec[["formula.list"]][[model_idx]]
  h_alpha_vars <- ps_spec[["h_alpha.list"]][[model_idx]]

  mu_vals <- rep(NA_real_, n_reps)
  se1_vals <- rep(NA_real_, n_reps)
  se2_vals <- rep(NA_real_, n_reps)
  alpha2_vals <- rep(NA_real_, n_reps)
  conv_vals <- rep(NA, n_reps)

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

      nr_inner <- function(start, W_hat) {
        alpha <- start
        for (k in 1:500) {
          g_mat <- Phi_alpha(alpha); G_vec <- colMeans(g_mat)
          Gamma_hat <- Gamma_fn(alpha)
          g_vec <- 2 * as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))
          if (max(abs(g_vec)) < 1e-8) break
          H_gn <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
          direction <- tryCatch(solve(H_gn, -g_vec), error = function(e) -g_vec)
          sn <- sqrt(sum(direction^2))
          if (sn > TRUST_RADIUS) direction <- direction * (TRUST_RADIUS / sn)
          obj_cur <- as.numeric(t(G_vec) %*% W_hat %*% G_vec)
          step <- 1.0; accepted <- FALSE
          for (ls in 1:30) {
            cand <- alpha + step * direction
            G_new <- colMeans(Phi_alpha(cand))
            obj_new <- as.numeric(t(G_new) %*% W_hat %*% G_new)
            if (is.finite(obj_new) && obj_new < obj_cur - 1e-4*step*sum(g_vec*direction)) {
              Gamma_c <- Gamma_fn(cand)
              H_c <- crossprod(Gamma_c, W_hat %*% Gamma_c)
              eig_c <- eigen(H_c, symmetric=TRUE, only.values=TRUE)$values
              cc <- max(eig_c)/max(min(eig_c), 1e-15)
              if (cc < COND_THRESH) { accepted <- TRUE; break }
            }
            step <- step * 0.5
          }
          if (accepted) { alpha <- cand } else { break }
        }
        alpha
      }

      estimates <- nr_inner(rep(0, p), diag(h_dim))
      best_grad <- Inf; best_est <- estimates; no_improve <- 0L; converged <- FALSE
      for (t in 1:200) {
        g_mat <- Phi_alpha(estimates); G_vec <- colMeans(g_mat)
        W_hat <- tryCatch(W_func(g_mat), error = function(e) diag(h_dim))
        Gamma_hat <- Gamma_fn(estimates)
        conv_grad <- 2*as.vector(crossprod(Gamma_hat, W_hat %*% G_vec))
        grad_norm <- max(abs(conv_grad))
        if (grad_norm < 1e-6) { converged <- TRUE; break }
        if (grad_norm < best_grad) { best_grad <- grad_norm; best_est <- estimates; no_improve <- 0L
        } else { no_improve <- no_improve+1L; if (no_improve >= 20L) { estimates <- best_est; break } }
        estimates <- nr_inner(estimates, W_hat)
      }

      pi_hat <- model_fn(estimates); mu_hat <- mean(r_vec*y_vec/pi_hat)
      g_mat <- Phi_alpha(estimates); G_vec <- colMeans(g_mat)
      W_hat <- tryCatch(W_func(g_mat), error = function(e) diag(h_dim))
      Gamma_hat <- Gamma_fn(estimates)
      H_mat <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
      H_inv <- tryCatch(solve(H_mat), error = function(e) NULL)
      if (!is.null(H_inv)) {
        GtW <- crossprod(Gamma_hat, W_hat)
        psi1 <- H_inv %*% GtW %*% t(g_mat)
        dmu <- colMeans(-r_vec*y_vec*(1-pi_hat)/pi_hat * design_mat)

        # SE1
        iid1 <- r_vec*y_vec/pi_hat - mu_hat - as.vector(t(dmu) %*% psi1)
        se1_vals[rep_i] <- sqrt(mean(iid1^2)/n)

        # SE2
        WG <- W_hat %*% G_vec
        gi_WG <- as.vector(g_mat %*% WG)
        Omega <- crossprod(g_mat)/n
        OmWG <- as.vector(Omega %*% WG)
        GtW_g <- GtW %*% t(g_mat)
        GtW_OmWG <- as.vector(GtW %*% OmWG)
        corr <- sweep(GtW_g * rep(gi_WG, each=p), 1, GtW_OmWG, "-")
        psi2 <- psi1 - H_inv %*% corr
        iid2 <- r_vec*y_vec/pi_hat - mu_hat - as.vector(t(dmu) %*% psi2)
        se2_vals[rep_i] <- sqrt(mean(iid2^2)/n)
      }
      mu_vals[rep_i] <- mu_hat; alpha2_vals[rep_i] <- estimates[2]; conv_vals[rep_i] <- converged
    }, error = function(e) NULL)
  }

  valid <- !is.na(mu_vals) & !is.na(se1_vals)
  esd <- sd(mu_vals[valid]); ese1 <- mean(se1_vals[valid]); ese2 <- mean(se2_vals[valid], na.rm=TRUE)
  cat(sprintf("\n=== %s, model %d ===\n", setting, model_idx))
  cat(sprintf("Valid: %d, Converged: %d/%d\n", sum(valid), sum(conv_vals[valid]), sum(valid)))
  cat(sprintf("Bias=%.4f, ESD=%.4f\n", mean(mu_vals[valid]) - mu_true, esd))
  cat(sprintf("ESE1=%.4f, ESE1/ESD=%.3f\n", ese1, ese1/esd))
  cat(sprintf("ESE2=%.4f, ESE2/ESD=%.3f\n", ese2, ese2/esd))
  cat(sprintf("alpha[2]: mean=%.3f, sd=%.3f, range=[%.3f, %.3f]\n",
      mean(alpha2_vals[valid]), sd(alpha2_vals[valid]),
      min(alpha2_vals[valid]), max(alpha2_vals[valid])))
}

cat(sprintf("Config: cond_threshold=%.0e, trust_radius=%.1f\n", COND_THRESH, TRUST_RADIUS))

# Setting4: all 3 models (scenario 9-3)
cat("\n========== Setting4 (scenario 9-3) ==========\n")
for (m in 1:3) {
  run_test("setting4", "miss50", "9-alt1", m)
}

# Setting3: model 2 (previous ESE2/ESD issue)
cat("\n========== Setting3 Model2 ==========\n")
run_test("setting3", "miss50", "9-alt1", 2)
