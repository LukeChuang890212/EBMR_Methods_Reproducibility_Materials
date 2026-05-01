setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})

n_val <- 2000
n_reps <- 200
COND_THRESH <- 1e8
TRUST_RADIUS <- 2.0

ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
mu_true <- get_mu_true("setting4")

formula_j <- ps_spec[["formula.list"]][[3]]
h_alpha_vars <- ps_spec[["h_alpha.list"]][[3]]

nr_inner_v2 <- function(start, W_hat, design_mat, r_vec, h_full, n, p,
                         model_fn, Phi_alpha, Gamma_fn) {
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

    # Try NR direction first, then fallback to steepest descent
    accepted <- FALSE
    for (try_dir in 1:2) {
      if (try_dir == 2) {
        # Fallback: steepest descent with trust region
        direction <- -g_vec
        sn <- sqrt(sum(direction^2))
        if (sn > TRUST_RADIUS) direction <- direction * (TRUST_RADIUS / sn)
      }
      step <- 1.0
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
      if (accepted) break
    }
    if (accepted) { alpha <- cand } else { break }
  }
  alpha
}

run_variant <- function(label, inner_fn) {
  mu_vals <- rep(NA_real_, n_reps)
  se2_vals <- rep(NA_real_, n_reps)
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

      estimates <- inner_fn(rep(0, p), diag(h_dim), design_mat, r_vec, h_full, n, p,
                            model_fn, Phi_alpha, Gamma_fn)
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
        estimates <- inner_fn(estimates, W_hat, design_mat, r_vec, h_full, n, p,
                              model_fn, Phi_alpha, Gamma_fn)
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
      mu_vals[rep_i] <- mu_hat; conv_vals[rep_i] <- converged
    }, error = function(e) NULL)
  }

  valid <- !is.na(mu_vals) & !is.na(se2_vals)
  esd <- sd(mu_vals[valid]); ese2 <- mean(se2_vals[valid])
  cat(sprintf("%s: ESE2/ESD=%.3f, Bias=%.4f, conv=%d/%d\n",
      label, ese2/esd, mean(mu_vals[valid]) - mu_true, sum(conv_vals[valid]), sum(valid)))
}

# Original: NR only
nr_inner_v1 <- function(start, W_hat, design_mat, r_vec, h_full, n, p,
                         model_fn, Phi_alpha, Gamma_fn) {
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

cat("Setting4, model 3:\n")
run_variant("v1 (NR only)", nr_inner_v1)
run_variant("v2 (NR + GD fallback)", nr_inner_v2)
