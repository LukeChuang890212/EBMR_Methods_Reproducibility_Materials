setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

n_val <- 2000
n_reps <- 200

ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
mu_true <- get_mu_true("setting4")

run_all_se <- function(model_idx, optimizer_name) {
  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL),
    optimizer = optimizer_name
  )
  W <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

  h_alpha_vars <- ps_spec[["h_alpha.list"]][[model_idx]]

  mu_vals <- rep(NA_real_, n_reps)
  # 4 SE variants:
  se_A <- rep(NA_real_, n_reps)  # Simple SE1: -(GtWG)^-1 GtW g_i, then r*y/pi - mu - dmu'*psi
  se_B <- rep(NA_real_, n_reps)  # Simple SE2: SE1 + CUE correction
  se_C <- rep(NA_real_, n_reps)  # Package psi1 -> mu_iid
  se_D <- rep(NA_real_, n_reps)  # Package psi2 -> mu_iid

  for (rep_i in 1:n_reps) {
    dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
    tryCatch({
      ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W)
      r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
      pi_hat <- ebmr$ps_fit.list[[1]]$fitted.values
      mu_hat <- mean(r_vec * y_vec / pi_hat)
      n <- n_val
      gmm_fit <- ebmr$ps_fit.list[[1]]$gmm_fit
      design_mat <- ebmr$ps_fit.list[[1]]$design_matrix
      p <- ncol(design_mat)

      link_type <- ebmr$ps_fit.list[[1]]$link_type
      if (!is.null(link_type) && link_type == "logistic_complement") {
        dot_pi <- -design_mat * (pi_hat * (1 - pi_hat))
      } else {
        dot_pi <- design_mat * (pi_hat * (1 - pi_hat))
      }
      ry_ps_inv2 <- as.vector(r_vec * y_vec * (pi_hat^(-2)))
      H_alpha <- colMeans(dot_pi * ry_ps_inv2)

      Gamma_hat <- gmm_fit$Gamma.hat
      g_mat <- gmm_fit$g.matrix
      W_hat <- gmm_fit$W.hat
      G_vec <- gmm_fit$eta_s
      GtW <- crossprod(Gamma_hat, W_hat)
      GtWG <- GtW %*% Gamma_hat
      H_inv <- tryCatch(solve(GtWG), error = function(e) NULL)

      if (!is.null(H_inv)) {
        # Simple psi for alpha
        psi1_s <- H_inv %*% GtW %*% t(g_mat)  # k x n
        dmu <- colMeans(-r_vec*y_vec*(1-pi_hat)/pi_hat * design_mat)

        # SE_A: simple SE1
        iid_A <- r_vec*y_vec/pi_hat - mu_hat - as.vector(t(dmu) %*% psi1_s)
        se_A[rep_i] <- sqrt(mean(iid_A^2)/n)

        # SE_B: simple SE2 (CUE correction)
        WG <- W_hat %*% G_vec
        gi_WG <- as.vector(g_mat %*% WG)
        Omega <- crossprod(g_mat)/n
        OmWG <- as.vector(Omega %*% WG)
        GtW_g <- GtW %*% t(g_mat)
        GtW_OmWG <- as.vector(GtW %*% OmWG)
        corr <- sweep(GtW_g * rep(gi_WG, each=p), 1, GtW_OmWG, "-")
        psi2_s <- psi1_s - H_inv %*% corr
        iid_B <- r_vec*y_vec/pi_hat - mu_hat - as.vector(t(dmu) %*% psi2_s)
        se_B[rep_i] <- sqrt(mean(iid_B^2)/n)
      }

      # SE_C: package psi1
      psi1_pkg <- gmm_fit$psi
      if (!any(is.na(psi1_pkg))) {
        mu_iid_C <- as.vector(t(r_vec/pi_hat*y_vec) - t(H_alpha) %*% psi1_pkg)
        se_C[rep_i] <- sqrt(var(mu_iid_C)/n)
      }

      # SE_D: package psi2
      psi2_pkg <- gmm_fit$psi2
      if (!is.null(psi2_pkg) && is.matrix(psi2_pkg) && !any(is.na(psi2_pkg))) {
        mu_iid_D <- as.vector(t(r_vec/pi_hat*y_vec) - t(H_alpha) %*% psi2_pkg)
        se_D[rep_i] <- sqrt(var(mu_iid_D)/n)
      }

      mu_vals[rep_i] <- mu_hat
    }, error = function(e) NULL)
  }

  valid <- !is.na(mu_vals)
  cat(sprintf("\n%s, model %d, %s:\n", "setting4", model_idx, optimizer_name))
  cat(sprintf("  Bias=%.4f, ESD=%.4f\n", mean(mu_vals[valid]) - mu_true, sd(mu_vals[valid])))
  for (nm in c("A","B","C","D")) {
    sv <- get(paste0("se_", nm))
    v <- valid & !is.na(sv)
    label <- c(A="Simple SE1", B="Simple SE2 (CUE)", C="Package psi1", D="Package psi2")[nm]
    if (sum(v) > 0) {
      cat(sprintf("  %-20s: ESE/ESD=%.3f  (valid=%d)\n", label, mean(sv[v])/sd(mu_vals[v]), sum(v)))
    }
  }
}

for (m in c(3, 2)) {
  for (opt in c("L-BFGS-B", "constrained_nr")) {
    run_all_se(m, opt)
  }
}
