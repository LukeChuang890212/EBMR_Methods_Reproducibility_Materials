setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
source("Simulation.r")
library(EBMRalgorithmFast4)
library(parallel)
library(foreach)
library(doSNOW)

mu_true <- get_mu_true("setting3")

ps_spec <- get_ps_spec("9-alt1")
single_ps_spec <- list(
  formula.list = ps_spec[["formula.list"]][2],
  h_alpha.list = ps_spec[["h_alpha.list"]][2],
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]]
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

run_diagnostic <- function(n_val, n_reps, mu_true) {
  cat(sprintf("\n###############################################\n"))
  cat(sprintf("### n = %d, %d reps (setting3.B1, model 2) ###\n", n_val, n_reps))
  cat(sprintf("###############################################\n"))
  cat(sprintf("mu_true = %.4f\n", mu_true))

  cores <- max(1, detectCores() - 2)
  cl <- makeCluster(cores)
  registerDoSNOW(cl)
  pb <- txtProgressBar(max = n_reps, style = 3)
  opts <- list(progress = function(i) setTxtProgressBar(pb, i))

  clusterExport(cl, c("setting3.B1", "single_ps_spec", "W_func"), envir = environment())

  results_raw <- foreach(rep_i = 1:n_reps,
                         .combine = cbind,
                         .options.snow = opts,
                         .packages = c("EBMRalgorithmFast4", "numDeriv"),
                         .errorhandling = "pass") %dopar% {

    set.seed(12345 + rep_i)
    dat <- setting3.B1(n_val)
    res_vec <- rep(NA, 8)

    tryCatch({
      ebmr <- EBMRAlgorithmFast4[["new"]]("y", single_ps_spec, dat, W_func)
      res <- ebmr[["EBMR_IPW"]](
        h_nu = function(dat) cbind(u1 = dat[["u1"]], u2 = dat[["u2"]], z1 = dat[["z1"]], z2 = dat[["z2"]]),
        se.fit = TRUE, type = "HT"
      )

      res_vec[1] <- res[["mu_ipw"]]
      res_vec[2] <- 0
      res_vec[3] <- res[["se_ipw"]]
      res_vec[4] <- 0

      gmm_fit <- ebmr[["ps_fit.list"]][[1]][["gmm_fit"]]
      Gamma.hat <- gmm_fit[["Gamma.hat"]]
      W.hat <- gmm_fit[["W.hat"]]
      g.matrix <- gmm_fit[["g.matrix"]]
      psi_alpha <- gmm_fit[["psi"]]
      eta_s <- gmm_fit[["eta_s"]]
      Q_full <- gmm_fit[["Q"]]

      n <- nrow(g.matrix)
      param_dim <- nrow(psi_alpha)

      GtW <- crossprod(Gamma.hat, W.hat)
      GtWG <- GtW %*% Gamma.hat
      H_0n <- tryCatch(-solve(GtWG), error = function(e) -solve(GtWG + 1e-8 * diag(nrow(GtWG))))
      t_g <- t(g.matrix)
      H_1n <- GtW %*% (t_g - eta_s)

      ps_fit <- ebmr[["ps_fit.list"]][[1]]
      pi_hat <- ps_fit[["fitted.values"]]
      r_vec <- dat[["r"]]
      y_vec <- dat[["y"]]
      design_mat <- ps_fit[["design_matrix"]]
      link_type <- ps_fit[["link_type"]]

      eps_mix <- 0.01
      scale_mix <- 0.99
      f_hat <- (pi_hat - eps_mix) / scale_mix

      ry_ps_inv2 <- as.vector(r_vec * y_vec * (pi_hat^(-2)))
      if (link_type == "logistic_complement") {
        dot_pi <- -scale_mix * design_mat * (f_hat * (1 - f_hat))
      } else if (link_type == "logistic") {
        dot_pi <- scale_mix * design_mat * (f_hat * (1 - f_hat))
      }
      H_alpha.w <- colMeans(dot_pi * ry_ps_inv2)
      base_term <- as.vector(r_vec / pi_hat * y_vec)

      psi_simple <- H_0n %*% GtW %*% t_g
      mu_iid_simple <- base_term - as.vector(t(H_alpha.w) %*% psi_simple)
      res_vec[5] <- sqrt(var(mu_iid_simple) / n)

      psi_H1only <- Q_full %*% H_1n
      mu_iid_H1 <- base_term - as.vector(t(H_alpha.w) %*% psi_H1only)
      res_vec[6] <- sqrt(var(mu_iid_H1) / n)

      W_eta <- as.vector(W.hat %*% eta_s)
      tGamma_W_eta <- as.vector(GtW %*% eta_s)
      GtW_g <- GtW %*% t_g
      g_W_eta <- as.vector(g.matrix %*% W_eta)
      H_2n_3 <- matrix(tGamma_W_eta, param_dim, n) - GtW_g * matrix(g_W_eta, param_dim, n, byrow = TRUE)
      Q_inv <- solve(Q_full)
      H_sum <- Q_inv %*% psi_alpha
      H_2n_2 <- H_sum - H_1n - H_2n_3

      psi_noMn <- H_0n %*% (H_1n + H_2n_2 + H_2n_3)
      mu_iid_noMn <- base_term - as.vector(t(H_alpha.w) %*% psi_noMn)
      res_vec[7] <- sqrt(var(mu_iid_noMn) / n)

      psi_noH2n3 <- Q_full %*% (H_1n + H_2n_2)
      mu_iid_noH2n3 <- base_term - as.vector(t(H_alpha.w) %*% psi_noH2n3)
      res_vec[8] <- sqrt(var(mu_iid_noH2n3) / n)

    }, error = function(e) {})

    res_vec
  }

  close(pb)
  stopCluster(cl)

  # Handle errors returned as list elements
  if (is.list(results_raw)) {
    results_raw <- do.call(cbind, lapply(results_raw, function(x) if (is.numeric(x)) x else rep(NA, 8)))
  }
  rownames(results_raw) <- c("mu_ipw", "mu_ipw.true", "se_ipw", "se_ipw.true",
                             "se_simple", "se_H1only", "se_noMn", "se_noH2n3")

  n_err <- sum(is.na(results_raw[1, ]))
  cat(sprintf("\nErrors: %d / %d\n", n_err, n_reps))

  # Apply clean_sim_result
  cleaned <- clean_sim_result(results_raw[1:4, ], multiplier = 3, verbose = TRUE)

  na_mask <- apply(results_raw[1:4, , drop = FALSE], 2, function(x) !any(is.na(x)))
  results_nona <- results_raw[, na_mask, drop = FALSE]
  mu_ipw <- results_nona[1, ]
  Q1 <- quantile(mu_ipw, 0.25)
  Q3 <- quantile(mu_ipw, 0.75)
  IQR_value <- Q3 - Q1
  is_outlier <- mu_ipw < (Q1 - 3 * IQR_value) | mu_ipw > (Q3 + 3 * IQR_value)
  max_remove <- floor(0.01 * length(mu_ipw))
  if (sum(is_outlier) > max_remove && max_remove > 0) {
    dist_from_median <- abs(mu_ipw - median(mu_ipw))
    outlier_idx <- which(is_outlier)
    keep_idx <- outlier_idx[order(dist_from_median[outlier_idx], decreasing = TRUE)[(max_remove + 1):length(outlier_idx)]]
    is_outlier[keep_idx] <- FALSE
  }
  rc <- results_nona[, !is_outlier, drop = FALSE]

  cat(sprintf("\nKept: %d (removed %d NA, %d outlier)\n", ncol(rc), sum(!na_mask), sum(is_outlier)))

  # Before cleaning
  mu_all <- results_nona[1, ]
  esd_all <- sd(mu_all)
  cat(sprintf("\n--- BEFORE cleaning ---\n"))
  cat(sprintf("N: %d, Bias: %.4f, ESD: %.4f\n", length(mu_all), mean(mu_all) - mu_true, esd_all))
  cat(sprintf("%-15s %10s\n", "SE variant", "Mean ESE"))
  for (nm in c("se_ipw", "se_simple", "se_H1only", "se_noMn", "se_noH2n3")) {
    idx <- which(rownames(results_nona) == nm)
    cat(sprintf("%-15s %10.4f\n", nm, mean(results_nona[idx, ])))
  }

  # After cleaning
  mu_c <- rc[1, ]
  esd_c <- sd(mu_c)
  cat(sprintf("\n--- AFTER cleaning ---\n"))
  cat(sprintf("N: %d, Bias: %.4f, ESD: %.4f\n\n", length(mu_c), mean(mu_c) - mu_true, esd_c))
  cat(sprintf("%-15s %10s %10s\n", "SE variant", "Mean ESE", "Mean/ESD"))
  for (nm in c("se_ipw", "se_simple", "se_H1only", "se_noMn", "se_noH2n3")) {
    idx <- which(rownames(rc) == nm)
    se_vals <- rc[idx, ]
    cat(sprintf("%-15s %10.4f %10.4f\n", nm, mean(se_vals), mean(se_vals)/esd_c))
  }

  ci_l <- rc[1, ] - 1.96 * rc[3, ]
  ci_u <- rc[1, ] + 1.96 * rc[3, ]
  cp <- mean(mu_true >= ci_l & mu_true <= ci_u)
  cat(sprintf("\nCP (se_ipw): %.4f\n", cp))

  ci_l2 <- rc[1, ] - 1.96 * rc[8, ]
  ci_u2 <- rc[1, ] + 1.96 * rc[8, ]
  cp2 <- mean(mu_true >= ci_l2 & mu_true <= ci_u2)
  cat(sprintf("CP (se_noH2n3): %.4f\n", cp2))
}

run_diagnostic(2000, 500, mu_true)
run_diagnostic(3000, 500, mu_true)
