## Investigate NAs in setting3, scenario 9-2, model 2, n=500, miss30
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

sim_file <- "Simulation_Results/EBMR_IPW_setting3-miss30-scenario9-2_2_n500_replicate1000_test59.RDS"
sim_result <- readRDS(sim_file)

# Find NA reps
mu_ipw <- sim_result[1, ]
se_ipw <- sim_result[3, ]
na_reps <- which(is.na(mu_ipw) | is.na(se_ipw))
cat(sprintf("NA reps (%d): %s\n\n", length(na_reps), paste(na_reps, collapse=", ")))

# Check: is it mu that's NA, se, or both?
mu_na <- which(is.na(mu_ipw))
se_na <- which(is.na(se_ipw))
cat(sprintf("mu_ipw NA: %d reps\n", length(mu_na)))
cat(sprintf("se_ipw NA: %d reps\n", length(se_na)))
cat(sprintf("Both NA: %d reps\n", length(intersect(mu_na, se_na))))
cat(sprintf("Only se NA: %d reps\n\n", length(setdiff(se_na, mu_na))))

# Refit NA reps to diagnose
ps_spec <- get_ps_spec("9")
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
n_val <- 500
all_data <- readRDS("Simulation_Data/Setting3.B2_n500_replicate1000.RDS")

m_idx <- 2
single_ps <- list(
  formula.list = list(ps_spec[["formula.list"]][[m_idx]]),
  h_alpha.list = list(ps_spec[["h_alpha.list"]][[m_idx]]),
  inv_link = ps_spec[["inv_link"]], outcome = ps_spec[["outcome"]]
)

cat("=== Refitting NA reps ===\n\n")
for (rep_i in na_reps) {
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]

  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
    pf <- ebmr$ps_fit.list[[1]]
    ps <- pf$fitted.values
    alpha <- pf$coefficients
    se_alpha <- pf$se

    # Check GMM fit details
    gmm <- pf$gmm_fit
    converged <- gmm$opt$converged
    obj <- gmm$opt$objective
    grad <- gmm$opt$final_grad_norm

    # Condition number
    g_mat <- (dat$r / ps - 1) * pf$h_x
    W_hat <- tryCatch(solve(crossprod(g_mat) / n_val), error = function(e) NULL)

    if (is.null(W_hat)) {
      cat(sprintf("Rep %3d: W_hat SINGULAR, alpha=(%s), PS=[%.4f,%.4f]\n",
          rep_i, paste(round(alpha, 3), collapse=", "), min(ps), max(ps)))
      next
    }

    Gamma <- crossprod(pf$h_x * (dat$r * (1-ps) / ps), pf$design_matrix) / n_val
    H <- crossprod(Gamma, W_hat %*% Gamma)
    eig <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    cond <- max(eig) / max(min(eig), 1e-15)

    # Check H_mat for SE (the CUE H_mat)
    h_dim <- ncol(pf$h_x); k <- ncol(pf$design_matrix)
    GtW <- crossprod(Gamma, W_hat)
    G_vec <- colMeans(g_mat)
    eta_s <- matrix(G_vec, h_dim, 1)
    GW_row <- as.vector(t(eta_s) %*% W_hat)
    # Quick check of H_mat conditioning
    cf <- dat$r * (1 - ps) / ps
    dc_X <- cf * pf$design_matrix
    R_mat <- matrix(0, h_dim * k, k)
    for (j in 1:k) {
      block <- crossprod(dc_X, pf$design_matrix[, j] * pf$h_x) / n_val
      for (l in 1:h_dim) R_mat[(l-1)*k + j, ] <- block[, l]
    }
    S_mat <- matrix(0, h_dim^2, k)
    for (j in 1:k) {
      dg_j <- cf * pf$design_matrix[, j] * pf$h_x
      cross <- (crossprod(dg_j, g_mat) + crossprod(g_mat, dg_j)) / n_val
      S_mat[, j] <- as.vector(cross)
    }
    GtWG <- crossprod(Gamma, W_hat %*% Gamma)
    GW_kron_Ik <- kronecker(matrix(GW_row, 1, h_dim), diag(k))
    GW_kron_GtW <- kronecker(matrix(GW_row, 1, h_dim), GtW)
    H_mat <- GtWG + GW_kron_Ik %*% R_mat - GW_kron_GtW %*% S_mat
    H_ev <- eigen(H_mat, symmetric = FALSE, only.values = TRUE)$values
    H_cond <- max(Mod(H_ev)) / max(min(Mod(H_ev)), 1e-15)

    cat(sprintf("Rep %3d: conv=%s, obj=%.1e, grad=%.1e, cond=%.1e, H_cond=%.1e\n",
        rep_i, converged, obj, grad, cond, H_cond))
    cat(sprintf("         alpha=(%s)\n", paste(round(alpha, 3), collapse=", ")))
    cat(sprintf("         SE=(%s)\n", paste(if(any(is.na(se_alpha))) "NA" else round(se_alpha, 3), collapse=", ")))
    cat(sprintf("         PS=[%.4f,%.4f], >0.99:%d, <0.01:%d\n",
        min(ps), max(ps), sum(ps > 0.99), sum(ps < 0.01)))

  }, error = function(e) {
    cat(sprintf("Rep %3d: ERROR: %s\n", rep_i, e$message))
  })
}
