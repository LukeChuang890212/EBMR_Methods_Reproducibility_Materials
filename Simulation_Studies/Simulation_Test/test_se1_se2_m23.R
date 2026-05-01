setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
  library(EBMRalgorithmFast4)
})

W_func  <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_nu_fn <- function(dat) cbind(u1=dat[["u1"]], u2=dat[["u2"]], z1=dat[["z1"]], z2=dat[["z2"]], u1_u2=dat[["u1"]]*dat[["u2"]])
n_val   <- 2000
n_reps  <- 200

# Scenario 9-3: ps_spec "9-alt1", models 2+3
ps_spec <- get_ps_spec("9-alt1")
ps_sub <- list(
  formula.list = ps_spec[["formula.list"]][2:3],
  h_alpha.list = ps_spec[["h_alpha.list"]][2:3],
  inv_link     = ps_spec[["inv_link"]],
  outcome      = ps_spec[["outcome"]]
)

# Data from setting4.B1
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
mu_true <- get_mu_true("setting4")

cat(sprintf("mu.true = %.4f\n", mu_true))
cat(sprintf("=== se1 vs se2, %d reps, models 2+3, scenario 9-3, setting4.B1 ===\n\n", n_reps))

mu_vals  <- rep(NA_real_, n_reps)
se1_vals <- rep(NA_real_, n_reps)
se2_vals <- rep(NA_real_, n_reps)
w2_vals  <- rep(NA_real_, n_reps)

n_err <- 0
for (i in 1:n_reps) {
  dat <- all_data[((i-1)*n_val + 1):(i*n_val), ]
  tryCatch({
    ebmr <- EBMRAlgorithmFast4[["new"]]("y", ps_sub, dat, W_func)
    res  <- ebmr[["EBMR_IPW"]](h_nu=h_nu_fn, type="HT", se.fit=TRUE)
    mu_vals[i]  <- res[["mu_ipw"]]
    se1_vals[i] <- res[["se_ipw"]]
    w_hat_ <- res[["w.hat"]]
    w2_vals[i] <- w_hat_[1]

    # --- se2-based se_ipw ---
    se2_vals[i] <- tryCatch({
    J2 <- 2
    r_   <- dat[["r"]]
    y_   <- dat[["y"]]
    nu_hat_ <- res[["nu.hat"]]
    ps_mat_ <- res[["ps.matrix"]]
    ens_ps_ <- as.vector(ps_mat_ %*% w_hat_)
    alpha_dim_ <- sapply(1:J2, function(j) length(ebmr[["ps_fit.list"]][[j]][["gmm_fit"]][["estimates"]]))

    # dot_pi
    dot_pi_full <- matrix(0, n_val, sum(alpha_dim_))
    for (j2 in 1:J2) {
      pi_j   <- ebmr[["ps_fit.list"]][[j2]][["fitted.values"]]
      des_j  <- ebmr[["ps_fit.list"]][[j2]][["design_matrix"]]
      lt_j   <- ebmr[["ps_fit.list"]][[j2]][["link_type"]]
      cols_j <- (sum(alpha_dim_[0:(j2-1)])+1):sum(alpha_dim_[1:j2])
      if (lt_j == "logistic") {
        dot_pi_full[, cols_j] <- des_j * (pi_j * (1 - pi_j))
      } else {
        dot_pi_full[, cols_j] <- -des_j * (pi_j * (1 - pi_j))
      }
    }

    ry_ps2_ <- as.vector(r_ * y_ * (ens_ps_^(-2)))
    H_alpha_ <- colMeans(t(t(dot_pi_full) * rep(w_hat_, alpha_dim_)) * ry_ps2_)

    # psi2 from each PS model
    psi2_list <- lapply(1:J2, function(j) ebmr[["ps_fit.list"]][[j]][["gmm_fit"]][["psi2"]])
    psi2_ok <- all(sapply(psi2_list, function(p) is.matrix(p) && !any(is.na(p))))

    if (psi2_ok) {
      psi2_alpha <- do.call(rbind, psi2_list)
      iid_ <- as.vector(r_ / ens_ps_ * y_) - as.vector(t(H_alpha_) %*% psi2_alpha)

      # nu adjustment (J > 1)
      if (J2 > 1) {
        dot_W_ <- function(nu) {
          nu <- as.vector(nu)
          (diag(2*nu)*sum(nu^2) - 2*nu %*% t(nu^2)) / (sum(nu^2)^2)
        }
        dot_W_nu_ <- dot_W_(nu_hat_)
        wH_nu_ <- colMeans(ps_mat_ %*% t(dot_W_nu_) * ry_ps2_)

        ens_fit <- res[["ensemble_fit"]]
        psi_nu_ <- ens_fit[["gmm_fit"]][["psi"]]
        Gamma_nu_ <- ens_fit[["gmm_fit"]][["Gamma.hat"]]
        W_nu_ <- ens_fit[["gmm_fit"]][["W.hat"]]

        h_nu_mat_ <- ens_fit[["h_x"]]  # augmented with intercept (5 cols, not 4)
        ps_nu_ <- as.vector(ps_mat_ %*% nu_hat_)
        r_ps2_nu_ <- as.vector(-r_ * (ps_nu_^(-2)))
        Phi_nu_alpha_ <- crossprod(cbind(h_nu_mat_) * r_ps2_nu_,
                                    t(t(dot_pi_full) * rep(nu_hat_, alpha_dim_))) / n_val
        GtW_nu_ <- crossprod(Gamma_nu_, W_nu_)
        dot_nu_ <- -solve(GtW_nu_ %*% Gamma_nu_) %*% GtW_nu_ %*% Phi_nu_alpha_

        iid_ <- iid_ - as.vector(
          (t(wH_nu_) %*% dot_nu_) %*% psi2_alpha + t(wH_nu_) %*% psi_nu_)
      }
      sqrt(var(iid_) / n_val)
    } else { NA_real_ }
    }, error = function(e) {
      if (i <= 3) cat(sprintf("  se2 error rep %d: %s\n", i, conditionMessage(e)))
      NA_real_
    })
  }, error = function(e) {
    n_err <<- n_err + 1
    if (n_err <= 5) cat(sprintf("  Error rep %d: %s\n", i, conditionMessage(e)))
  })
  if (i %% 50 == 0) cat(sprintf("  %d / %d done (errors: %d)\n", i, n_reps, n_err))
}

cat(sprintf("\nErrors: %d / %d\n", n_err, n_reps))

valid <- !is.na(mu_vals)
mu_v  <- mu_vals[valid]
se1_v <- se1_vals[valid]
se2_v <- se2_vals[valid]

esd  <- sd(mu_v)
ese1 <- mean(se1_v)
ese2 <- mean(se2_v, na.rm = TRUE)

cat(sprintf("\n=== RESULTS (models 2+3, scenario 9-3, setting4.B1) ===\n"))
cat(sprintf("Valid reps: %d\n", sum(valid)))
cat(sprintf("Bias: %.4f\n", mean(mu_v) - mu_true))
cat(sprintf("ESD:  %.4f\n", esd))
cat(sprintf("ESE1 (original):  %.4f  (ESE1/ESD = %.3f)\n", ese1, ese1/esd))
cat(sprintf("ESE2 (corrected): %.4f  (ESE2/ESD = %.3f)\n", ese2, ese2/esd))

ratio <- se2_v / se1_v
cat(sprintf("\nSE2/SE1 ratio: Mean=%.3f, Median=%.3f\n",
    mean(ratio, na.rm=TRUE), median(ratio, na.rm=TRUE)))

ci1_l <- mu_v - 1.96 * se1_v; ci1_u <- mu_v + 1.96 * se1_v
cp1   <- mean(mu_true >= ci1_l & mu_true <= ci1_u)
ci2_l <- mu_v - 1.96 * se2_v; ci2_u <- mu_v + 1.96 * se2_v
cp2   <- mean(mu_true >= ci2_l & mu_true <= ci2_u, na.rm = TRUE)

cat(sprintf("\nCP (se1): %.3f\n", cp1))
cat(sprintf("CP (se2): %.3f\n", cp2))
cat(sprintf("\nw2 (model 2 weight): mean=%.3f, range=[%.3f, %.3f]\n",
    mean(w2_vals[valid]), min(w2_vals[valid]), max(w2_vals[valid])))
