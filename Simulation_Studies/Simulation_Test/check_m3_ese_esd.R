setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
  devtools::load_all("../EBMRalgorithmFast4")
})

W_func  <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_nu_fn <- function(dat) cbind(u1=dat[["u1"]], u2=dat[["u2"]], z1=dat[["z1"]], z2=dat[["z2"]], u1_u2=dat[["u1"]]*dat[["u2"]])
n_val   <- 2000
n_reps  <- 200

ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
mu_true <- get_mu_true("setting4")

ps_sub3 <- list(
  formula.list = ps_spec[["formula.list"]][3],
  h_alpha.list = ps_spec[["h_alpha.list"]][3],
  inv_link     = ps_spec[["inv_link"]],
  outcome      = ps_spec[["outcome"]]
)

mu_vals   <- rep(NA_real_, n_reps)
se1_vals  <- rep(NA_real_, n_reps)
se2_vals  <- rep(NA_real_, n_reps)
conv_vals <- rep(NA, n_reps)
grad_vals <- rep(NA_real_, n_reps)
obj_vals  <- rep(NA_real_, n_reps)
iter_vals <- rep(NA_integer_, n_reps)
alpha1_vals <- rep(NA_real_, n_reps)
alpha2_vals <- rep(NA_real_, n_reps)

for (i in 1:n_reps) {
  dat <- all_data[((i-1)*n_val + 1):(i*n_val), ]
  tryCatch({
    ebmr <- EBMRAlgorithmFast4[["new"]]("y", ps_sub3, dat, W_func)
    res  <- ebmr[["EBMR_IPW"]](h_nu=h_nu_fn, type="HT", se.fit=TRUE)
    mu_vals[i]  <- res[["mu_ipw"]]
    se1_vals[i] <- res[["se_ipw"]]
    gf <- ebmr[["ps_fit.list"]][[1]][["gmm_fit"]]
    conv_vals[i] <- gf[["opt"]][["converged"]]
    grad_vals[i] <- gf[["opt"]][["final_grad_norm"]]
    obj_vals[i]  <- gf[["opt"]][["objective"]]
    iter_vals[i] <- gf[["opt"]][["iterations"]]
    alpha1_vals[i] <- gf[["estimates"]][1]
    alpha2_vals[i] <- gf[["estimates"]][2]

    # se2-based se_ipw: replace psi1 with psi2 in influence function (J=1, no nu adjustment)
    psi2 <- gf[["psi2"]]
    if (is.matrix(psi2) && !any(is.na(psi2))) {
      r_   <- dat[["r"]]; y_ <- dat[["y"]]
      w_hat <- res[["w.hat"]]
      ens_ps <- as.vector(res[["ps.matrix"]] %*% w_hat)
      pi_j   <- ebmr[["ps_fit.list"]][[1]][["fitted.values"]]
      des_j  <- ebmr[["ps_fit.list"]][[1]][["design_matrix"]]
      lt_j   <- ebmr[["ps_fit.list"]][[1]][["link_type"]]
      if (lt_j == "logistic") {
        dot_pi <- des_j * (pi_j * (1 - pi_j))
      } else {
        dot_pi <- -des_j * (pi_j * (1 - pi_j))
      }
      ry_ps2 <- as.vector(r_ * y_ * (ens_ps^(-2)))
      H_alpha <- colMeans(dot_pi * w_hat * ry_ps2)
      iid2 <- as.vector(r_ / ens_ps * y_) - as.vector(t(H_alpha) %*% psi2)
      se2_vals[i] <- sqrt(var(iid2) / n_val)
    }
  }, error = function(e) NULL)
}

valid <- !is.na(mu_vals)
conv <- conv_vals[valid]
mu_v <- mu_vals[valid]
se1_v <- se1_vals[valid]
se2_v <- se2_vals[valid]

cat(sprintf("mu.true = %.4f\n", mu_true))
cat(sprintf("Total valid: %d, converged: %d, not converged: %d\n", sum(valid), sum(conv), sum(!conv)))
cat(sprintf("se2 available: %d\n\n", sum(!is.na(se2_v))))

# Overall
esd <- sd(mu_v)
ese1 <- mean(se1_v)
ese2 <- mean(se2_v, na.rm=TRUE)
ci1_l <- mu_v - 1.96*se1_v; ci1_u <- mu_v + 1.96*se1_v
cp1 <- mean(mu_true >= ci1_l & mu_true <= ci1_u)
ci2_l <- mu_v - 1.96*se2_v; ci2_u <- mu_v + 1.96*se2_v
cp2 <- mean(mu_true >= ci2_l & mu_true <= ci2_u, na.rm=TRUE)
cat(sprintf("=== Overall ===\n"))
cat(sprintf("Bias=%.4f, ESD=%.4f\n", mean(mu_v)-mu_true, esd))
cat(sprintf("ESE1=%.4f, ESE1/ESD=%.3f, CP1=%.3f\n", ese1, ese1/esd, cp1))
cat(sprintf("ESE2=%.4f, ESE2/ESD=%.3f, CP2=%.3f\n\n", ese2, ese2/esd, cp2))

# Split by convergence
for (label in c("Converged", "NOT converged")) {
  sel <- if (label == "Converged") conv else !conv
  n_sel <- sum(sel)
  mu_s <- mu_vals[valid][sel]
  se1_s <- se1_vals[valid][sel]
  se2_s <- se2_vals[valid][sel]
  esd_s <- sd(mu_s)
  ese1_s <- mean(se1_s)
  ese2_s <- mean(se2_s, na.rm=TRUE)
  ci1_l <- mu_s - 1.96*se1_s; ci1_u <- mu_s + 1.96*se1_s
  cp1_s <- mean(mu_true >= ci1_l & mu_true <= ci1_u)
  ci2_l <- mu_s - 1.96*se2_s; ci2_u <- mu_s + 1.96*se2_s
  cp2_s <- mean(mu_true >= ci2_l & mu_true <= ci2_u, na.rm=TRUE)
  cat(sprintf("=== %s reps (n=%d, se2 avail=%d) ===\n", label, n_sel, sum(!is.na(se2_s))))
  cat(sprintf("Bias=%.4f, ESD=%.4f\n", mean(mu_s)-mu_true, esd_s))
  cat(sprintf("ESE1=%.4f, ESE1/ESD=%.3f, CP1=%.3f\n", ese1_s, ese1_s/esd_s, cp1_s))
  cat(sprintf("ESE2=%.4f, ESE2/ESD=%.3f, CP2=%.3f\n\n", ese2_s, ese2_s/esd_s, cp2_s))
}

# Alpha and convergence diagnostics
cat(sprintf("=== alpha[2] (y coef) ===\n"))
cat(sprintf("Converged:     mean=%.4f, sd=%.4f\n", mean(alpha2_vals[valid][conv]), sd(alpha2_vals[valid][conv])))
cat(sprintf("Not converged: mean=%.4f, sd=%.4f\n\n", mean(alpha2_vals[valid][!conv]), sd(alpha2_vals[valid][!conv])))

cat(sprintf("=== grad norm ===\n"))
cat(sprintf("Converged:     mean=%.2e, max=%.2e\n", mean(grad_vals[valid][conv]), max(grad_vals[valid][conv])))
cat(sprintf("Not converged: mean=%.2e, max=%.2e\n", mean(grad_vals[valid][!conv]), max(grad_vals[valid][!conv])))
