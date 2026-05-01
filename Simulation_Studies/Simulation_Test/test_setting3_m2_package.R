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
data_file <- misspecified_model_all_data_file.list[["setting3"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
mu_true <- get_mu_true("setting3")

# Model 2 only
ps_sub2 <- list(
  formula.list = ps_spec[["formula.list"]][2],
  h_alpha.list = ps_spec[["h_alpha.list"]][2],
  inv_link     = ps_spec[["inv_link"]],
  outcome      = ps_spec[["outcome"]]
)

cat(sprintf("Formula: %s\n", deparse(ps_spec[["formula.list"]][[2]])))
cat(sprintf("h_alpha: %s\n", paste(ps_spec[["h_alpha.list"]][[2]], collapse=", ")))
cat(sprintf("mu.true = %.4f\n\n", mu_true))

mu_vals   <- rep(NA_real_, n_reps)
se1_vals  <- rep(NA_real_, n_reps)
se2_vals  <- rep(NA_real_, n_reps)
conv_vals <- rep(NA, n_reps)
grad_vals <- rep(NA_real_, n_reps)
alpha2_vals <- rep(NA_real_, n_reps)

for (i in 1:n_reps) {
  dat <- all_data[((i-1)*n_val + 1):(i*n_val), ]
  tryCatch({
    ebmr <- EBMRAlgorithmFast4[["new"]]("y", ps_sub2, dat, W_func)
    res  <- ebmr[["EBMR_IPW"]](h_nu=h_nu_fn, type="HT", se.fit=TRUE)
    mu_vals[i]  <- res[["mu_ipw"]]
    se1_vals[i] <- res[["se_ipw"]]
    gf <- ebmr[["ps_fit.list"]][[1]][["gmm_fit"]]
    conv_vals[i] <- gf[["opt"]][["converged"]]
    grad_vals[i] <- gf[["opt"]][["final_grad_norm"]]
    alpha2_vals[i] <- gf[["estimates"]][2]

    # se2 via package psi2
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

cat(sprintf("Total valid: %d, converged: %d, not converged: %d\n", sum(valid), sum(conv), sum(!conv)))
cat(sprintf("se2 available: %d\n\n", sum(!is.na(se2_v))))

esd <- sd(mu_v)
ese1 <- mean(se1_v)
ese2 <- mean(se2_v, na.rm=TRUE)
cat(sprintf("=== Overall ===\n"))
cat(sprintf("Bias=%.4f, ESD=%.4f\n", mean(mu_v)-mu_true, esd))
cat(sprintf("ESE1=%.4f, ESE1/ESD=%.3f\n", ese1, ese1/esd))
cat(sprintf("ESE2=%.4f, ESE2/ESD=%.3f\n\n", ese2, ese2/esd))

for (label in c("Converged", "NOT converged")) {
  sel <- if (label == "Converged") conv else !conv
  n_sel <- sum(sel)
  if (n_sel == 0) next
  mu_s <- mu_vals[valid][sel]
  se1_s <- se1_vals[valid][sel]
  se2_s <- se2_vals[valid][sel]
  esd_s <- sd(mu_s)
  ese1_s <- mean(se1_s)
  ese2_s <- mean(se2_s, na.rm=TRUE)
  cat(sprintf("=== %s reps (n=%d) ===\n", label, n_sel))
  cat(sprintf("Bias=%.4f, ESD=%.4f\n", mean(mu_s)-mu_true, esd_s))
  cat(sprintf("ESE1=%.4f, ESE1/ESD=%.3f\n", ese1_s, ese1_s/esd_s))
  cat(sprintf("ESE2=%.4f, ESE2/ESD=%.3f\n\n", ese2_s, ese2_s/esd_s))
}

cat(sprintf("=== alpha[2] ===\n"))
cat(sprintf("mean=%.4f, sd=%.4f, range=[%.4f, %.4f]\n",
    mean(alpha2_vals[valid]), sd(alpha2_vals[valid]),
    min(alpha2_vals[valid]), max(alpha2_vals[valid])))
