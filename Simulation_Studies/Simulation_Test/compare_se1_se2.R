setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
  library(EBMRalgorithmFast4)
})

W_func  <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_nu_fn <- function(dat) cbind(u1=dat[["u1"]], u2=dat[["u2"]], z1=dat[["z1"]], z2=dat[["z2"]])
n_val   <- 2000
n_reps  <- 200

# Setting3, scenario 9-alt1, model 2 only
ps_spec <- get_ps_spec("9-alt1")
single_ps_spec <- list(
  formula.list = ps_spec[["formula.list"]][2],
  h_alpha.list = ps_spec[["h_alpha.list"]][2],
  inv_link     = ps_spec[["inv_link"]],
  outcome      = ps_spec[["outcome"]]
)
mu_true <- get_mu_true("setting3")
all_data <- readRDS(misspecified_model_all_data_file.list[["setting3"]][["miss50"]][[1]])

cat(sprintf("mu.true = %.4f\n", mu_true))
cat(sprintf("=== Comparing se1 vs se2 (%d reps, model 2 alone, setting3) ===\n\n", n_reps))

mu_vals  <- numeric(n_reps)
se1_vals <- numeric(n_reps)
se2_vals <- numeric(n_reps)

n_err <- 0
for (i in 1:n_reps) {
  dat <- all_data[((i-1)*n_val + 1):(i*n_val), ]
  tryCatch({
    ebmr <- EBMRAlgorithmFast4[["new"]]("y", single_ps_spec, dat, W_func)
    res  <- ebmr[["EBMR_IPW"]](h_nu=h_nu_fn, type="HT", se.fit=TRUE)
    mu_vals[i]  <- res[["mu_ipw"]]
    se1_vals[i] <- res[["se_ipw"]]  # uses psi1 (original)

    # Compute se_ipw using psi2 instead of psi1
    gmm_fit  <- ebmr[["ps_fit.list"]][[1]][["gmm_fit"]]
    psi2_mat <- gmm_fit[["psi2"]]  # alpha_dim x n
    if (!is.null(psi2_mat) && is.matrix(psi2_mat) && !any(is.na(psi2_mat))) {
      r_   <- dat[["r"]]
      y_   <- dat[["y"]]
      pi_  <- ebmr[["ps_fit.list"]][[1]][["fitted.values"]]
      des_ <- ebmr[["ps_fit.list"]][[1]][["design_matrix"]]

      # Determine link type
      lt <- ebmr[["ps_fit.list"]][[1]][["link_type"]]
      if (lt == "logistic") {
        dot_pi_ <- des_ * (pi_ * (1 - pi_))
      } else if (lt == "logistic_complement") {
        dot_pi_ <- -des_ * (pi_ * (1 - pi_))
      } else {
        stop("Unsupported link type: ", lt)
      }

      ry_ps2_  <- as.vector(r_ * y_ * (pi_^(-2)))
      H_alpha_ <- colMeans(dot_pi_ * ry_ps2_)
      iid_     <- as.vector(r_ / pi_ * y_) - as.vector(t(H_alpha_) %*% psi2_mat)
      se2_vals[i] <- sqrt(var(iid_) / n_val)
    } else {
      se2_vals[i] <- NA
    }
  }, error = function(e) {
    n_err <<- n_err + 1
    mu_vals[i]  <<- NA
    se1_vals[i] <<- NA
    se2_vals[i] <<- NA
    if (n_err <= 5) cat(sprintf("  Error rep %d: %s\n", i, conditionMessage(e)))
  })
  if (i %% 50 == 0) cat(sprintf("  %d / %d done (errors: %d)\n", i, n_reps, n_err))
}

cat(sprintf("\nErrors: %d / %d\n", n_err, n_reps))

# Results
valid <- !is.na(mu_vals)
mu_v  <- mu_vals[valid]
se1_v <- se1_vals[valid]
se2_v <- se2_vals[valid]

esd  <- sd(mu_v)
ese1 <- mean(se1_v)
ese2 <- mean(se2_v, na.rm = TRUE)

cat(sprintf("\n=== RESULTS ===\n"))
cat(sprintf("Valid reps: %d\n", sum(valid)))
cat(sprintf("Bias: %.4f\n", mean(mu_v) - mu_true))
cat(sprintf("ESD:  %.4f\n", esd))
cat(sprintf("ESE1 (original):  %.4f  (ESE1/ESD = %.3f)\n", ese1, ese1/esd))
cat(sprintf("ESE2 (corrected): %.4f  (ESE2/ESD = %.3f)\n", ese2, ese2/esd))

# Paired comparison
ratio <- se2_v / se1_v
cat(sprintf("\nSE2/SE1 ratio: Mean=%.3f, Median=%.3f, Min=%.3f, Max=%.3f\n",
    mean(ratio, na.rm=TRUE), median(ratio, na.rm=TRUE),
    min(ratio, na.rm=TRUE), max(ratio, na.rm=TRUE)))

# Coverage probability
ci1_l <- mu_v - 1.96 * se1_v
ci1_u <- mu_v + 1.96 * se1_v
cp1   <- mean(mu_true >= ci1_l & mu_true <= ci1_u)

ci2_l <- mu_v - 1.96 * se2_v
ci2_u <- mu_v + 1.96 * se2_v
cp2   <- mean(mu_true >= ci2_l & mu_true <= ci2_u, na.rm = TRUE)

cat(sprintf("\nCP (se1): %.3f\n", cp1))
cat(sprintf("CP (se2): %.3f\n", cp2))
