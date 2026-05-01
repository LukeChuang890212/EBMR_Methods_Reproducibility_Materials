setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Data_Generation.r")
source("config/scenarios.R")
library(EBMRalgorithmFast4)
source("Basic_setup.r")

data_file <- misspecified_model_all_data_file.list$setting3$miss50[[1]]
all_data <- readRDS(data_file)
n_val <- 2000

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

# Compare model 1 (correct) vs model 2 (misspecified) on the SAME rep
ps_spec <- get_ps_spec("9-alt1")

for (model_idx in c(1, 2)) {
  subset_ps_spec <- list(
    formula.list = ps_spec$formula.list[model_idx],
    h_alpha.list = ps_spec$h_alpha.list[model_idx],
    inv_link = ps_spec$inv_link,
    outcome = ps_spec$outcome
  )

  cat(sprintf("\n===== MODEL %d (PS spec %d) =====\n", model_idx, model_idx))
  cat(sprintf("Formula: %s\n", deparse(ps_spec$formula.list[[model_idx]])))

  # Check 10 reps
  for (rep_i in c(1, 2, 3, 437, 929)) {
    dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]

    tryCatch({
      ebmr <- EBMRAlgorithmFast4$new("y", subset_ps_spec, dat, W_func)
      ps_fit <- ebmr$ps_fit.list[[1]]
      ps_fitted <- ps_fit$fitted.values
      r <- as.numeric(dat$r)
      y <- dat$y
      n <- length(y)

      # H_alpha_w concentration
      h_alpha_mat <- as.matrix(ps_fit$design_matrix)
      dot_pi <- h_alpha_mat * ps_fitted * (1 - ps_fitted)
      ry_ps_inv2 <- r * y * ps_fitted^(-2)

      # Contribution of top obs to H_alpha_w
      contrib_norm <- apply(dot_pi * ry_ps_inv2, 1, function(x) sqrt(sum(x^2))) / n
      total_norm <- sqrt(sum(colSums(dot_pi * ry_ps_inv2)^2)) / n
      top1_pct <- max(contrib_norm) / total_norm * 100

      # SE comparison
      res_ipw <- ebmr$EBMR_IPW(
        h_nu = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2),
        se.fit = TRUE, type = "HT"
      )
      base <- r / ps_fitted * y
      se_naive <- sqrt(var(base)/n)

      cat(sprintf("  Rep %3d: min(ps)=%.5f, max(r/ps)=%6.1f, top1_H%%=%.1f%%, se=%.4f, se_naive=%.4f, ratio=%.2f\n",
                  rep_i, min(ps_fitted), max(r/ps_fitted), top1_pct,
                  res_ipw$se_ipw, se_naive, res_ipw$se_ipw/se_naive))
    }, error = function(e) {
      cat(sprintf("  Rep %3d: ERROR: %s\n", rep_i, conditionMessage(e)))
    })
  }
}

cat("\nDone!\n")
