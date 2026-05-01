#------------------------------------------------------------------------------#
# Detailed SE Investigation - Why is analytical SE wrong?
#------------------------------------------------------------------------------#

library(EBMRalgorithmFast)
library(tidyverse)
library(numDeriv)

source("MHD_functions.R")

#------------------------------------------------------------------------------#
# Setup
#------------------------------------------------------------------------------#
original_data <- read.csv("data_application.csv")
percent <- original_data$Percentage
n <- 2486
class_n <- round(n * percent / 100)
dat <- gen_data(original_data, class_n, n)

cat("================================================================================\n")
cat("              DETAILED SE INVESTIGATION                                         \n")
cat("================================================================================\n\n")

#------------------------------------------------------------------------------#
# Test both health groups and both MNAR models
#------------------------------------------------------------------------------#

full_ps_specifications <- list(
  formula.list = list(
    r ~ o(teacher_report) + father,
    r ~ o(teacher_report) + parent_report
  ),
  h_x_names.list = list(
    c("father", "parent_report"),
    c("father", "parent_report")
  ),
  inv_link = function(eta) 1 / (1 + exp(eta))
)

W <- function(g.matrix) {
  solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
}

init.list <- list(
  list(c(-0.3834341, 1.07945126, -0.1318883),
       c(-4.2670736, 5.9582123, -0.5972129)),
  list(c(0.09019492, -1.511548, -0.3946616),
       c(-1.969056, 3.5823611, -1.37232734))
)

model_names <- c("Model 1 (o(y) + father)", "Model 2 (o(y) + parent_report)")

set.seed(12345)
n_boot <- 500

results_summary <- data.frame()

for (health_val in 0:1) {
  subdat <- dat[dat$health == health_val, ]
  n_sub <- nrow(subdat)

  cat("================================================================================\n")
  cat("Health =", health_val, ", n =", n_sub, "\n")
  cat("================================================================================\n\n")

  for (model_idx in 1:2) {
    cat(model_names[model_idx], "\n")
    cat(paste(rep("-", 60), collapse = ""), "\n")

    # Fit model
    ps_specifications <- list(
      formula.list = full_ps_specifications$formula.list[model_idx],
      h_x_names.list = full_ps_specifications$h_x_names.list[model_idx],
      alpha_init.list = list(init.list[[health_val + 1]][[model_idx]]),
      inv_link = full_ps_specifications$inv_link
    )

    ebmr <- EBMRAlgorithmFast$new("teacher_report", ps_specifications, subdat, W)
    ps_fit <- ebmr$ps_fit.list[[1]]

    alpha_hat <- ps_fit$coefficients
    se_analytical <- ps_fit$se

    cat("Point estimates: ", round(alpha_hat, 4), "\n")
    cat("Analytical SE:   ", round(se_analytical, 4), "\n")

    # Bootstrap
    r_vec <- subdat$r
    y_vec <- subdat$teacher_report
    if (model_idx == 1) {
      x_vec <- subdat$father
    } else {
      x_vec <- subdat$parent_report
    }
    h_x <- cbind(1, subdat$father, subdat$parent_report)
    design_matrix <- cbind(1, y_vec, x_vec)
    inv_link <- function(eta) 1 / (1 + exp(eta))

    boot_estimates <- matrix(NA, n_boot, 3)

    for (b in 1:n_boot) {
      idx <- sample(1:n_sub, replace = TRUE)
      boot_dat <- subdat[idx, ]

      r_boot <- boot_dat$r
      y_boot <- boot_dat$teacher_report
      if (model_idx == 1) {
        x_boot <- boot_dat$father
      } else {
        x_boot <- boot_dat$parent_report
      }
      h_x_boot <- cbind(1, boot_dat$father, boot_dat$parent_report)
      design_boot <- cbind(1, y_boot, x_boot)

      obj_boot <- function(param) {
        eta <- design_boot %*% param
        pi_vec <- inv_link(eta)
        rw <- r_boot / as.vector(pi_vec)
        g.mat <- (rw - 1) * h_x_boot
        G <- colMeans(g.mat)
        W_mat <- tryCatch(solve(crossprod(g.mat) / n_sub), error = function(e) diag(3))
        return(as.numeric(t(G) %*% W_mat %*% G))
      }

      tryCatch({
        opt <- optim(alpha_hat, obj_boot, method = "L-BFGS-B",
                     lower = rep(-10, 3), upper = rep(10, 3),
                     control = list(maxit = 1000))
        if (opt$convergence == 0) {
          boot_estimates[b, ] <- opt$par
        }
      }, error = function(e) {})
    }

    boot_valid <- boot_estimates[complete.cases(boot_estimates), ]
    se_boot <- apply(boot_valid, 2, sd)

    cat("Bootstrap SE:    ", round(se_boot, 4), "(n_valid =", nrow(boot_valid), ")\n")
    cat("Ratio (Ana/Boot):", round(se_analytical / se_boot, 3), "\n")

    # z-values and p-values
    z_analytical <- alpha_hat / se_analytical
    z_boot <- alpha_hat / se_boot
    p_analytical <- 2 * pnorm(-abs(z_analytical))
    p_boot <- 2 * pnorm(-abs(z_boot))

    cat("\nCoefficient on y (alpha_1):\n")
    cat(sprintf("  Analytical: z = %.3f, p = %.4f\n", z_analytical[2], p_analytical[2]))
    cat(sprintf("  Bootstrap:  z = %.3f, p = %.4f\n", z_boot[2], p_boot[2]))

    results_summary <- rbind(results_summary, data.frame(
      health = health_val,
      model = model_idx,
      alpha_y = alpha_hat[2],
      se_analytical = se_analytical[2],
      se_boot = se_boot[2],
      ratio = se_analytical[2] / se_boot[2],
      p_analytical = p_analytical[2],
      p_boot = p_boot[2]
    ))

    cat("\n\n")
  }
}

#------------------------------------------------------------------------------#
# Summary table
#------------------------------------------------------------------------------#
cat("================================================================================\n")
cat("SUMMARY: COEFFICIENT ON Y ACROSS ALL MODELS\n")
cat("================================================================================\n\n")

print(results_summary, row.names = FALSE)

cat("\n\nKey finding: Analytical SE underestimates the true SE by a factor of",
    round(mean(results_summary$ratio), 2), "on average.\n")

cat("\nUsing bootstrap SE, the p-values for the coefficient on y are:\n")
for (i in 1:nrow(results_summary)) {
  with(results_summary[i, ], {
    cat(sprintf("  Health=%d, Model=%d: alpha_y = %.3f, p = %.4f %s\n",
                health, model, alpha_y, p_boot,
                ifelse(p_boot < 0.05, "*", ifelse(p_boot < 0.1, ".", ""))))
  })
}

cat("\n================================================================================\n")
