setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
library(EBMRalgorithmFast4)

data_file <- misspecified_model_all_data_file.list$setting3$miss50[[1]]
all_data <- readRDS(data_file)
n_val <- 2000

ps_spec <- get_ps_spec("9-alt1")
subset_ps_spec <- list(
  formula.list = ps_spec$formula.list[2],
  h_alpha.list = ps_spec$h_alpha.list[2],
  inv_link = ps_spec$inv_link,
  outcome = ps_spec$outcome
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

# Identify outlier reps from the fresh run
# Use the worst outliers from the old results as reference
outlier_reps <- c(929, 438, 430, 221, 180, 781, 754, 615, 711, 587)
normal_reps <- c(1, 2, 3, 4, 5)

cat(sprintf("%-6s %-8s %-6s %-10s %-12s %-10s %-10s %-10s %-10s %-40s\n",
            "Rep", "Type", "Iter", "GradNorm", "Objective", "Converged",
            "mu_ipw", "min(ps)", "max(r/ps)", "Alpha"))

for (rep_i in c(outlier_reps, normal_reps)) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]
  label <- if (rep_i %in% outlier_reps) "OUT" else "NORM"

  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", subset_ps_spec, dat, W_func)
    fit <- ebmr$ps_fit.list[[1]]
    gmm <- fit$gmm_fit
    ps <- fit$fitted.values
    r_vec <- as.numeric(dat$r)

    res <- ebmr$EBMR_IPW(
      h_nu = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2),
      se.fit = TRUE, type = "HT"
    )

    cat(sprintf("%-6d %-8s %-6d %-10.2e %-12.6e %-10s %-10.4f %-10.6f %-10.1f [%s]\n",
                rep_i, label,
                gmm$opt$iterations,
                gmm$opt$final_grad_norm,
                gmm$opt$objective,
                gmm$opt$converged,
                res$mu_ipw,
                min(ps), max(r_vec/ps),
                paste(round(fit$coefficients, 4), collapse=", ")))

    # Check if this is a degenerate solution: try glm as reference
    glm_fit <- glm(r ~ dat$u1 + dat$u2 + dat$z2 + dat$y,
                   family = binomial(link = "logit"), data = dat)
    glm_coef <- coef(glm_fit)
    # glm uses r ~ u1 + u2 + z2 + y but our model is r ~ y + u1 + z2
    # Need to match the covariate order
    cat(sprintf("       GLM coef: [%s]\n", paste(round(glm_coef, 4), collapse=", ")))
    cat(sprintf("       GMM coef: [%s]\n", paste(round(fit$coefficients, 4), collapse=", ")))

    # Compare objectives: GMM at its solution vs GMM at glm solution
    # Reconstruct to check
    ps_at_gmm <- ps
    glm_ps <- glm_fit$fitted.values

    # Moment conditions at GMM solution
    h_x <- model.matrix(~ dat$y + dat$u1 + dat$z2)
    g_gmm <- (r_vec / ps_at_gmm - 1) * h_x
    G_gmm <- colMeans(g_gmm)
    W_mat <- solve(crossprod(g_gmm) / n_val)
    obj_gmm <- as.numeric(t(G_gmm) %*% W_mat %*% G_gmm)

    # Moment conditions at GLM solution (need to reorder)
    g_glm <- (r_vec / glm_ps - 1) * h_x
    G_glm <- colMeans(g_glm)
    obj_glm <- as.numeric(t(G_glm) %*% W_mat %*% G_glm)

    cat(sprintf("       Obj at GMM sol: %.6e, Obj at GLM sol: %.6e\n", obj_gmm, obj_glm))

  }, error = function(e) {
    cat(sprintf("%-6d %-8s ERROR: %s\n", rep_i, label, conditionMessage(e)))
  })
  cat("\n")
}

cat("Done!\n")
