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

# Test on a few outlier reps and normal reps
outlier_reps <- c(929, 438, 430, 221, 180)
normal_reps <- c(1, 2, 3, 4, 5)

cat(sprintf("%-6s %-8s %-10s %-10s %-12s %-10s\n",
            "Rep", "Type", "mu_ipw", "se_ipw", "Objective", "Iter"))

for (rep_i in c(outlier_reps, normal_reps)) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]
  label <- if (rep_i %in% outlier_reps) "OUT" else "NORM"

  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", subset_ps_spec, dat, W_func)
    fit <- ebmr$ps_fit.list[[1]]
    res <- ebmr$EBMR_IPW(
      h_nu = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2),
      se.fit = TRUE, type = "HT"
    )
    cat(sprintf("%-6d %-8s %-10.4f %-10.4f %-12.6e %-10d [%s]\n",
                rep_i, label, res$mu_ipw, res$se_ipw,
                fit$gmm_fit$opt$objective, fit$gmm_fit$opt$iterations,
                paste(round(fit$coefficients, 4), collapse = ", ")))
  }, error = function(e) {
    cat(sprintf("%-6d %-8s ERROR: %s\n", rep_i, label, conditionMessage(e)))
  })
}

cat("\nDone!\n")
