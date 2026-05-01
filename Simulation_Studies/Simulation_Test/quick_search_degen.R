setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

n_val <- 2000
ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

model_idx <- 3

# Get a normal rep's solution as init
dat1 <- all_data[1:n_val, ]
single_ps <- list(
  formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
  h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]],
  alpha_init.list = list(NULL),
  optimizer = "constrained_nr"
)
ebmr1 <- EBMRAlgorithmFast4$new("y", single_ps, dat1, W_fn)
normal_alpha <- ebmr1$ps_fit.list[[1]]$gmm_fit$estimates
cat(sprintf("Normal rep 1 alpha: %s, cond=%.1e\n",
    paste(round(normal_alpha, 3), collapse=", "),
    ebmr1$ps_fit.list[[1]]$gmm_fit$opt$solution_cond))

# Test 3 degenerate reps: init from normal_alpha
for (rep_i in c(35, 68, 166)) {
  cat(sprintf("\n=== Rep %d ===\n", rep_i))
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]

  # Run with normal_alpha as init
  single_ps$alpha_init.list <- list(normal_alpha)
  single_ps$optimizer <- "constrained_nr"
  ebmr_cnr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
  g_cnr <- ebmr_cnr$ps_fit.list[[1]]$gmm_fit
  cat(sprintf("  CNR from normal init: grad=%.2e, cond=%.2e, alpha=%s\n",
      g_cnr$opt$final_grad_norm, g_cnr$opt$solution_cond,
      paste(round(g_cnr$estimates, 3), collapse=", ")))

  single_ps$optimizer <- "L-BFGS-B"
  single_ps$alpha_init.list <- list(normal_alpha)
  ebmr_lb <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
  g_lb <- ebmr_lb$ps_fit.list[[1]]$gmm_fit
  cat(sprintf("  LB  from normal init: grad=%.2e, cond=%.2e, alpha=%s\n",
      g_lb$opt$final_grad_norm, g_lb$opt$solution_cond,
      paste(round(g_lb$estimates, 3), collapse=", ")))

  # Also try from zero
  single_ps$alpha_init.list <- list(NULL)
  single_ps$optimizer <- "L-BFGS-B"
  ebmr_lb0 <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
  g_lb0 <- ebmr_lb0$ps_fit.list[[1]]$gmm_fit
  cat(sprintf("  LB  from zero:        grad=%.2e, cond=%.2e, alpha=%s\n",
      g_lb0$opt$final_grad_norm, g_lb0$opt$solution_cond,
      paste(round(g_lb0$estimates, 3), collapse=", ")))

  single_ps$alpha_init.list <- list(NULL)
}
