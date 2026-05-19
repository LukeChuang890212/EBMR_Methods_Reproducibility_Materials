setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EPS", quiet = TRUE)

n_val <- 2000
n_reps <- 200
ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

model_idx <- 3

# Step 1: Identify the degenerate reps
degen_reps <- c()
good_reps <- c()
for (rep_i in 1:n_reps) {
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL),
    optimizer = "constrained_nr"
  )
  tryCatch({
    ebmr <- EPS$new("y", single_ps, dat, W_fn)
    sc <- ebmr$ps_fit.list[[1]]$gmm_fit$opt$solution_cond
    if (!is.na(sc) && sc >= 1e8) {
      degen_reps <- c(degen_reps, rep_i)
    } else {
      good_reps <- c(good_reps, rep_i)
    }
  }, error = function(e) { degen_reps <<- c(degen_reps, rep_i) })
}
cat(sprintf("Degenerate reps (%d): %s\n", length(degen_reps), paste(degen_reps, collapse=", ")))

# Step 2: For each degenerate rep, try many random starts with constrained_nr
# Also compare with L-BFGS-B and check the objective landscape
cat("\n=== Diagnosing degenerate reps ===\n")
for (rep_i in degen_reps[1:min(5, length(degen_reps))]) {
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
  r_vec <- dat[["r"]]

  # Default constrained_nr result
  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL),
    optimizer = "constrained_nr"
  )
  ebmr_cnr <- EPS$new("y", single_ps, dat, W_fn)
  gmm_cnr <- ebmr_cnr$ps_fit.list[[1]]$gmm_fit

  # L-BFGS-B result
  single_ps$optimizer <- "L-BFGS-B"
  ebmr_lb <- EPS$new("y", single_ps, dat, W_fn)
  gmm_lb <- ebmr_lb$ps_fit.list[[1]]$gmm_fit

  cat(sprintf("\nRep %d:\n", rep_i))
  cat(sprintf("  CNR:     obj=%.6f, grad=%.2e, cond=%.2e, alpha=%s\n",
      gmm_cnr$opt$objective, gmm_cnr$opt$final_grad_norm, gmm_cnr$opt$solution_cond,
      paste(round(gmm_cnr$estimates, 4), collapse=", ")))
  cat(sprintf("  L-BFGS-B: obj=%.6f, grad=%.2e, cond=%.2e, alpha=%s\n",
      gmm_lb$opt$objective, gmm_lb$opt$final_grad_norm, gmm_lb$opt$solution_cond,
      paste(round(gmm_lb$estimates, 4), collapse=", ")))

  # Check pi_hat range for both
  pi_cnr <- ebmr_cnr$ps_fit.list[[1]]$fitted.values
  pi_lb <- ebmr_lb$ps_fit.list[[1]]$fitted.values
  cat(sprintf("  CNR pi:  [%.6f, %.6f], mean=%.4f\n", min(pi_cnr), max(pi_cnr), mean(pi_cnr)))
  cat(sprintf("  LB  pi:  [%.6f, %.6f], mean=%.4f\n", min(pi_lb), max(pi_lb), mean(pi_lb)))
  cat(sprintf("  r mean=%.4f\n", mean(r_vec)))

  # Try many random starts with constrained_nr (no fallback)
  # We want to see if ANY non-degenerate solution exists
  best_nondegen_obj <- Inf
  best_nondegen_grad <- Inf
  best_nondegen_alpha <- NULL
  n_tries <- 20
  for (trial in 1:n_tries) {
    single_ps$optimizer <- "constrained_nr"
    single_ps$alpha_init.list <- list(rnorm(length(gmm_cnr$estimates)) * 0.5)
    tryCatch({
      ebmr_t <- EPS$new("y", single_ps, dat, W_fn)
      gmm_t <- ebmr_t$ps_fit.list[[1]]$gmm_fit
      sc <- gmm_t$opt$solution_cond
      if (!is.na(sc) && sc < 1e8 && gmm_t$opt$final_grad_norm < best_nondegen_grad) {
        best_nondegen_obj <- gmm_t$opt$objective
        best_nondegen_grad <- gmm_t$opt$final_grad_norm
        best_nondegen_alpha <- gmm_t$estimates
      }
    }, error = function(e) NULL)
  }
  single_ps$alpha_init.list <- list(NULL)  # reset

  if (!is.null(best_nondegen_alpha)) {
    cat(sprintf("  Best non-degen (20 tries): obj=%.6f, grad=%.2e, alpha=%s\n",
        best_nondegen_obj, best_nondegen_grad,
        paste(round(best_nondegen_alpha, 4), collapse=", ")))
  } else {
    cat("  No non-degenerate solution found in 20 random starts\n")
  }
}

# Step 3: Compare a good rep to understand the difference
cat("\n=== Good rep comparison ===\n")
rep_good <- good_reps[1]
dat <- all_data[((rep_good-1)*n_val + 1):(rep_good*n_val), ]
single_ps <- list(
  formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
  h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]],
  alpha_init.list = list(NULL),
  optimizer = "constrained_nr"
)
ebmr_g <- EPS$new("y", single_ps, dat, W_fn)
gmm_g <- ebmr_g$ps_fit.list[[1]]$gmm_fit
pi_g <- ebmr_g$ps_fit.list[[1]]$fitted.values
cat(sprintf("Good rep %d: obj=%.6f, grad=%.2e, cond=%.2e\n",
    rep_good, gmm_g$opt$objective, gmm_g$opt$final_grad_norm, gmm_g$opt$solution_cond))
cat(sprintf("  alpha=%s\n", paste(round(gmm_g$estimates, 4), collapse=", ")))
cat(sprintf("  pi: [%.6f, %.6f], mean=%.4f\n", min(pi_g), max(pi_g), mean(pi_g)))
