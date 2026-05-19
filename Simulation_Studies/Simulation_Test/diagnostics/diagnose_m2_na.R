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

for (opt in c("L-BFGS-B", "constrained_nr")) {
  cat(sprintf("\n=== M2 %s: diagnosing NAs ===\n", opt))
  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[2]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[2]]),
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL),
    optimizer = opt
  )

  for (rep_i in 1:n_reps) {
    dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
    tryCatch({
      ebmr <- EPS$new("y", single_ps, dat, W_fn)
      gmm_fit <- ebmr$ps_fit.list[[1]]$gmm_fit
      psi2 <- gmm_fit$psi2
      mu_ok <- !is.na(mean(dat[["r"]] / ebmr$ps_fit.list[[1]]$fitted.values * dat[["y"]]))
      se2_ok <- !is.null(psi2) && is.matrix(psi2) && !any(is.na(psi2))
      if (!mu_ok || !se2_ok) {
        cat(sprintf("  Rep %d: mu_ok=%s, se2_ok=%s, converged=%s, grad=%.2e, cond=%.2e\n",
            rep_i, mu_ok, se2_ok, gmm_fit$opt$converged,
            gmm_fit$opt$final_grad_norm, gmm_fit$opt$solution_cond))
        if (!se2_ok && !is.null(psi2)) {
          cat(sprintf("    psi2 class=%s, is.matrix=%s, any_na=%s\n",
              class(psi2), is.matrix(psi2), any(is.na(psi2))))
        }
        if (!se2_ok && is.null(psi2)) cat("    psi2 is NULL\n")
        if (!se2_ok && identical(psi2, NA)) cat("    psi2 is NA\n")
      }
    }, error = function(e) {
      cat(sprintf("  Rep %d: ERROR: %s\n", rep_i, e$message))
    })
  }
}
