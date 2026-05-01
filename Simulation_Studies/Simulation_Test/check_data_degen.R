setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

n_val <- 2000
n_reps <- 200
ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
model_idx <- 3

# First identify degenerate reps
cond_vals <- rep(NA_real_, n_reps)
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
    ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
    cond_vals[rep_i] <- ebmr$ps_fit.list[[1]]$gmm_fit$opt$solution_cond
  }, error = function(e) NULL)
}

degen <- which(!is.na(cond_vals) & cond_vals >= 1e8)
normal <- which(!is.na(cond_vals) & cond_vals < 1e8)
cat(sprintf("Degenerate reps: %d, Normal reps: %d\n", length(degen), length(normal)))

# Compare data characteristics
cat("\n=== Data characteristics: Degenerate vs Normal ===\n")
cat(sprintf("  %20s  %12s  %12s  %10s\n", "Statistic", "Degen(mean)", "Normal(mean)", "p-value"))

collect <- function(stat_fn, label) {
  vals_d <- sapply(degen, function(i) {
    dat <- all_data[((i-1)*n_val + 1):(i*n_val), ]
    stat_fn(dat)
  })
  vals_n <- sapply(normal, function(i) {
    dat <- all_data[((i-1)*n_val + 1):(i*n_val), ]
    stat_fn(dat)
  })
  tt <- tryCatch(t.test(vals_d, vals_n), error = function(e) list(p.value = NA))
  cat(sprintf("  %20s  %12.4f  %12.4f  %10.4f\n",
      label, mean(vals_d), mean(vals_n), tt$p.value))
  invisible(list(degen = vals_d, normal = vals_n))
}

# Check data columns
cat(sprintf("\nColumns in data: %s\n", paste(names(all_data[1:n_val, ]), collapse=", ")))

collect(function(d) mean(d[["r"]]), "mean(r)")
collect(function(d) mean(d[["y"]]), "mean(y)")
collect(function(d) mean(d[["y"]][d[["r"]]==1]), "mean(y|r=1)")
collect(function(d) sd(d[["y"]]), "sd(y)")
collect(function(d) sd(d[["y"]][d[["r"]]==1]), "sd(y|r=1)")

# Check covariates
cov_names <- setdiff(names(all_data[1:n_val, ]), c("r", "y"))
for (cn in cov_names[1:min(length(cov_names), 8)]) {
  collect(function(d) mean(d[[cn]]), sprintf("mean(%s)", cn))
  collect(function(d) sd(d[[cn]]), sprintf("sd(%s)", cn))
}

# Check covariate-missingness relationship
cat("\n=== Covariate-missingness correlation ===\n")
for (cn in cov_names[1:min(length(cov_names), 5)]) {
  collect(function(d) cor(d[[cn]], d[["r"]]), sprintf("cor(%s,r)", cn))
}

# Check the PS model design matrix columns
cat("\n=== PS model design matrix properties ===\n")
dat1 <- all_data[1:n_val, ]
single_ps <- list(
  formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
  h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]],
  alpha_init.list = list(NULL),
  optimizer = "L-BFGS-B"
)
ebmr1 <- EBMRAlgorithmFast4$new("y", single_ps, dat1, W_fn)
dm_cols <- colnames(ebmr1$ps_fit.list[[1]]$design_matrix)
hx_cols <- colnames(ebmr1$ps_fit.list[[1]]$h_x)
cat(sprintf("  Design matrix columns: %s\n", paste(dm_cols, collapse=", ")))
cat(sprintf("  h(x) columns: %s\n", paste(hx_cols, collapse=", ")))

# For each rep, compute design matrix stats
cat("\n=== Design matrix column means ===\n")
for (j in seq_along(dm_cols)) {
  collect(function(d) {
    ebmr_tmp <- EBMRAlgorithmFast4$new("y", single_ps, d, W_fn)
    mean(ebmr_tmp$ps_fit.list[[1]]$design_matrix[, j])
  }, sprintf("dm_col%d_mean", j))
}

# Check true propensity score (if available) in degenerate vs normal
cat("\n=== True missingness rate comparison ===\n")
r_vals <- collect(function(d) mean(d[["r"]]), "miss_rate")
cat(sprintf("\nDegenerate miss rates: %s\n",
    paste(round(r_vals$degen, 3), collapse=", ")))
cat(sprintf("Normal miss rate range: [%.3f, %.3f]\n",
    min(r_vals$normal), max(r_vals$normal)))
