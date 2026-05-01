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
dat <- all_data[1:n_val, ]
single_ps <- list(
  formula.list = list(ps_spec[["formula.list"]][[model_idx]]),
  h_alpha.list = list(ps_spec[["h_alpha.list"]][[model_idx]]),
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]],
  alpha_init.list = list(NULL),
  optimizer = "constrained_nr"
)
ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
gmm_fit <- ebmr$ps_fit.list[[1]]$gmm_fit
design_mat <- ebmr$ps_fit.list[[1]]$design_matrix
alpha_hat <- gmm_fit$estimates
pi_hat <- ebmr$ps_fit.list[[1]]$fitted.values
n <- n_val; r_vec <- dat[["r"]]

# Package g
g_pkg <- gmm_fit$g.matrix
cat(sprintf("Package g.matrix: dim = %d x %d\n", nrow(g_pkg), ncol(g_pkg)))

# Package h_x
h_x_pkg <- gmm_fit$h_x
cat(sprintf("Package h_x: dim = %d x %d\n", nrow(h_x_pkg), ncol(h_x_pkg)))

# My h_x
h_alpha_vars <- ps_spec[["h_alpha.list"]][[model_idx]]
cat(sprintf("h_alpha_vars: %s\n", paste(h_alpha_vars, collapse=", ")))
h_x_mine <- cbind(1, as.matrix(dat[, h_alpha_vars, drop = FALSE]))
cat(sprintf("My h_x: dim = %d x %d\n", nrow(h_x_mine), ncol(h_x_mine)))

# Compare h_x
cat(sprintf("h_x max diff: %.6e\n", max(abs(h_x_pkg - h_x_mine))))

# Compare pi_hat computation
eta <- as.vector(design_mat %*% alpha_hat)
pi_mine <- plogis(-eta)
cat(sprintf("pi max diff: %.6e\n", max(abs(pi_hat - pi_mine))))

# My g
g_mine <- (r_vec / pi_mine - 1) * h_x_mine
cat(sprintf("g max diff (mine vs pkg): %.6e\n", max(abs(g_mine - g_pkg))))

# Print first few rows of g
cat("\nFirst 5 rows of g_pkg:\n")
print(head(g_pkg, 5))
cat("\nFirst 5 rows of g_mine:\n")
print(head(g_mine, 5))

# Check column names
cat("\nPackage h_x colnames: ", paste(colnames(h_x_pkg), collapse=", "), "\n")
cat("My h_x colnames: ", paste(colnames(h_x_mine), collapse=", "), "\n")

# Check if h_x columns are in different order
if (ncol(h_x_pkg) == ncol(h_x_mine)) {
  for (i in 1:ncol(h_x_pkg)) {
    for (j in 1:ncol(h_x_mine)) {
      if (max(abs(h_x_pkg[, i] - h_x_mine[, j])) < 1e-10) {
        cat(sprintf("h_x_pkg col %d matches h_x_mine col %d\n", i, j))
      }
    }
  }
}
