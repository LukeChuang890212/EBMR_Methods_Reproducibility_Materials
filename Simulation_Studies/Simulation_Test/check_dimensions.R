## Quick check: dimensions of design_mat and h_x for each model
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

n_val <- 2000
ps_spec <- get_ps_spec("9-alt1")
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

settings <- list(
  list(name = "S4 M3", setting = "setting4", m_idx = 3),
  list(name = "S4 M2", setting = "setting4", m_idx = 2),
  list(name = "S3 M2", setting = "setting3", m_idx = 2)
)

for (cfg in settings) {
  data_file <- misspecified_model_all_data_file.list[[cfg$setting]][["miss50"]][[1]]
  all_data <- readRDS(data_file)
  dat <- all_data[1:n_val, ]
  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[cfg$m_idx]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[cfg$m_idx]]),
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL),
    optimizer = "L-BFGS-B"
  )
  ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
  dm <- ebmr$ps_fit.list[[1]]$design_matrix
  hx <- ebmr$ps_fit.list[[1]]$h_x
  cat(sprintf("%s: design_mat=%dx%d (%s), h_x=%dx%d (%s)\n",
      cfg$name, nrow(dm), ncol(dm), paste(colnames(dm), collapse=","),
      nrow(hx), ncol(hx), paste(colnames(hx), collapse=",")))
  cat(sprintf("  k=%d params, h_dim=%d moments → overid=%d\n",
      ncol(dm), ncol(hx), ncol(hx) - ncol(dm)))
  cat(sprintf("  Without intercept: h_dim=%d moments → overid=%d\n",
      ncol(hx)-1, ncol(hx)-1 - ncol(dm)))
  cat(sprintf("  Formula: %s\n", deparse(ps_spec[["formula.list"]][[cfg$m_idx]])))
  cat(sprintf("  h_alpha: %s\n\n",
      if(is.function(ps_spec[["h_alpha.list"]][[cfg$m_idx]])) "function"
      else paste(ps_spec[["h_alpha.list"]][[cfg$m_idx]], collapse=",")))
}
