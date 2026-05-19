setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
source("Simulation.r")
library(EPS)

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

rep_i <- 1
dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]

ebmr <- EPS$new("y", subset_ps_spec, dat, W_func)
res <- ebmr$EBMR_IPW(
  h_nu = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2),
  se.fit = TRUE, type = "HT"
)
cat(sprintf("mu_ipw = %.4f, se_ipw = %.4f\n", res$mu_ipw, res$se_ipw))

# Check what fields are available
fit <- ebmr$ps_fit.list[[1]]
cat("ps_fit fields: ", paste(names(fit), collapse=", "), "\n")

# Check for ensemble_ps
cat("Has ensemble_ps field: ", "ensemble_ps" %in% names(res), "\n")
cat("res fields: ", paste(names(res), collapse=", "), "\n")
