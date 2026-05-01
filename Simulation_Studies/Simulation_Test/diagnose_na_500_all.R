setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Data_Generation.r")
source("config/scenarios.R")
library(EBMRalgorithmFast4)
source("Basic_setup.r")

data_file <- misspecified_model_all_data_file.list$setting3$miss50[[2]]
all_data <- readRDS(data_file)
n_val <- 500

ps_spec <- get_ps_spec("9-alt1")
subset_ps_spec <- list(
  formula.list = ps_spec$formula.list[2],
  h_alpha.list = ps_spec$h_alpha.list[2],
  inv_link = ps_spec$inv_link,
  outcome = ps_spec$outcome
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

na_reps <- c(93, 140, 327, 475, 577, 632, 667, 675, 683, 734, 821, 863, 890, 915, 970, 984)

cat(sprintf("%-6s %-10s %-10s %-8s\n", "Rep", "mu_ipw", "se_ipw", "Status"))
for (rep_i in na_reps) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]
  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", subset_ps_spec, dat, W_func)
    res <- ebmr$EBMR_IPW(
      h_nu = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2),
      se.fit = TRUE, type = "HT"
    )
    cat(sprintf("%-6d %-10.4f %-10.4f OK\n", rep_i, res$mu_ipw, res$se_ipw))
  }, error = function(e) {
    cat(sprintf("%-6d %-10s %-10s ERROR: %s\n", rep_i, "NA", "NA", conditionMessage(e)))
  })
}

cat("\nDone!\n")
