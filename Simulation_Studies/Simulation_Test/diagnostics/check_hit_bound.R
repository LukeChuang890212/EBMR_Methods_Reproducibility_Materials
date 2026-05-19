setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
  devtools::load_all("../EPS")
})
W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
n_val <- 2000; n_reps <- 200

ps_spec <- get_ps_spec("9-alt1")
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)

ps_sub3 <- list(
  formula.list = ps_spec[["formula.list"]][3],
  h_alpha.list = ps_spec[["h_alpha.list"]][3],
  inv_link     = ps_spec[["inv_link"]],
  outcome      = ps_spec[["outcome"]]
)

n_hit <- 0
for (i in 1:n_reps) {
  dat <- all_data[((i-1)*n_val + 1):(i*n_val), ]
  tryCatch({
    ebmr <- EPS[["new"]]("y", ps_sub3, dat, W_func)
    alpha <- ebmr[["ps_fit.list"]][[1]][["gmm_fit"]][["estimates"]]
    if (any(abs(alpha) > 9.99)) {
      n_hit <- n_hit + 1
      cat(sprintf("Rep %3d: alpha=(%s)\n", i, paste(round(alpha, 4), collapse=", ")))
    }
  }, error = function(e) NULL)
}
cat(sprintf("\nTotal reps hitting bound (|alpha| > 9.99): %d / %d\n", n_hit, n_reps))
