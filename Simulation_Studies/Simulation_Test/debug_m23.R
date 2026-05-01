setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
  library(EBMRalgorithmFast4)
})

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_nu_fn <- function(dat) cbind(u1=dat[["u1"]], u2=dat[["u2"]], z1=dat[["z1"]], z2=dat[["z2"]])
n_val <- 2000

ps_spec <- get_ps_spec("9-alt1")
ps_sub <- list(
  formula.list = ps_spec[["formula.list"]][2:3],
  h_alpha.list = ps_spec[["h_alpha.list"]][2:3],
  inv_link     = ps_spec[["inv_link"]],
  outcome      = ps_spec[["outcome"]]
)

data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
cat(sprintf("Data file: %s\n", data_file))
all_data <- readRDS(data_file)
dat <- all_data[1:n_val, ]

cat("Creating EBMR...\n")
ebmr <- EBMRAlgorithmFast4[["new"]]("y", ps_sub, dat, W_func)
cat("PS models fitted OK\n")

# Check ps_fit dimensions
for (j in 1:2) {
  gf <- ebmr[["ps_fit.list"]][[j]][["gmm_fit"]]
  cat(sprintf("Model %d: alpha_dim=%d, psi dim=%s, psi2 dim=%s\n",
    j, length(gf[["estimates"]]),
    paste(dim(gf[["psi"]]), collapse="x"),
    if (is.matrix(gf[["psi2"]])) paste(dim(gf[["psi2"]]), collapse="x") else "NULL/NA"))
}

cat("\nRunning EBMR_IPW...\n")
tryCatch({
  res <- ebmr[["EBMR_IPW"]](h_nu=h_nu_fn, type="HT", se.fit=TRUE)
  cat(sprintf("mu_ipw=%.4f, se_ipw=%.4f\n", res[["mu_ipw"]], res[["se_ipw"]]))
  cat(sprintf("w.hat = %s\n", paste(round(res[["w.hat"]], 4), collapse=", ")))
}, error = function(e) {
  cat(sprintf("Error: %s\n", conditionMessage(e)))
})
