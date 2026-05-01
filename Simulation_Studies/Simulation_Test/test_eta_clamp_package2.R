## Quick test: single rep with package to check if eta clamping works
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

ps_spec <- get_ps_spec("9")
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_nu <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1u2=dat$u1*dat$u2)

all_data <- readRDS("Simulation_Data/Setting3.B2_n500_replicate1000.RDS")

# Test on the reps that previously gave "non-finite value" errors
na_reps <- c(100, 139, 315, 351, 412, 417, 423, 443, 556, 585, 602, 822, 847, 852, 968)

m_idx <- 2
single_ps <- list(
  formula.list = list(ps_spec[["formula.list"]][[m_idx]]),
  h_alpha.list = list(ps_spec[["h_alpha.list"]][[m_idx]]),
  inv_link = ps_spec[["inv_link"]], outcome = ps_spec[["outcome"]]
)

cat("Testing previously-NA reps with eta clamping in package:\n\n")
for (rep_i in na_reps) {
  dat <- all_data[((rep_i-1)*500 + 1):(rep_i*500), ]
  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
    result <- ebmr$EBMR_IPW(h_nu = h_nu, true_ps = NULL, type = "HT")
    cat(sprintf("Rep %3d: mu=%.4f, se=%s\n", rep_i, result$mu_ipw,
        if(is.na(result$se_ipw)) "NA" else sprintf("%.4f", result$se_ipw)))
  }, error = function(e) {
    cat(sprintf("Rep %3d: ERROR: %s\n", rep_i, e$message))
  })
}
