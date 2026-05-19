## Check how similar M2 and M3 fitted PS values are
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EPS", quiet = TRUE)

ps_spec_base <- get_ps_spec("9")
h_alpha_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2)
W_sm <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]]
all_data <- NULL
for (fi in seq_along(data_file)) {
  if (file.exists(data_file[[fi]])) {
    test <- readRDS(data_file[[fi]])
    if (nrow(test) / 1000 >= 2) { all_data <- test; break }
  }
}

nn <- 2000
ps_spec_m23 <- list(
  formula.list = ps_spec_base$formula.list[2:3],
  h_alpha.list = list(h_alpha_fn, h_alpha_fn),
  inv_link = ps_spec_base$inv_link,
  outcome = ps_spec_base$outcome
)

cat("=== M2 vs M3 fitted PS comparison ===\n\n")
cat(sprintf("%-4s | %-10s | %-10s | %-10s | %-10s | %-10s\n",
    "Rep", "cor", "mean|diff|", "max|diff|", "mean_M2", "mean_M3"))
cat(strrep("-", 65), "\n")

cors <- c(); mean_diffs <- c()
for (i in 1:20) {
  dat <- all_data[((i-1)*nn + 1):(i*nn), ]
  ebmr <- EPS$new("y", ps_spec_m23, dat, W_sm)
  ps_m2 <- ebmr$ps_fit.list[[1]]$fitted.values
  ps_m3 <- ebmr$ps_fit.list[[2]]$fitted.values

  cr <- cor(ps_m2, ps_m3)
  md <- mean(abs(ps_m2 - ps_m3))
  mx <- max(abs(ps_m2 - ps_m3))
  cors <- c(cors, cr); mean_diffs <- c(mean_diffs, md)

  cat(sprintf("%4d | %.6f | %.6f   | %.6f   | %.6f   | %.6f\n",
      i, cr, md, mx, mean(ps_m2), mean(ps_m3)))
}

cat(sprintf("\nOverall: mean cor=%.6f, mean |diff|=%.6f\n", mean(cors), mean(mean_diffs)))

# Also check the formulas
cat("\nM2 formula:", deparse(ps_spec_base$formula.list[[2]]), "\n")
cat("M3 formula:", deparse(ps_spec_base$formula.list[[3]]), "\n")
