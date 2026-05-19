## Check alpha stability across reps
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

ps_spec_base <- get_ps_spec("7")
h_alpha_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2)
W_sm <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

data_file <- correct_model_all_data_file.list[["setting3"]][["miss50"]]
all_data <- NULL
for (fi in seq_along(data_file)) {
  if (file.exists(data_file[[fi]])) {
    test <- readRDS(data_file[[fi]])
    if (nrow(test) / 1000 == 2000) { all_data <- test; break }
  }
}

ps_spec_m2 <- list(
  formula.list = list(ps_spec_base$formula.list[[2]]),
  h_alpha.list = list(h_alpha_fn),
  inv_link = ps_spec_base$inv_link,
  outcome = ps_spec_base$outcome
)

nn <- 2000
cat("M2:", deparse(ps_spec_base$formula.list[[2]]), "\n\n")

# Collect alpha, PS min, mu for 100 reps
alpha_mat <- matrix(NA, 100, 4)
ps_min_vec <- numeric(100)
mu_vec <- numeric(100)

for (i in 1:100) {
  dat <- all_data[((i-1)*nn + 1):(i*nn), ]
  ebmr <- EBMRAlgorithmFast4$new("y", ps_spec_m2, dat, W_sm)
  alpha_mat[i,] <- ebmr$ps_fit.list[[1]]$coefficients
  ps_min_vec[i] <- min(ebmr$ps_fit.list[[1]]$fitted.values)
  mu_vec[i] <- mean(dat$r * dat$y / ebmr$ps_fit.list[[1]]$fitted.values)
}

cat("Alpha distribution across 100 reps:\n")
for (k in 1:4) {
  cat(sprintf("  alpha%d: mean=%7.4f sd=%7.4f min=%7.4f max=%7.4f\n",
      k, mean(alpha_mat[,k]), sd(alpha_mat[,k]), min(alpha_mat[,k]), max(alpha_mat[,k])))
}

cat(sprintf("\nPS_min: mean=%.4f sd=%.4f min=%.4f max=%.4f\n",
    mean(ps_min_vec), sd(ps_min_vec), min(ps_min_vec), max(ps_min_vec)))
cat(sprintf("mu_ipw: mean=%.4f sd=%.4f min=%.4f max=%.4f\n",
    mean(mu_vec), sd(mu_vec), min(mu_vec), max(mu_vec)))

# Identify extreme reps
extreme <- which(abs(mu_vec - mean(mu_vec)) > 2 * sd(mu_vec))
cat(sprintf("\nExtreme reps (|mu - mean| > 2*sd): %d\n", length(extreme)))
if (length(extreme) > 0) {
  for (i in extreme) {
    cat(sprintf("  Rep %3d: alpha=(%s) PS_min=%.4f mu=%.4f\n",
        i, paste(round(alpha_mat[i,], 4), collapse=", "), ps_min_vec[i], mu_vec[i]))
  }
}

# Show normal reps for comparison
normal <- which(abs(mu_vec - mean(mu_vec)) < 0.5 * sd(mu_vec))
cat(sprintf("\nNormal reps (first 5):\n"))
for (i in head(normal, 5)) {
  cat(sprintf("  Rep %3d: alpha=(%s) PS_min=%.4f mu=%.4f\n",
      i, paste(round(alpha_mat[i,], 4), collapse=", "), ps_min_vec[i], mu_vec[i]))
}
