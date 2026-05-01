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

mu_vals <- cond_vals <- rep(NA_real_, n_reps)
for (rep_i in 1:n_reps) {
  dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
  single_ps <- list(
    formula.list = list(ps_spec[["formula.list"]][[3]]),
    h_alpha.list = list(ps_spec[["h_alpha.list"]][[3]]),
    inv_link = ps_spec[["inv_link"]],
    outcome = ps_spec[["outcome"]],
    alpha_init.list = list(NULL),
    optimizer = "constrained_nr"
  )
  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
    pi_hat <- ebmr$ps_fit.list[[1]]$fitted.values
    mu_vals[rep_i] <- mean(dat[["r"]] * dat[["y"]] / pi_hat)
    cond_vals[rep_i] <- ebmr$ps_fit.list[[1]]$gmm_fit$opt$solution_cond
  }, error = function(e) NULL)
}

valid <- !is.na(mu_vals)
mu <- mu_vals[valid]

Q1 <- quantile(mu, 0.25)
Q3 <- quantile(mu, 0.75)
IQR_val <- Q3 - Q1
lower <- Q1 - 3 * IQR_val
upper <- Q3 + 3 * IQR_val

cat(sprintf("Q1=%.4f, Q3=%.4f, IQR=%.4f\n", Q1, Q3, IQR_val))
cat(sprintf("3IQR bounds: [%.4f, %.4f]\n", lower, upper))
cat(sprintf("# below lower: %d\n", sum(mu < lower)))
cat(sprintf("# above upper: %d\n", sum(mu > upper)))

# Also check 1.5IQR
lower15 <- Q1 - 1.5 * IQR_val
upper15 <- Q3 + 1.5 * IQR_val
cat(sprintf("\n1.5IQR bounds: [%.4f, %.4f]\n", lower15, upper15))
cat(sprintf("# below lower: %d\n", sum(mu < lower15)))
cat(sprintf("# above upper: %d\n", sum(mu > upper15)))

# Show the degenerate reps
degen <- which(valid)[cond_vals[valid] >= 1e8]
cat(sprintf("\nDegenerate reps (cond>=1e8): %d\n", length(degen)))
cat(sprintf("Their mu values: %s\n", paste(round(mu_vals[degen], 4), collapse=", ")))
cat(sprintf("Min normal mu: %.4f\n", min(mu_vals[valid & cond_vals < 1e8], na.rm=TRUE)))
cat(sprintf("Max degenerate mu: %.4f\n", max(mu_vals[degen])))
