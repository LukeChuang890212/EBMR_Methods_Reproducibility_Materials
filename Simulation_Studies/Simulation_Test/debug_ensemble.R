setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

ps_spec_base <- get_ps_spec("9-alt1")
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

h_alpha_rich <- function(dat) cbind(
  u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2,
  u1u2 = dat$u1*dat$u2, u1z1 = dat$u1*dat$z1, u1z2 = dat$u1*dat$z2,
  u2z1 = dat$u2*dat$z1, u2z2 = dat$u2*dat$z2, z1z2 = dat$z1*dat$z2
)
ps_spec <- list(
  formula.list = ps_spec_base$formula.list,
  h_alpha.list = list(h_alpha_rich, h_alpha_rich, h_alpha_rich),
  inv_link = ps_spec_base$inv_link,
  outcome = ps_spec_base$outcome
)
## Try multiple h_nu options
h_nu_rich <- function(dat) cbind(
  1,
  u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2,
  u1u2 = dat$u1*dat$u2, u1z1 = dat$u1*dat$z1, u1z2 = dat$u1*dat$z2,
  u2z1 = dat$u2*dat$z1, u2z2 = dat$u2*dat$z2, z1z2 = dat$z1*dat$z2,
  u1u2z1 = dat$u1*dat$u2*dat$z1
)
h_nu_simple <- function(dat) cbind(
  u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2,
  u1u2 = dat$u1*dat$u2
)
h_nu_fn <- h_nu_rich

data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
alpha_true <- misspecified_model_alpha.true.list[["setting4"]][["miss50"]][[1]]
ps_model_true <- function(dat, alpha.true) {
  X <- cbind(rep(1, nrow(dat)), dat$y, dat$u1, dat$u2)
  1 / (1 + exp(X %*% alpha.true))
}

dat <- all_data[1:2000, ]
cat("Creating EBMR object...\n")
tryCatch({
  ebmr <- EBMRAlgorithmFast4$new("y", ps_spec, dat, W_fn)
  cat("EBMR object created successfully.\n")
  cat("Number of PS fits:", length(ebmr$ps_fit.list), "\n")
  for (j in 1:length(ebmr$ps_fit.list)) {
    pf <- ebmr$ps_fit.list[[j]]
    cat(sprintf("  Model %d: k=%d, h_dim=%d, coef=%s\n",
        j, ncol(pf$design_matrix), ncol(pf$h_x),
        paste(round(pf$coefficients, 3), collapse=", ")))
  }
}, error = function(e) {
  cat("ERROR creating EBMR:", e$message, "\n")
})

## Diagnose h_nu rank issue
cat("\nDiagnosing h_nu W matrix rank...\n")
h_nu_mat <- h_nu_fn(dat)
cat(sprintf("h_nu dimensions: %d x %d\n", nrow(h_nu_mat), ncol(h_nu_mat)))
cat(sprintf("h_nu rank: %d\n", qr(h_nu_mat)$rank))

# Check the ensemble g_mat: (r/ps_nu - 1) * h_nu, where ps_nu = ps.matrix %*% nu
ps2 <- ebmr$ps_fit.list[[2]]$fitted.values
ps3 <- ebmr$ps_fit.list[[3]]$fitted.values
ps_mat <- cbind(ps2, ps3)
cat(sprintf("PS2 range: [%.6f, %.6f]\n", min(ps2), max(ps2)))
cat(sprintf("PS3 range: [%.6f, %.6f]\n", min(ps3), max(ps3)))

# With nu = (0.5, 0.5), ps_nu = average
r <- dat$r
nu_init <- c(0.5, 0.5)
ps_nu <- ps_mat %*% nu_init
g_nu <- as.vector(r / ps_nu - 1) * h_nu_mat
cat(sprintf("g_nu dimensions: %d x %d\n", nrow(g_nu), ncol(g_nu)))
W_nu_cross <- crossprod(g_nu) / nrow(g_nu)
cat(sprintf("W_nu_cross rank: %d (of %d)\n", qr(W_nu_cross)$rank, ncol(W_nu_cross)))
sv <- svd(W_nu_cross)$d
cat(sprintf("W_nu singular values: %s\n", paste(format(sv, digits=3, scientific=TRUE), collapse=", ")))

## Also check: does u1*z1 create linear dependency with u1 and z1?
## u1 is binary, z1 is binary => u1*z1 is binary
## Check if intercept + u1 + z1 + u1*z1 are linearly dependent
test_mat <- cbind(1, dat$u1, dat$z1, dat$u1*dat$z1)
cat(sprintf("\n(1, u1, z1, u1*z1) rank: %d (of %d)\n", qr(test_mat)$rank, ncol(test_mat)))

cat("\nRunning EBMR_IPW with model_indices=c(2,3)...\n")
tryCatch({
  result <- ebmr$EBMR_IPW(
    h_nu = h_nu_fn,
    model_indices = c(2, 3),
    true_ps = ps_model_true(dat, alpha_true),
    type = "HT"
  )
  cat("Result:\n")
  print(unlist(result[1:4]))
}, error = function(e) {
  cat("ERROR in EBMR_IPW:", e$message, "\n")
})
