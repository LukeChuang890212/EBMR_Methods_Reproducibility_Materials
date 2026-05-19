## Check if W converges in nu estimation for mu_011
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EPS", quiet = TRUE)

ps_spec_base <- get_ps_spec("9-alt1")
h_alpha_fn <- function(dat) cbind(
  u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2)
W_fn <- function(g.matrix) solve(var(g.matrix))
h_nu_fn <- function(dat) cbind(
  u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2,
  u1u2 = dat$u1*dat$u2, z1z2 = dat$z1*dat$z2,
  u1z1 = dat$u1*dat$z1, z1u2 = dat$z1*dat$u2,
  u1z2 = dat$u1*dat$z2, u2z2 = dat$u2*dat$z2)

data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]]
all_data <- NULL
for (fi in seq_along(data_file)) {
  if (file.exists(data_file[[fi]])) {
    test <- readRDS(data_file[[fi]])
    if (nrow(test) / 1000 == 2000) { all_data <- test; break }
  }
}

ps_spec_sub <- list(
  formula.list = ps_spec_base$formula.list[2:3],
  h_alpha.list = list(h_alpha_fn, h_alpha_fn),
  inv_link = ps_spec_base$inv_link,
  outcome = ps_spec_base$outcome
)

cat("=== Checking W convergence in nu estimation for mu_011 ===\n\n")

# Manually replicate nu estimation to track W changes
for (rep_i in 1:5) {
  dat <- all_data[((rep_i-1)*2000 + 1):(rep_i*2000), ]
  ebmr <- EPS$new("y", ps_spec_sub, dat, W_fn)

  # Get PS matrix for models 2,3
  ps_mat <- do.call(cbind, lapply(ebmr$ps_fit.list, function(pf) pf$fitted.values))
  r_vec <- as.vector(dat$r)
  n <- 2000
  h_x <- cbind(1, h_nu_fn(dat))
  h_dim <- ncol(h_x)
  J <- 2

  # Replicate iterative CUE for nu
  g_nu <- function(nu) {
    ps_nu <- as.vector(ps_mat %*% nu)
    as.vector(r_vec / ps_nu - 1) * h_x
  }
  G_nu <- function(nu) colMeans(g_nu(nu))

  nu <- rep(1/J, J)
  cat(sprintf("--- Rep %d ---\n", rep_i))

  # Step 1: W=I
  obj_I <- function(nu) { G <- G_nu(nu); sum(G^2) }
  opt1 <- optim(nu, obj_I, method = "L-BFGS-B", control = list(maxit = 1000))
  nu <- opt1$par

  # Track W changes
  prev_W_norm <- 0
  for (t in 1:30) {
    g_mat <- g_nu(nu)
    W_hat <- tryCatch(solve(var(g_mat)), error = function(e) diag(h_dim))
    W_norm <- norm(W_hat, "F")
    W_change <- abs(W_norm - prev_W_norm) / max(prev_W_norm, 1e-10)

    obj_fn <- function(nu_) { G <- matrix(G_nu(nu_), ncol=1); as.numeric(t(G) %*% W_hat %*% G) }
    opt_t <- optim(nu, obj_fn, method = "L-BFGS-B", control = list(maxit = 1000))
    nu_new <- opt_t$par
    nu_change <- max(abs(nu_new - nu))
    nu <- nu_new
    w <- nu^2 / sum(nu^2)

    if (t <= 10 || t %% 5 == 0) {
      cat(sprintf("  iter %2d: nu=(%.4f,%.4f) w=(%.4f,%.4f) |dnu|=%.4e |dW|=%.4e\n",
          t, nu[1], nu[2], w[1], w[2], nu_change, W_change))
    }
    prev_W_norm <- W_norm
  }
  cat("\n")
}
