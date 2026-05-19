## Check damped nu point estimates for 20 reps
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

ps_spec_base <- get_ps_spec("9")
h_alpha_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2)
W_sm <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_nu_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1u2=dat$u1*dat$u2)

data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]]
all_data <- readRDS(data_file[[1]])
nn <- 2000

ps_spec_m23 <- list(
  formula.list = ps_spec_base$formula.list[2:3],
  h_alpha.list = list(h_alpha_fn, h_alpha_fn),
  inv_link = ps_spec_base$inv_link,
  outcome = ps_spec_base$outcome
)

cat("=== Damped W nu vs package CUE nu, 20 reps ===\n\n")
cat(sprintf("%-5s | %-25s | %-25s\n", "Rep", "Damped W", "Package CUE"))
cat(strrep("-", 60), "\n")

for (i in 1:20) {
  dat <- all_data[((i-1)*nn + 1):(i*nn), ]
  ebmr <- EBMRAlgorithmFast4$new("y", ps_spec_m23, dat, W_sm)
  ps_mat <- do.call(cbind, lapply(ebmr$ps_fit.list, function(pf) pf$fitted.values))
  r_vec <- as.vector(dat$r)
  h_x <- cbind(1, h_nu_fn(dat))
  h_dim <- ncol(h_x); n <- nn

  g_nu <- function(nu) {
    w <- nu^2 / sum(nu^2)
    ps_nu <- as.vector(ps_mat %*% w)
    as.vector(r_vec / ps_nu - 1) * h_x
  }
  G_nu <- function(nu) colMeans(g_nu(nu))

  # Damped W
  nu <- c(0.5, 0.5)
  lambda <- 0.5; W_prev <- diag(h_dim); prev_gn <- Inf
  obj_I <- function(nu) { G <- G_nu(nu); sum(G^2) }
  opt1 <- optim(nu, obj_I, method = "L-BFGS-B", control = list(maxit = 500))
  nu <- opt1$par
  g_mat <- g_nu(nu)
  W_prev <- tryCatch(solve(t(g_mat) %*% g_mat / n), error = function(e) diag(h_dim))

  for (t in 1:200) {
    g_mat <- g_nu(nu)
    W_new <- tryCatch(solve(t(g_mat) %*% g_mat / n), error = function(e) diag(h_dim))
    W_hat <- lambda * W_prev + (1 - lambda) * W_new
    W_prev <- W_hat
    nu_old <- nu
    obj_fn <- function(nu) { G <- matrix(G_nu(nu), ncol=1); as.numeric(t(G) %*% W_hat %*% G) }
    opt <- optim(nu, obj_fn, method = "L-BFGS-B", control = list(maxit = 500))
    nu <- opt$par
    gn <- max(abs(nu - nu_old))
    if (gn >= prev_gn) { lambda <- min(lambda + 0.1, 0.95) } else { lambda <- max(lambda - 0.02, 0.0) }
    prev_gn <- gn
    if (gn < 1e-6) break
  }
  w_damped <- nu^2 / sum(nu^2)

  # Package CUE
  ipw <- ebmr$EBMR_IPW(h_nu_fn, model_indices = 1:2, se.fit = FALSE)

  cat(sprintf("%5d | w=(%.3f,%.3f)         | w=(%.3f,%.3f)\n",
      i, w_damped[1], w_damped[2], ipw$w.hat[1], ipw$w.hat[2]))
}
