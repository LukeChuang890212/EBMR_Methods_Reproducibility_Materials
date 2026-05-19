## Compare CUE objective at the two nu basins for mu_011
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

ps_spec_base <- get_ps_spec("9-alt1")
h_alpha_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2)
W_fn <- function(g.matrix) solve(var(g.matrix))
h_nu_fn <- function(dat) cbind(
  u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2,
  u1u2=dat$u1*dat$u2, z1z2=dat$z1*dat$z2,
  u1z1=dat$u1*dat$z1, z1u2=dat$z1*dat$u2,
  u1z2=dat$u1*dat$z2, u2z2=dat$u2*dat$z2)

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

cat("=== Comparing two nu basins for mu_011 ===\n\n")

for (rep_i in 1:10) {
  dat <- all_data[((rep_i-1)*2000 + 1):(rep_i*2000), ]
  ebmr <- EBMRAlgorithmFast4$new("y", ps_spec_sub, dat, W_fn)

  ps_mat <- do.call(cbind, lapply(ebmr$ps_fit.list, function(pf) pf$fitted.values))
  r_vec <- as.vector(dat$r)
  n <- 2000
  h_x <- cbind(1, h_nu_fn(dat))
  h_dim <- ncol(h_x)

  g_nu <- function(nu) {
    ps_nu <- as.vector(ps_mat %*% nu)
    as.vector(r_vec / ps_nu - 1) * h_x
  }
  G_nu <- function(nu) colMeans(g_nu(nu))
  cue_obj <- function(nu) {
    g_mat <- g_nu(nu)
    G <- colMeans(g_mat)
    W <- tryCatch(solve(var(g_mat)), error = function(e) diag(h_dim))
    as.numeric(t(G) %*% W %*% G)
  }

  # Try from two inits: (1,0) and (0,1)
  run_from <- function(init) {
    nu <- init
    for (t in 1:50) {
      g_mat <- g_nu(nu)
      W_hat <- tryCatch(solve(var(g_mat)), error = function(e) diag(h_dim))
      obj_fn <- function(nu_) { G <- matrix(G_nu(nu_), ncol=1); as.numeric(t(G) %*% W_hat %*% G) }
      opt <- optim(nu, obj_fn, method = "L-BFGS-B", control = list(maxit = 1000))
      nu <- opt$par
    }
    w <- nu^2 / sum(nu^2)
    obj <- cue_obj(nu)
    mu <- mean(r_vec * dat$y / as.vector(ps_mat %*% w))
    list(nu = nu, w = w, obj = obj, mu = mu)
  }

  res_M2 <- run_from(c(1, 0.01))  # start near M2
  res_M3 <- run_from(c(0.01, 1))  # start near M3

  cat(sprintf("Rep %2d: M2-dom: w=(%.4f,%.4f) obj=%.6e mu=%.4f | M3-dom: w=(%.4f,%.4f) obj=%.6e mu=%.4f | obj_ratio=%.3f\n",
      rep_i,
      res_M2$w[1], res_M2$w[2], res_M2$obj, res_M2$mu,
      res_M3$w[1], res_M3$w[2], res_M3$obj, res_M3$mu,
      res_M2$obj / res_M3$obj))
}
