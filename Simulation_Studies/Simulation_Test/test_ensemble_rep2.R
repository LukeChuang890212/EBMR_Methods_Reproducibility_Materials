setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
source("Basic_setup.r"); source("Data_Generation.r")
source("config/scenarios.R"); source("Simulation.r")
library(EBMRalgorithmFast4)

ps_spec <- get_ps_spec("9-alt1")
ps_spec_13 <- list(
  formula.list = ps_spec[["formula.list"]][c(1,3)],
  h_alpha.list = ps_spec[["h_alpha.list"]][c(1,3)],
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]]
)
W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

set.seed(12347)   # rep 2
dat    <- setting4.B1(2000)
ebmr   <- EBMRAlgorithmFast4$new("y", ps_spec_13, dat, W_func)

# Get PS fits first (no SE for speed)
h_nu <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1_u2=dat$u1*dat$u2)
res  <- ebmr$EBMR_IPW(h_nu=h_nu, type="HT", se.fit=FALSE)

cat("Model 1: alpha_norm=", round(sqrt(sum(ebmr$ps_fit.list[[1]]$coefficients^2)),3),
    " hess_pd=", ebmr$ps_fit.list[[1]]$gmm_fit$opt$hessian_pd,
    " obj=", round(ebmr$ps_fit.list[[1]]$gmm_fit$opt$objective,6), "\n")
cat("Model 2: alpha_norm=", round(sqrt(sum(ebmr$ps_fit.list[[2]]$coefficients^2)),3),
    " hess_pd=", ebmr$ps_fit.list[[2]]$gmm_fit$opt$hessian_pd,
    " obj=", round(ebmr$ps_fit.list[[2]]$gmm_fit$opt$objective,6), "\n")
cat("w.hat:", round(res$w.hat, 4), "\n")
cat("mu_ipw:", round(res$mu_ipw, 4), "\n\n")

# Now manually run ensemble step with diagnostics
ps.matrix <- do.call(cbind, lapply(ebmr$ps_fit.list, function(f) f$fitted.values))
cat("ps.matrix range: [", round(min(ps.matrix),4), ",", round(max(ps.matrix),4), "]\n")
cat("Combined PS at uniform w:", round(range(ps.matrix %*% c(0.5,0.5)), 4), "\n")

# Scan nu2 keeping nu1=1 fixed: what does CUE obj look like?
n   <- nrow(dat)
r   <- as.integer(dat$y != 0 | !is.na(dat$y))  # use r from ebmr
r   <- as.vector(ebmr$.__enclos_env__$private$r)
h_x2 <- h_nu(dat)
h_x  <- cbind(1, h_x2)
h_dim <- ncol(h_x)

Phi_nu <- function(param) {
  ps_nu <- as.vector(ps.matrix %*% param)
  rw    <- r / ps_nu
  (rw - 1) * h_x
}
W_cue <- function(g.mat) solve(crossprod(g.mat) / nrow(g.mat))

cat("\n=== CUE objective scan over nu=(nu1, 1) ===\n")
cat(sprintf("%-8s  %10s\n", "nu1", "CUE obj"))
for (nu1 in c(0.001, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10)) {
  param <- c(nu1, 1)
  g_mat <- Phi_nu(param)
  G     <- colMeans(g_mat)
  W_hat <- tryCatch(W_cue(g_mat), error=function(e) diag(h_dim))
  obj   <- as.numeric(t(G) %*% W_hat %*% G)
  w2    <- 1 / (nu1^2 + 1)
  cat(sprintf("%-8.4f  %10.6f  (w2=%.4f)\n", nu1, obj, w2))
}

cat("\n=== CUE objective scan over nu2 (nu1=1 fixed) ===\n")
cat(sprintf("%-8s  %10s\n", "nu2", "CUE obj"))
for (nu2 in c(0.01, 0.1, 0.3, 0.5, 0.7, 1, 2, 3, 5, 7, 10)) {
  param <- c(1, nu2)
  g_mat <- Phi_nu(param)
  G     <- colMeans(g_mat)
  W_hat <- tryCatch(W_cue(g_mat), error=function(e) diag(h_dim))
  obj   <- as.numeric(t(G) %*% W_hat %*% G)
  w2    <- nu2^2 / (1 + nu2^2)
  cat(sprintf("%-8.3f  %10.6f  (w2=%.3f)\n", nu2, obj, w2))
}
