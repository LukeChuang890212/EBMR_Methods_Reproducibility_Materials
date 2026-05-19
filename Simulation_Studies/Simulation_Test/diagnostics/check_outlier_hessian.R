setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
  library(EPS)
})

W_func  <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_nu_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1_u2=dat$u1*dat$u2)

ps_spec <- get_ps_spec("9-alt1")
ps_sub  <- list(formula.list = ps_spec[["formula.list"]][2],
                h_alpha.list = ps_spec[["h_alpha.list"]][2],
                inv_link     = ps_spec[["inv_link"]],
                outcome      = ps_spec[["outcome"]])
mu_true <- get_mu_true("setting3")
all_data <- readRDS("Simulation_Data/Setting3.B1_n2000_replicate1000.RDS")
n_val <- 2000

# Outlier reps identified from diagnose_outliers_m2.R (|mu - mu_true| > 3*ESD over 500 reps)
outlier_reps <- c(125, 154, 180, 216, 221, 232, 255, 263, 285, 308, 314, 353, 420, 430, 438)
# PD=TRUE: 216, 221, 232, 255, 285, 353, 420, 430
# PD=FALSE: 125, 154, 180, 263, 308, 314, 438

formula_2 <- ps_spec[["formula.list"]][[2]]
eta_max   <- 20

Q_fn <- function(alpha, dm, r_col, h_x) {
  n     <- nrow(dm)
  h_dim <- ncol(h_x)
  eta   <- pmin(pmax(as.vector(dm %*% alpha), -eta_max), eta_max)
  pi_v  <- 1 / (1 + exp(eta))
  gm    <- (r_col / pi_v - 1) * h_x
  G     <- matrix(colMeans(gm), h_dim, 1)
  W     <- tryCatch(solve(t(gm) %*% gm / n), error = function(e) diag(h_dim))
  as.numeric(crossprod(G, W %*% G))
}

hess_eigs <- function(alpha, dm, r_col, h_x) {
  p   <- length(alpha)
  eps <- 1e-4 * max(1, sqrt(sum(alpha^2)))
  f0  <- Q_fn(alpha, dm, r_col, h_x)
  fp  <- vapply(seq_len(p), function(j) { x <- alpha; x[j] <- x[j]+eps; Q_fn(x, dm, r_col, h_x) }, numeric(1))
  fm  <- vapply(seq_len(p), function(j) { x <- alpha; x[j] <- x[j]-eps; Q_fn(x, dm, r_col, h_x) }, numeric(1))
  H   <- matrix(0, p, p)
  for (j in seq_len(p)) H[j,j] <- (fp[j] - 2*f0 + fm[j]) / eps^2
  if (p > 1) for (j in seq_len(p-1)) for (k in (j+1):p) {
    x <- alpha; x[j] <- x[j]+eps; x[k] <- x[k]+eps
    H[j,k] <- H[k,j] <- (Q_fn(x, dm, r_col, h_x) - fp[j] - fp[k] + f0) / eps^2
  }
  sort(eigen(H, symmetric=TRUE, only.values=TRUE)$values)
}

cat(sprintf("%-5s %-8s %-8s %-10s %-10s %-10s %-10s %-8s %-8s %s\n",
    "Rep", "hpd", "mu", "eig1", "eig2", "eig3", "eig4", "ratio", "norm_a", "alpha"))

for (i in outlier_reps) {
  set.seed(12345 + i)
  dat  <- all_data[((i-1)*n_val + 1):(i*n_val), ]
  ebmr <- EPS$new("y", ps_sub, dat, W_func)
  res  <- ebmr$EBMR_IPW(h_nu=h_nu_fn, type="HT", se.fit=FALSE)
  fit  <- ebmr$ps_fit.list[[1]]$gmm_fit
  hpd  <- isTRUE(fit$opt$hessian_pd)
  alpha <- fit$estimates

  dm    <- model.matrix(formula_2, data=dat)
  h_x   <- cbind(1, as.matrix(dat[c("u1","u2","z1","z2")]))
  r_col <- dat$r

  eigs  <- hess_eigs(alpha, dm, r_col, h_x)
  ratio <- min(eigs) / max(abs(eigs))
  anorm <- sqrt(sum(alpha^2))

  cat(sprintf("%-5d %-8s %8.3f %10.4f %10.4f %10.4f %10.4f %8.5f %8.3f [%s]\n",
      i, hpd, res$mu_ipw,
      eigs[1], eigs[2], eigs[3], eigs[4],
      ratio, anorm,
      paste(round(alpha, 2), collapse=", ")))
}
