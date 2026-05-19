setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
  library(EBMRalgorithmFast4)
})

W_func  <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

ps_spec <- get_ps_spec("9-alt1")
ps_sub  <- list(formula.list = ps_spec[["formula.list"]][2],
                h_alpha.list = ps_spec[["h_alpha.list"]][2],
                inv_link     = ps_spec[["inv_link"]],
                outcome      = ps_spec[["outcome"]])

all_data <- readRDS("Simulation_Data/Setting3.B1_n2000_replicate1000.RDS")
n_val <- 2000
n_reps <- 100

cat(sprintf("%-6s %-7s %-10s %-10s %-10s %-10s %-8s\n",
    "Rep", "hess_pd", "eig1", "eig2", "eig3", "eig4", "min/max"))

pd_count <- 0
for (i in 1:n_reps) {
  set.seed(12345 + i)
  dat <- all_data[((i-1)*n_val + 1):(i*n_val), ]

  formula_2 <- ps_spec[["formula.list"]][[2]]
  dm <- model.matrix(formula_2, data=dat)
  h_x <- cbind(1, as.matrix(dat[c("u1","u2","z1","z2")]))
  n <- nrow(dat); h_dim <- ncol(h_x)
  r_col <- dat$r; eta_max <- 20

  Q_full <- function(alpha) {
    eta <- pmin(pmax(as.vector(dm %*% alpha), -eta_max), eta_max)
    pi_vec <- 1 / (1 + exp(eta))
    gm <- (r_col / pi_vec - 1) * h_x
    G <- matrix(colMeans(gm), h_dim, 1)
    W <- tryCatch(solve(t(gm) %*% gm / n), error = function(e) diag(h_dim))
    as.numeric(crossprod(G, W %*% G))
  }

  ebmr <- EBMRAlgorithmFast4$new("y", ps_sub, dat, W_func)
  h_nu_fn <- function(d) cbind(u1=d$u1, u2=d$u2, z1=d$z1, z2=d$z2, u1_u2=d$u1*d$u2)
  ebmr$EBMR_IPW(h_nu=h_nu_fn, type="HT", se.fit=FALSE)
  estimates <- ebmr$ps_fit.list[[1]]$gmm_fit$estimates
  fd_pd <- isTRUE(ebmr$ps_fit.list[[1]]$gmm_fit$opt$hessian_pd)
  if (fd_pd) pd_count <- pd_count + 1

  # Compute eigenvalues at solution
  p <- length(estimates)
  eps <- 1e-4 * max(1, sqrt(sum(estimates^2)))
  f0 <- Q_full(estimates)
  fp <- vapply(seq_len(p), function(j) { x <- estimates; x[j] <- x[j]+eps; Q_full(x) }, numeric(1))
  fm <- vapply(seq_len(p), function(j) { x <- estimates; x[j] <- x[j]-eps; Q_full(x) }, numeric(1))
  H <- matrix(0, p, p)
  for (j in seq_len(p)) H[j,j] <- (fp[j] - 2*f0 + fm[j]) / eps^2
  if (p > 1) for (j in seq_len(p-1)) for (k in (j+1):p) {
    x <- estimates; x[j] <- x[j]+eps; x[k] <- x[k]+eps
    H[j,k] <- H[k,j] <- (Q_full(x) - fp[j] - fp[k] + f0) / eps^2
  }
  eigs <- sort(eigen(H, symmetric=TRUE, only.values=TRUE)$values)
  ratio <- min(eigs) / max(abs(eigs))

  if (!fd_pd) {
    cat(sprintf("rep%3d  %-7s %10.4f %10.4f %10.4f %10.4f %8.5f\n",
        i, fd_pd, eigs[1], eigs[2], eigs[3], eigs[4], ratio))
  }
}

cat(sprintf("\nPD count: %d/%d\n", pd_count, n_reps))
cat("\nNow check threshold sensitivity:\n")
cat("Current threshold: min_eig > 2e-3 * max_abs_eig AND min_eig > 0\n")
cat("A looser threshold (just min_eig > 0, any positive) would classify all\n")
cat("cases with positive min eigenvalue as PD, regardless of ratio.\n")
