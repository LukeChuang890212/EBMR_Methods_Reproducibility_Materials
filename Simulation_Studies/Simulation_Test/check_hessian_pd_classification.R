setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
  library(EBMRalgorithmFast4)
  library(numDeriv)
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

cat("Comparing FD hessian (current) vs Richardson hessian for PD classification\n")
cat(sprintf("%-6s %-12s %-12s %-10s %-10s\n", "Rep", "FD_min_eig", "Rich_min_eig", "FD_pd", "Rich_pd"))

fd_pd_count <- rich_pd_count <- 0

for (i in 1:n_reps) {
  set.seed(12345 + i)
  dat <- all_data[((i-1)*n_val + 1):(i*n_val), ]

  ebmr <- EBMRAlgorithmFast4$new("y", ps_sub, dat, W_func)
  h_nu_fn <- function(d) cbind(u1=d$u1, u2=d$u2, z1=d$z1, z2=d$z2, u1_u2=d$u1*d$u2)
  ebmr$EBMR_IPW(h_nu=h_nu_fn, type="HT", se.fit=FALSE)

  fit <- ebmr$ps_fit.list[[1]]$gmm_fit
  estimates <- fit$estimates
  fd_pd <- isTRUE(fit$opt$hessian_pd)
  if (fd_pd) fd_pd_count <- fd_pd_count + 1

  # Rebuild Q_full for Richardson
  r_col <- dat$r
  formula_2 <- ps_spec[["formula.list"]][[2]]
  dm <- model.matrix(formula_2, data=dat)
  h_x <- cbind(1, as.matrix(dat[c("u1","u2","z1","z2")]))
  n <- nrow(dat); h_dim <- ncol(h_x)
  eta_max <- 20

  Q_full <- function(alpha) {
    eta <- pmin(pmax(as.vector(dm %*% alpha), -eta_max), eta_max)
    pi_vec <- 1 / (1 + exp(eta))
    gm <- (r_col / pi_vec - 1) * h_x
    G <- matrix(colMeans(gm), h_dim, 1)
    W <- tryCatch(solve(t(gm) %*% gm / n), error = function(e) diag(h_dim))
    as.numeric(crossprod(G, W %*% G))
  }

  H_rich <- tryCatch(hessian(Q_full, estimates), error = function(e) matrix(NA, 4, 4))
  eig_rich <- tryCatch(eigen(H_rich, symmetric=TRUE, only.values=TRUE)$values, error=function(e) rep(NA, 4))
  min_rich <- min(eig_rich)
  max_rich <- max(abs(eig_rich))
  rich_pd <- !is.na(min_rich) && min_rich > 2e-3 * max_rich && min_rich > 0
  if (rich_pd) rich_pd_count <- rich_pd_count + 1

  if (fd_pd != rich_pd) {
    cat(sprintf("rep%3d  FD=%-5s  Rich=%-5s  (DISAGREE)\n", i, fd_pd, rich_pd))
  }
}

cat(sprintf("\nFD hessian:        PD in %d/%d reps\n", fd_pd_count, n_reps))
cat(sprintf("Richardson hessian: PD in %d/%d reps\n", rich_pd_count, n_reps))
