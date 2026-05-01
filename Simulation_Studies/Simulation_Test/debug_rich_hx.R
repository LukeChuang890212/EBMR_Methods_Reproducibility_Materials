setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({ source("Basic_setup.r"); source("Data_Generation.r"); source("config/scenarios.R"); source("Simulation.r") })
devtools::load_all("../EBMRalgorithmFast4", quiet=TRUE)
n_val <- 2000; ps_spec <- get_ps_spec("9-alt1")
W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
compute_pi_fn <- function(eta) plogis(-eta)
data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
dat <- all_data[1:n_val, ]
single_ps <- list(formula.list=list(ps_spec[["formula.list"]][[3]]),
  h_alpha.list=list(ps_spec[["h_alpha.list"]][[3]]),
  inv_link=ps_spec[["inv_link"]], outcome=ps_spec[["outcome"]],
  alpha_init.list=list(NULL), optimizer="L-BFGS-B")
ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn)
h_x_orig <- ebmr$ps_fit.list[[1]]$h_x
h_x_rich <- cbind(h_x_orig, u1sq=dat$u1^2, u2sq=dat$u2^2, z1sq=dat$z1^2, z2sq=dat$z2^2, u1u2=dat$u1*dat$u2)
design_mat <- ebmr$ps_fit.list[[1]]$design_matrix
cat("h_x_rich dim:", dim(h_x_rich), "\n")

# Run the optimizer for one rep
alpha_hat <- c(0.5, -2, 0.3, 0.1)  # moderate alpha
pi_hat <- compute_pi_fn(as.vector(design_mat %*% alpha_hat))
r_vec <- dat$r; y_vec <- dat$y; n <- n_val
g_mat <- (r_vec / pi_hat - 1) * h_x_rich
h <- ncol(h_x_rich); k <- ncol(design_mat)

# Check W
Omega <- crossprod(g_mat)/n
cat("Omega cond:", kappa(Omega), "\n")
W_hat <- tryCatch(solve(Omega), error=function(e) { cat("ERROR:", e$message, "\n"); NULL })
if (is.null(W_hat)) { cat("W_hat failed\n"); q() }

# Check H_mat
cf <- r_vec * (1-pi_hat)/pi_hat
Gamma_hat <- crossprod(h_x_rich * cf, design_mat)/n
GtWG <- crossprod(Gamma_hat, W_hat %*% Gamma_hat)
GtW <- crossprod(Gamma_hat, W_hat)
G_vec <- colMeans(g_mat)
eta_s <- matrix(G_vec, h, 1)
WG <- as.vector(W_hat %*% eta_s)
GW_row <- as.vector(t(eta_s) %*% W_hat)

dc_X <- cf * design_mat
R_mat <- matrix(0, h*k, k)
for (j in 1:k) {
  block <- crossprod(dc_X, design_mat[,j]*h_x_rich)/n
  for (l in 1:h) R_mat[(l-1)*k+j,] <- block[,l]
}
S_mat <- matrix(0, h^2, k)
for (j in 1:k) {
  dg_j <- cf*design_mat[,j]*h_x_rich
  cross <- (crossprod(dg_j,g_mat)+crossprod(g_mat,dg_j))/n
  S_mat[,j] <- as.vector(cross)
}
GW_kron_Ik <- kronecker(matrix(GW_row, 1, h), diag(k))
GW_kron_GtW <- kronecker(matrix(GW_row, 1, h), GtW)
H_mat <- GtWG + GW_kron_Ik %*% R_mat - GW_kron_GtW %*% S_mat

cat("H_mat:\n"); print(round(H_mat, 6))
ev <- eigen(H_mat, symmetric=FALSE, only.values=TRUE)$values
cat("H_mat eigenvalues:", ev, "\n")
H_cond <- max(Mod(ev))/max(min(Mod(ev)), 1e-15)
cat("H_mat cond:", H_cond, "\n")
cat("H_cond > 1e15?", H_cond > 1e15, "\n")
