setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
  library(EBMRalgorithmFast4)
})

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_nu_fn <- function(dat) cbind(u1=dat[["u1"]], u2=dat[["u2"]], z1=dat[["z1"]], z2=dat[["z2"]])
n_val <- 2000

ps_spec <- get_ps_spec("9-alt1")
ps_sub <- list(
  formula.list = ps_spec[["formula.list"]][2:3],
  h_alpha.list = ps_spec[["h_alpha.list"]][2:3],
  inv_link     = ps_spec[["inv_link"]],
  outcome      = ps_spec[["outcome"]]
)

data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
dat <- all_data[1:n_val, ]

ebmr <- EBMRAlgorithmFast4[["new"]]("y", ps_sub, dat, W_func)
res  <- ebmr[["EBMR_IPW"]](h_nu=h_nu_fn, type="HT", se.fit=TRUE)

J2 <- 2
r_ <- dat[["r"]]; y_ <- dat[["y"]]
nu_hat_ <- res[["nu.hat"]]; w_hat_ <- res[["w.hat"]]
ps_mat_ <- res[["ps.matrix"]]
ens_ps_ <- as.vector(ps_mat_ %*% w_hat_)
alpha_dim_ <- sapply(1:J2, function(j) length(ebmr[["ps_fit.list"]][[j]][["gmm_fit"]][["estimates"]]))

cat(sprintf("alpha_dim: %s\n", paste(alpha_dim_, collapse=",")))
cat(sprintf("nu_hat: %s\n", paste(round(nu_hat_, 4), collapse=",")))
cat(sprintf("w_hat: %s\n", paste(round(w_hat_, 4), collapse=",")))

# dot_pi
dot_pi_full <- matrix(0, n_val, sum(alpha_dim_))
for (j2 in 1:J2) {
  pi_j <- ebmr[["ps_fit.list"]][[j2]][["fitted.values"]]
  des_j <- ebmr[["ps_fit.list"]][[j2]][["design_matrix"]]
  lt_j <- ebmr[["ps_fit.list"]][[j2]][["link_type"]]
  cols_j <- (sum(alpha_dim_[0:(j2-1)])+1):sum(alpha_dim_[1:j2])
  cat(sprintf("Model %d: link=%s, cols=%s\n", j2, lt_j, paste(cols_j, collapse=",")))
  if (lt_j == "logistic") {
    dot_pi_full[, cols_j] <- des_j * (pi_j * (1 - pi_j))
  } else {
    dot_pi_full[, cols_j] <- -des_j * (pi_j * (1 - pi_j))
  }
}

ry_ps2_ <- as.vector(r_ * y_ * (ens_ps_^(-2)))
H_alpha_ <- colMeans(t(t(dot_pi_full) * rep(w_hat_, alpha_dim_)) * ry_ps2_)
cat(sprintf("H_alpha dim: %d\n", length(H_alpha_)))

psi2_list <- lapply(1:J2, function(j) ebmr[["ps_fit.list"]][[j]][["gmm_fit"]][["psi2"]])
for (j2 in 1:J2) cat(sprintf("psi2 model %d: %s\n", j2, paste(dim(psi2_list[[j2]]), collapse="x")))

psi2_alpha <- do.call(rbind, psi2_list)
cat(sprintf("psi2_alpha: %s\n", paste(dim(psi2_alpha), collapse="x")))

iid_ <- as.vector(r_ / ens_ps_ * y_) - as.vector(t(H_alpha_) %*% psi2_alpha)
cat(sprintf("iid_ length: %d\n", length(iid_)))

# nu adjustment
dot_W_ <- function(nu) {
  nu <- as.vector(nu)
  (diag(2*nu)*sum(nu^2) - 2*nu %*% t(nu^2)) / (sum(nu^2)^2)
}
dot_W_nu_ <- dot_W_(nu_hat_)
cat(sprintf("dot_W_nu dim: %s\n", paste(dim(dot_W_nu_), collapse="x")))

wH_nu_ <- colMeans(ps_mat_ %*% t(dot_W_nu_) * ry_ps2_)
cat(sprintf("wH_nu length: %d\n", length(wH_nu_)))

ens_fit <- res[["ensemble_fit"]]
psi_nu_ <- ens_fit[["gmm_fit"]][["psi"]]
Gamma_nu_ <- ens_fit[["gmm_fit"]][["Gamma.hat"]]
W_nu_ <- ens_fit[["gmm_fit"]][["W.hat"]]

cat(sprintf("psi_nu: %s\n", paste(dim(psi_nu_), collapse="x")))
cat(sprintf("Gamma_nu: %s\n", paste(dim(Gamma_nu_), collapse="x")))
cat(sprintf("W_nu: %s\n", paste(dim(W_nu_), collapse="x")))

h_nu_mat_ <- h_nu_fn(dat)
ps_nu_ <- as.vector(ps_mat_ %*% nu_hat_)
r_ps2_nu_ <- as.vector(-r_ * (ps_nu_^(-2)))

cat(sprintf("h_nu_mat dim: %s\n", paste(dim(h_nu_mat_), collapse="x")))
cat(sprintf("dot_pi_full dim: %s\n", paste(dim(dot_pi_full), collapse="x")))
cat(sprintf("nu_hat: %s, alpha_dim: %s\n", paste(round(nu_hat_,4), collapse=","), paste(alpha_dim_, collapse=",")))

# This line likely fails
tryCatch({
  Phi_nu_alpha_ <- crossprod(cbind(h_nu_mat_) * r_ps2_nu_,
                              t(t(dot_pi_full) * rep(nu_hat_, alpha_dim_))) / n_val
  cat(sprintf("Phi_nu_alpha dim: %s\n", paste(dim(Phi_nu_alpha_), collapse="x")))
}, error = function(e) cat(sprintf("Error at Phi_nu_alpha: %s\n", conditionMessage(e))))

tryCatch({
  GtW_nu_ <- crossprod(Gamma_nu_, W_nu_)
  cat(sprintf("GtW_nu dim: %s\n", paste(dim(GtW_nu_), collapse="x")))
  dot_nu_ <- -solve(GtW_nu_ %*% Gamma_nu_) %*% GtW_nu_ %*% Phi_nu_alpha_
  cat(sprintf("dot_nu dim: %s\n", paste(dim(dot_nu_), collapse="x")))
}, error = function(e) cat(sprintf("Error at dot_nu: %s\n", conditionMessage(e))))

tryCatch({
  adj1 <- (t(wH_nu_) %*% dot_nu_) %*% psi2_alpha
  cat(sprintf("adj1 dim: %s\n", paste(dim(adj1), collapse="x")))
}, error = function(e) cat(sprintf("Error at adj1: %s\n", conditionMessage(e))))

tryCatch({
  adj2 <- t(wH_nu_) %*% psi_nu_
  cat(sprintf("adj2 dim: %s\n", paste(dim(adj2), collapse="x")))
}, error = function(e) cat(sprintf("Error at adj2: %s\n", conditionMessage(e))))
