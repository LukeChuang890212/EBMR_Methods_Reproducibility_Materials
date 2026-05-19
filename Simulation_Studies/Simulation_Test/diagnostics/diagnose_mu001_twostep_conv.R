## Check if alpha converges in two-step GMM for mu_001
## i.e., does step 2 alpha differ much from step 1 alpha?
## And what happens if we do a 3rd step?
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

ps_spec_base <- get_ps_spec("9-alt1")
h_alpha_fn <- function(dat) cbind(
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

cat("=== Multi-step GMM: tracking alpha across steps ===\n\n")

for (rep_i in 1:5) {
  dat <- all_data[((rep_i-1)*2000 + 1):(rep_i*2000), ]
  design_mat <- model.matrix(ps_spec_base$formula.list[[3]], dat)
  r_vec <- as.vector(dat$r)
  n <- 2000
  h_x <- cbind(1, h_alpha_fn(dat))
  h_dim <- ncol(h_x)
  alpha_dim <- ncol(design_mat)

  compute_pi <- function(eta) 1/(1+exp(eta))

  g_fn <- function(alpha) {
    eta <- as.vector(design_mat %*% alpha)
    eta <- pmin(pmax(eta, -20), 20)
    pi_vec <- compute_pi(eta)
    as.vector(r_vec/pi_vec - 1) * h_x
  }
  G_fn <- function(alpha) colMeans(g_fn(alpha))

  cat(sprintf("--- Rep %d ---\n", rep_i))

  # Step 1: W = I
  obj1 <- function(a) { G <- G_fn(a); sum(G^2) }
  opt1 <- optim(rep(0, alpha_dim), obj1, method = "L-BFGS-B",
                control = list(maxit = 1000))
  alpha1 <- opt1$par
  G1 <- G_fn(alpha1)
  cat(sprintf("  Step 1 (W=I):     alpha=(%s), max|G|=%.4e\n",
      paste(round(alpha1, 4), collapse=", "), max(abs(G1))))

  # Step 2: W = solve(var(g(alpha1)))
  g1 <- g_fn(alpha1)
  W2 <- solve(var(g1))
  obj2 <- function(a) { G <- matrix(G_fn(a), ncol=1); as.numeric(t(G) %*% W2 %*% G) }
  opt2 <- optim(alpha1, obj2, method = "L-BFGS-B", control = list(maxit = 1000))
  alpha2 <- opt2$par
  G2 <- G_fn(alpha2)
  cat(sprintf("  Step 2 (W=var1):  alpha=(%s), max|G|=%.4e, |alpha2-alpha1|=%.4f\n",
      paste(round(alpha2, 4), collapse=", "), max(abs(G2)), max(abs(alpha2-alpha1))))

  # Step 3: W = solve(var(g(alpha2)))
  g2 <- g_fn(alpha2)
  W3 <- solve(var(g2))
  obj3 <- function(a) { G <- matrix(G_fn(a), ncol=1); as.numeric(t(G) %*% W3 %*% G) }
  opt3 <- optim(alpha2, obj3, method = "L-BFGS-B", control = list(maxit = 1000))
  alpha3 <- opt3$par
  G3 <- G_fn(alpha3)
  cat(sprintf("  Step 3 (W=var2):  alpha=(%s), max|G|=%.4e, |alpha3-alpha2|=%.4f\n",
      paste(round(alpha3, 4), collapse=", "), max(abs(G3)), max(abs(alpha3-alpha2))))

  # Step 4: W = solve(var(g(alpha3)))
  g3 <- g_fn(alpha3)
  W4 <- solve(var(g3))
  obj4 <- function(a) { G <- matrix(G_fn(a), ncol=1); as.numeric(t(G) %*% W4 %*% G) }
  opt4 <- optim(alpha3, obj4, method = "L-BFGS-B", control = list(maxit = 1000))
  alpha4 <- opt4$par
  G4 <- G_fn(alpha4)
  cat(sprintf("  Step 4 (W=var3):  alpha=(%s), max|G|=%.4e, |alpha4-alpha3|=%.4f\n",
      paste(round(alpha4, 4), collapse=", "), max(abs(G4)), max(abs(alpha4-alpha3))))

  # Step 5
  g4 <- g_fn(alpha4)
  W5 <- solve(var(g4))
  obj5 <- function(a) { G <- matrix(G_fn(a), ncol=1); as.numeric(t(G) %*% W5 %*% G) }
  opt5 <- optim(alpha4, obj5, method = "L-BFGS-B", control = list(maxit = 1000))
  alpha5 <- opt5$par
  G5 <- G_fn(alpha5)
  cat(sprintf("  Step 5 (W=var4):  alpha=(%s), max|G|=%.4e, |alpha5-alpha4|=%.4f\n\n",
      paste(round(alpha5, 4), collapse=", "), max(abs(G5)), max(abs(alpha5-alpha4))))
}
