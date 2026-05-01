## Idea 2: Joint selection-outcome model
## Stacked moments:
##   g1 = (R/pi(X,Y;alpha) - 1) * h(X)           [IPW, dim=5]
##   g2 = R * (Y - m(X;alpha,beta)) * X_out       [outcome w/ offset, dim=5]
## where m(X;alpha,beta) = expit(log(pi1/pi0) + X_out'beta)
## Total: 10 moments, 9 parameters (alpha:4, beta:5), overid=1
##
## WHY THIS BREAKS DEGENERACY:
## At degenerate alpha (alpha_y -> -inf), pi1 -> 1, so
## log(pi1/pi0) -> -log(pi0) -> +large. Then m(X) -> expit(+inf) -> 1
## for all X. But Y=0 respondents exist, so g2 = R*(0 - 1)*X ≠ 0.
## No finite beta can fix this because the offset overwhelms X'beta.
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

n_val <- 2000; n_reps <- 1000
ps_spec <- get_ps_spec("9-alt1")
W_fn_raw <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
compute_pi_fn <- function(eta) plogis(-eta)
tol <- 1e-6; max_outer <- 100L; cond_limit <- 1e8

run_joint_gmm <- function(dat, design_mat, h_x, link_type) {
  r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
  n <- nrow(dat); k_alpha <- ncol(design_mat); h_dim <- ncol(h_x)
  X_out <- cbind(1, dat$u1, dat$u2, dat$z1, dat$z2)  # outcome model covariates
  k_beta <- ncol(X_out)
  k_total <- k_alpha + k_beta  # 4 + 5 = 9
  g_total_dim <- h_dim + k_beta  # 5 + 5 = 10

  # Compute pi1(X;alpha) and pi0(X;alpha) for all obs
  compute_pi1_pi0 <- function(alpha) {
    # design_mat has columns: intercept, y, covariates
    # For pi1: set y=1; for pi0: set y=0
    # Find which column is y
    dm1 <- design_mat; dm0 <- design_mat
    dm1[, 2] <- 1; dm0[, 2] <- 0  # column 2 is y
    pi1 <- compute_pi_fn(as.vector(dm1 %*% alpha))
    pi0 <- compute_pi_fn(as.vector(dm0 %*% alpha))
    list(pi1 = pi1, pi0 = pi0)
  }

  # Stacked moment function g(theta) = [g1; g2]
  g_stacked <- function(theta) {
    alpha <- theta[1:k_alpha]
    beta <- theta[(k_alpha+1):k_total]
    pi_v <- compute_pi_fn(as.vector(design_mat %*% alpha))
    g1 <- (r_vec / pi_v - 1) * h_x  # n x h_dim

    pp <- compute_pi1_pi0(alpha)
    # Offset: log(pi1/pi0), clamp to avoid Inf
    log_ratio <- log(pmax(pp$pi1, 1e-15)) - log(pmax(pp$pi0, 1e-15))
    m_val <- plogis(log_ratio + as.vector(X_out %*% beta))
    g2 <- r_vec * (y_vec - m_val) * X_out  # n x k_beta

    cbind(g1, g2)  # n x g_total_dim
  }

  # Gradient of GMM objective w.r.t. theta (numerical)
  gmm_grad_num <- function(theta, W_hat) {
    G0 <- colMeans(g_stacked(theta))
    grad <- numeric(k_total)
    eps <- 1e-5
    for (j in 1:k_total) {
      theta_p <- theta; theta_p[j] <- theta_p[j] + eps
      theta_m <- theta; theta_m[j] <- theta_m[j] - eps
      Gp <- colMeans(g_stacked(theta_p))
      Gm <- colMeans(g_stacked(theta_m))
      Gamma_j <- (Gp - Gm) / (2 * eps)
      grad[j] <- 2 * sum(Gamma_j * (W_hat %*% G0))
    }
    grad
  }

  # Objective function
  gmm_obj <- function(theta, W_hat) {
    G <- colMeans(g_stacked(theta))
    as.numeric(crossprod(G, W_hat %*% G))
  }

  # Initialize: alpha=0, beta from logistic regression of Y on X among R=1
  resp <- r_vec == 1
  beta_init <- tryCatch({
    fit <- glm(y_vec[resp] ~ X_out[resp, -1], family = binomial)
    as.vector(coef(fit))
  }, error = function(e) rep(0, k_beta))
  if (any(is.na(beta_init))) beta_init[is.na(beta_init)] <- 0

  theta <- c(rep(0, k_alpha), beta_init)

  # Step 1: W = I
  W_hat <- diag(g_total_dim)
  opt <- optim(theta, function(th) gmm_obj(th, W_hat),
               gr = function(th) gmm_grad_num(th, W_hat),
               method = "L-BFGS-B", control = list(maxit = 500))
  theta <- opt$par

  # Step 2: Iterative GMM
  for (t in 1:max_outer) {
    g_mat <- g_stacked(theta)
    G_vec <- colMeans(g_mat)
    W_hat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(g_total_dim))

    # Check gradient norm (only for alpha part)
    grad_full <- gmm_grad_num(theta, W_hat)
    grad_norm <- max(abs(grad_full[1:k_alpha]))
    if (grad_norm < tol) break

    opt <- optim(theta, function(th) gmm_obj(th, W_hat),
                 gr = function(th) gmm_grad_num(th, W_hat),
                 method = "L-BFGS-B", control = list(maxit = 500))
    theta <- opt$par
  }

  alpha_hat <- theta[1:k_alpha]
  beta_hat <- theta[(k_alpha+1):k_total]
  pi_final <- compute_pi_fn(as.vector(design_mat %*% alpha_hat))

  # Forward cond (using only g1 moments)
  g1_mat <- (r_vec / pi_final - 1) * h_x
  W1 <- tryCatch(solve(crossprod(g1_mat) / n), error = function(e) diag(h_dim))
  cf <- r_vec * (1 - pi_final) / pi_final
  Gamma1 <- crossprod(h_x * cf, design_mat) / n
  H1 <- crossprod(Gamma1, W1 %*% Gamma1)
  eig1 <- eigen(H1, symmetric = TRUE, only.values = TRUE)$values
  fwd_cond <- max(eig1) / max(min(eig1), 1e-15)

  list(alpha = alpha_hat, beta = beta_hat, cond = fwd_cond, pi_hat = pi_final,
       grad_norm = max(abs(grad_full[1:k_alpha])))
}

# SE: use sandwich formula with full stacked moments
compute_se_joint <- function(dat, alpha_hat, beta_hat, design_mat, h_x, link_type, n) {
  r_vec <- dat[["r"]]; y_vec <- dat[["y"]]
  pi_hat <- compute_pi_fn(as.vector(design_mat %*% alpha_hat))
  mu_hat <- mean(r_vec * y_vec / pi_hat)
  k_alpha <- length(alpha_hat); k_beta <- length(beta_hat)
  k_total <- k_alpha + k_beta; h_dim <- ncol(h_x)
  X_out <- cbind(1, dat$u1, dat$u2, dat$z1, dat$z2)
  g_total_dim <- h_dim + k_beta
  theta <- c(alpha_hat, beta_hat)

  # Compute stacked moments
  dm1 <- design_mat; dm0 <- design_mat; dm1[,2] <- 1; dm0[,2] <- 0
  pi1 <- compute_pi_fn(as.vector(dm1 %*% alpha_hat))
  pi0 <- compute_pi_fn(as.vector(dm0 %*% alpha_hat))
  log_ratio <- log(pmax(pi1, 1e-15)) - log(pmax(pi0, 1e-15))
  m_val <- plogis(log_ratio + as.vector(X_out %*% beta_hat))
  g1 <- (r_vec / pi_hat - 1) * h_x
  g2 <- r_vec * (y_vec - m_val) * X_out
  g_mat <- cbind(g1, g2)

  W_hat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) NULL)
  if (is.null(W_hat)) return(list(mu = mu_hat, se = NA))

  # Numerical Jacobian Gamma = E[dg/dtheta']
  eps <- 1e-5
  Gamma <- matrix(0, g_total_dim, k_total)
  g_fn_full <- function(th) {
    a <- th[1:k_alpha]; b <- th[(k_alpha+1):k_total]
    pv <- compute_pi_fn(as.vector(design_mat %*% a))
    dm1t <- design_mat; dm0t <- design_mat; dm1t[,2] <- 1; dm0t[,2] <- 0
    p1 <- compute_pi_fn(as.vector(dm1t %*% a)); p0 <- compute_pi_fn(as.vector(dm0t %*% a))
    lr <- log(pmax(p1,1e-15)) - log(pmax(p0,1e-15))
    mv <- plogis(lr + as.vector(X_out %*% b))
    cbind((r_vec/pv - 1)*h_x, r_vec*(y_vec - mv)*X_out)
  }
  for (j in 1:k_total) {
    th_p <- theta; th_p[j] <- th_p[j] + eps
    th_m <- theta; th_m[j] <- th_m[j] - eps
    Gamma[, j] <- (colMeans(g_fn_full(th_p)) - colMeans(g_fn_full(th_m))) / (2*eps)
  }

  H_mat <- crossprod(Gamma, W_hat %*% Gamma)
  H_cond <- tryCatch({
    ev <- eigen(H_mat, symmetric = FALSE, only.values = TRUE)$values
    max(Mod(ev)) / max(min(Mod(ev)), 1e-15)
  }, error = function(e) Inf)
  if (!is.finite(H_cond) || H_cond > 1e15) return(list(mu = mu_hat, se = NA))
  H_inv <- tryCatch(solve(H_mat), error = function(e) NULL)
  if (is.null(H_inv)) return(list(mu = mu_hat, se = NA))

  # psi_theta = -H_inv %*% Gamma' %*% W %*% g_i' for each i
  GtW <- crossprod(Gamma, W_hat)  # k_total x g_total_dim
  psi_theta <- -H_inv %*% GtW %*% t(g_mat)  # k_total x n

  # d mu / d alpha (same as before)
  if (!is.null(link_type) && link_type == "logistic_complement") {
    dot_pi <- -design_mat * (pi_hat * (1 - pi_hat))
  } else { dot_pi <- design_mat * (pi_hat * (1 - pi_hat)) }
  ry_ps_inv2 <- as.vector(r_vec * y_vec * (pi_hat^(-2)))
  H_alpha <- colMeans(dot_pi * ry_ps_inv2)  # length k_alpha

  # d mu / d theta = (d mu / d alpha, 0)
  dmu_dtheta <- c(H_alpha, rep(0, k_beta))

  # mu influence: psi_mu_i = R_i*Y_i/pi_i - dmu_dtheta' %*% psi_theta_i
  mu_iid <- as.vector(r_vec * y_vec / pi_hat) - as.vector(t(dmu_dtheta) %*% psi_theta)
  se <- sqrt(var(mu_iid) / n)
  list(mu = mu_hat, se = se)
}

settings <- list(
  list(name = "S4 M3", setting = "setting4", m_idx = 3),
  list(name = "S4 M2", setting = "setting4", m_idx = 2),
  list(name = "S3 M2", setting = "setting3", m_idx = 2)
)
cat("=== Idea 2: Joint selection-outcome model (stacked moments) ===\n\n")

for (cfg in settings) {
  cat(sprintf("========== %s ==========\n", cfg$name))
  data_file <- misspecified_model_all_data_file.list[[cfg$setting]][["miss50"]][[1]]
  all_data <- readRDS(data_file)
  mu_true <- get_mu_true(cfg$setting)

  cat("Pre-building design matrices...\n")
  design_info <- vector("list", n_reps)
  for (rep_i in 1:n_reps) {
    dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
    single_ps <- list(
      formula.list = list(ps_spec[["formula.list"]][[cfg$m_idx]]),
      h_alpha.list = list(ps_spec[["h_alpha.list"]][[cfg$m_idx]]),
      inv_link = ps_spec[["inv_link"]], outcome = ps_spec[["outcome"]],
      alpha_init.list = list(NULL), optimizer = "L-BFGS-B")
    tryCatch({
      ebmr <- EBMRAlgorithmFast4$new("y", single_ps, dat, W_fn_raw)
      design_info[[rep_i]] <- list(
        design_mat = ebmr$ps_fit.list[[1]]$design_matrix,
        h_x = ebmr$ps_fit.list[[1]]$h_x,
        link_type = ebmr$ps_fit.list[[1]]$link_type)
    }, error = function(e) NULL)
  }
  cat("Done.\n")

  mu_vals <- se_vals <- cond_vals <- grad_vals <- rep(NA_real_, n_reps)
  for (rep_i in 1:n_reps) {
    if (rep_i %% 200 == 0) cat(sprintf("  rep %d...\n", rep_i))
    if (is.null(design_info[[rep_i]])) next
    dat <- all_data[((rep_i-1)*n_val + 1):(rep_i*n_val), ]
    di <- design_info[[rep_i]]
    tryCatch({
      res <- run_joint_gmm(dat, di$design_mat, di$h_x, di$link_type)
      cond_vals[rep_i] <- res$cond; grad_vals[rep_i] <- res$grad_norm
      se_res <- compute_se_joint(dat, res$alpha, res$beta, di$design_mat, di$h_x, di$link_type, n_val)
      mu_vals[rep_i] <- se_res$mu; se_vals[rep_i] <- se_res$se
    }, error = function(e) cat(sprintf("  Rep %d ERROR: %s\n", rep_i, e$message)))
  }
  valid <- !is.na(mu_vals) & !is.na(se_vals)
  n_degen <- sum(cond_vals[valid] >= cond_limit, na.rm = TRUE)
  n_conv <- sum(grad_vals[valid] < tol, na.rm = TRUE)
  cat(sprintf("  Valid: %d/%d, Degen: %d, Conv: %d\n", sum(valid), n_reps, n_degen, n_conv))
  cat(sprintf("  Bias     = %.4f\n", mean(mu_vals[valid]) - mu_true))
  cat(sprintf("  ESD      = %.4f\n", sd(mu_vals[valid])))
  cat(sprintf("  ESE      = %.4f\n", mean(se_vals[valid])))
  cat(sprintf("  ESE/ESD  = %.3f\n", mean(se_vals[valid]) / sd(mu_vals[valid])))
  ci_lo <- mu_vals[valid] - 1.96 * se_vals[valid]
  ci_hi <- mu_vals[valid] + 1.96 * se_vals[valid]
  cat(sprintf("  CP       = %.3f\n", mean(mu_true >= ci_lo & mu_true <= ci_hi)))
  cat("\n")
}
