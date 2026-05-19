## Test: Damped W update for nu estimation
## W_new = (1-lambda)*W_old + lambda*W_current
## This slows W adaptation, preventing drift to degenerate solution
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Mental_Health_Data_Application")
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)
library(Matrix); library(numDeriv)
source("MHD_functions.R")

original_data <- read.csv("data_application.csv")
percent <- original_data$Percentage
n <- 2486
class_n <- round(n * percent / 100)
dat <- gen_data(original_data, class_n, n)
dat$y <- dat$teacher_report

W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_alpha <- c("health", "father", "parent_report", "fp", "fh", "hp", "fhp")
h_nu_fn <- function(data) cbind(health=data$health, father=data$father, parent_report=data$parent_report,
                                  fp=data$fp, fh=data$fh, hp=data$hp, fhp=data$fhp)

ps_spec <- list(
  formula.list = list(
    r ~ teacher_report + health + father + teacher_report:health + teacher_report:father,
    r ~ teacher_report + health + parent_report + teacher_report:health + teacher_report:parent_report,
    r ~ teacher_report + parent_report + father + teacher_report:parent_report + teacher_report:father
  ),
  h_alpha.list = list(h_alpha, h_alpha, h_alpha),
  outcome = "teacher_report", inv_link = function(eta) 1/(1+exp(eta)), optimizer = "L-BFGS-B"
)

## Damped W nu solver
solve_nu_damped <- function(ps_mat, h_x, r_vec, nn, lambda = 0.5, max_outer = 200) {
  J <- ncol(ps_mat); h_dim <- ncol(h_x)

  obj_fixed <- function(nu, W) {
    ps_nu <- as.vector(ps_mat %*% nu)
    g_mat <- as.vector(r_vec / ps_nu - 1) * h_x
    G <- colMeans(g_mat)
    as.numeric(t(G) %*% W %*% G)
  }
  grad_fixed <- function(nu, W) {
    ps_nu <- as.vector(ps_mat %*% nu)
    g_mat <- as.vector(r_vec / ps_nu - 1) * h_x
    G <- matrix(colMeans(g_mat), h_dim, 1)
    Gamma <- crossprod(h_x * (-r_vec/(ps_nu^2)), ps_mat) / nn
    2 * as.vector(crossprod(Gamma, W %*% G))
  }

  # Step 1: W = I
  W_hat <- diag(h_dim)
  res <- optim(rep(1/J, J), function(x) obj_fixed(x, W_hat), function(x) grad_fixed(x, W_hat),
               method = "L-BFGS-B", control = list(maxit = 1000))
  estimates <- res$par

  # Iterate with damped W
  for (t in 1:max_outer) {
    g_mat <- (r_vec / as.vector(ps_mat %*% estimates) - 1) * h_x
    W_new <- tryCatch(W_fn(g_mat), error = function(e) diag(h_dim))
    # Damped update
    W_hat <- (1 - lambda) * W_hat + lambda * W_new

    gn <- max(abs(grad_fixed(estimates, W_hat)))
    if (gn < 1e-6) break

    res <- optim(estimates, function(x) obj_fixed(x, W_hat), function(x) grad_fixed(x, W_hat),
                 method = "L-BFGS-B", control = list(maxit = 1000))
    estimates <- res$par
  }

  w <- estimates^2 / sum(estimates^2)
  list(nu = estimates, w = w)
}

anal_se <- 0.0333  # from package

B <- 300
lambdas <- c(0.1, 0.2, 0.3, 0.5, 1.0)

cat("=== Damped W update for nu (B=300) ===\n\n")

for (lam in lambdas) {
  mus <- rep(NA, B)
  for (b in 1:B) {
    set.seed(12345 + b)
    idx <- sample(1:n, n, replace = TRUE)
    dat_b <- dat[idx, ]
    tryCatch({
      ebmr_b <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat_b, W_fn)
      ps_mat <- do.call(cbind, lapply(ebmr_b$ps_fit.list, function(pf) pf$fitted.values))
      h_x <- cbind(1, h_nu_fn(dat_b))
      n_b <- nrow(dat_b)
      nu_res <- solve_nu_damped(ps_mat, h_x, dat_b$r, n_b, lambda = lam)
      ps_nu <- as.vector(ps_mat %*% nu_res$nu)
      mus[b] <- mean(dat_b$r * dat_b$teacher_report / ps_nu)
    }, error = function(e) {})
  }
  boot_se <- sd(mus, na.rm = TRUE)
  cat(sprintf("  lambda=%.1f: boot_se=%.4f, ratio=%.3f\n", lam, boot_se, anal_se / boot_se))
}
