## Fix: Multi-start the inner loop of nu GMM (within the package's iterative scheme)
## For each W update, try multiple starting points and pick the one with lowest objective
## Use the Fast4 package for alpha estimation, custom nu solver
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

## Multi-start inner-loop nu GMM
## At each W iteration, solve from multiple starts, keep the best
solve_nu_multistart_inner <- function(ps_mat, h_x, r_vec, nn, max_outer = 200) {
  J <- ncol(ps_mat); h_dim <- ncol(h_x)

  # Starting points for inner loop
  inner_starts <- list(
    rep(1/J, J),
    c(1, rep(0.001, J-1))
  )
  for (k in 2:J) { s <- rep(0.001, J); s[k] <- 1; inner_starts[[k+1]] <- s }

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

  solve_inner <- function(init, W) {
    tryCatch({
      res <- optim(init, function(x) obj_fixed(x, W), function(x) grad_fixed(x, W),
                   method = "L-BFGS-B", control = list(maxit = 1000))
      list(par = res$par, value = res$value)
    }, error = function(e) list(par = init, value = Inf))
  }

  # Step 1: W = I, multi-start
  W_hat <- diag(h_dim)
  best_val <- Inf; best_est <- rep(1/J, J)
  for (init in inner_starts) {
    r <- solve_inner(init, W_hat)
    if (r$value < best_val) { best_val <- r$value; best_est <- r$par }
  }
  estimates <- best_est

  # Step 2: Iterate W with multi-start inner
  for (t in 1:max_outer) {
    g_mat <- (r_vec / as.vector(ps_mat %*% estimates) - 1) * h_x
    W_hat <- tryCatch(W_fn(g_mat), error = function(e) diag(h_dim))
    gn <- max(abs(grad_fixed(estimates, W_hat)))
    if (gn < 1e-6) break

    # Multi-start: try from current + all corner starts
    all_starts <- c(list(estimates), inner_starts)
    best_val <- Inf; best_est <- estimates
    for (init in all_starts) {
      r <- solve_inner(init, W_hat)
      if (r$value < best_val) { best_val <- r$value; best_est <- r$par }
    }
    estimates <- best_est
  }

  w <- estimates^2 / sum(estimates^2)
  list(nu = estimates, w = w)
}

# Analytical SE
ebmr_pt <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat, W_fn)
res_pt <- ebmr_pt$EBMR_IPW(h_nu = h_nu_fn, true_ps = NULL)
anal_se <- res_pt$se_ipw
cat(sprintf("Analytical SE: %.4f\n\n", anal_se))

B <- 300
cat(sprintf("=== Multi-start inner loop (B=%d) ===\n\n", B))

mus_fix <- mus_default <- rep(NA, B)

for (b in 1:B) {
  if (b %% 50 == 0) cat(sprintf("  rep %d...\n", b))
  set.seed(12345 + b)
  idx <- sample(1:n, n, replace = TRUE)
  dat_b <- dat[idx, ]

  tryCatch({
    ebmr_b <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat_b, W_fn)

    # Default
    res_d <- ebmr_b$EBMR_IPW(h_nu = h_nu_fn, true_ps = NULL, se.fit = FALSE)
    mus_default[b] <- res_d$mu_ipw

    # Multi-start inner
    ps_mat <- do.call(cbind, lapply(ebmr_b$ps_fit.list, function(pf) pf$fitted.values))
    h_x <- cbind(1, h_nu_fn(dat_b))
    n_b <- nrow(dat_b)
    nu_res <- solve_nu_multistart_inner(ps_mat, h_x, dat_b$r, n_b)
    ps_nu <- as.vector(ps_mat %*% nu_res$nu)
    mus_fix[b] <- mean(dat_b$r * dat_b$teacher_report / ps_nu)
  }, error = function(e) {})
}

cat(sprintf("\nDefault:          boot_se=%.4f, ratio=%.3f\n",
    sd(mus_default, na.rm=TRUE), anal_se / sd(mus_default, na.rm=TRUE)))
cat(sprintf("Multistart inner: boot_se=%.4f, ratio=%.3f\n",
    sd(mus_fix, na.rm=TRUE), anal_se / sd(mus_fix, na.rm=TRUE)))
cat(sprintf("Agree: %d/%d (%.1f%%)\n",
    sum(abs(mus_default - mus_fix) < 0.001, na.rm=TRUE),
    sum(!is.na(mus_default) & !is.na(mus_fix)),
    100*mean(abs(mus_default - mus_fix) < 0.001, na.rm=TRUE)))
