## Approach: Track CUE objective across all outer iterations, return the best nu seen
## The iterative GMM may visit the M1-dominated solution early (near W=I)
## but drift away. Capture it before it drifts.
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

## Nu GMM that tracks best CUE objective across iterations
solve_nu_best_cue <- function(ps_mat, h_x, r_vec, nn, max_outer = 200) {
  J <- ncol(ps_mat); h_dim <- ncol(h_x)

  cue_obj <- function(nu) {
    ps_nu <- as.vector(ps_mat %*% nu)
    g_mat <- as.vector(r_vec / ps_nu - 1) * h_x
    G <- colMeans(g_mat)
    W <- tryCatch(solve(crossprod(g_mat) / nn), error = function(e) diag(h_dim))
    as.numeric(t(G) %*% W %*% G)
  }

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

  # Track best CUE
  best_cue <- cue_obj(estimates)
  best_nu <- estimates

  # Iterate
  for (t in 1:max_outer) {
    g_mat <- (r_vec / as.vector(ps_mat %*% estimates) - 1) * h_x
    W_hat <- tryCatch(W_fn(g_mat), error = function(e) diag(h_dim))
    gn <- max(abs(grad_fixed(estimates, W_hat)))
    if (gn < 1e-6) break

    res <- optim(estimates, function(x) obj_fixed(x, W_hat), function(x) grad_fixed(x, W_hat),
                 method = "L-BFGS-B", control = list(maxit = 1000))
    estimates <- res$par

    # Check CUE at current iterate
    cur_cue <- tryCatch(cue_obj(estimates), error = function(e) Inf)
    if (cur_cue < best_cue) {
      best_cue <- cur_cue
      best_nu <- estimates
    }
  }

  w <- best_nu^2 / sum(best_nu^2)
  list(nu = best_nu, w = w, cue = best_cue)
}

anal_se <- ebmr_pt <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat, W_fn)
res_pt <- ebmr_pt$EBMR_IPW(h_nu = h_nu_fn, true_ps = NULL)
anal_se <- res_pt$se_ipw

B <- 300
cat(sprintf("=== Best-CUE-across-iterations (B=%d) ===\n\n", B))

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

    # Best CUE across iterations
    ps_mat <- do.call(cbind, lapply(ebmr_b$ps_fit.list, function(pf) pf$fitted.values))
    h_x <- cbind(1, h_nu_fn(dat_b))
    n_b <- nrow(dat_b)
    nu_res <- solve_nu_best_cue(ps_mat, h_x, dat_b$r, n_b)
    ps_nu <- as.vector(ps_mat %*% nu_res$nu)
    mus_fix[b] <- mean(dat_b$r * dat_b$teacher_report / ps_nu)
  }, error = function(e) {})
}

cat(sprintf("\nDefault:          boot_se=%.4f, ratio=%.3f\n",
    sd(mus_default, na.rm=TRUE), anal_se / sd(mus_default, na.rm=TRUE)))
cat(sprintf("Best-CUE-iterate: boot_se=%.4f, ratio=%.3f\n",
    sd(mus_fix, na.rm=TRUE), anal_se / sd(mus_fix, na.rm=TRUE)))
