## Test: Fix inner optimizer for J=3 nu
## Strategy A: Start from M1-dominated point (1, 0.001, 0.001) instead of (1/3, 1/3, 1/3)
## Strategy B: Global optimizer (NLOPT_GN_DIRECT)
## Strategy C: Start from each corner, pick best
## Strategy D: Constrain nu1 >= 0.5 (force M1 dominance)
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Mental_Health_Data_Application")
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)
library(Matrix); library(numDeriv); library(nloptr)
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

## Custom nu solver with W=I (since W updates don't matter)
solve_nu <- function(ps_mat, h_x, r_vec, nn, init, algo = "NLOPT_LD_LBFGS",
                     lower = NULL, upper = NULL) {
  J <- ncol(ps_mat); h_dim <- ncol(h_x)
  W_hat <- diag(h_dim)
  obj_f <- function(nu) {
    ps_nu <- as.vector(ps_mat %*% nu)
    g_mat <- as.vector(r_vec / ps_nu - 1) * h_x
    G <- colMeans(g_mat)
    as.numeric(t(G) %*% G)  # W=I
  }
  grad_f <- function(nu) {
    ps_nu <- as.vector(ps_mat %*% nu)
    g_mat <- as.vector(r_vec / ps_nu - 1) * h_x
    G <- matrix(colMeans(g_mat), h_dim, 1)
    Gamma <- crossprod(h_x * (-r_vec/(ps_nu^2)), ps_mat) / nn
    2 * as.vector(crossprod(Gamma, G))
  }

  lb <- if (!is.null(lower)) lower else rep(-Inf, J)
  ub <- if (!is.null(upper)) upper else rep(Inf, J)

  res <- nloptr(x0 = init, eval_f = obj_f, eval_grad_f = grad_f,
                lb = lb, ub = ub,
                opts = list(algorithm = algo, maxeval = 5000, xtol_rel = 1e-12))
  w <- res$solution^2 / sum(res$solution^2)
  list(nu = res$solution, w = w, obj = res$objective)
}

B <- 300
cat("=== Testing inner optimizer strategies for J=3 nu (B=300) ===\n\n")

# Get analytical SE from package
ebmr_pt <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat, W_fn)
res_pt <- ebmr_pt$EBMR_IPW(h_nu = h_nu_fn, true_ps = NULL)
anal_se <- res_pt$se_ipw

strategies <- list(
  list(name = "A: M1 init (1,0.001,0.001)", init = c(1, 0.001, 0.001)),
  list(name = "B: Uniform init (1/3,1/3,1/3)", init = c(1/3, 1/3, 1/3)),
  list(name = "C: Corner multistart", init = NULL),  # special handling
  list(name = "D: Package default", init = NULL)       # use package
)

for (strat in strategies) {
  cat(sprintf("--- %s ---\n", strat$name))
  mus <- rep(NA, B)

  for (b in 1:B) {
    if (b %% 100 == 0) cat(sprintf("  rep %d...\n", b))
    set.seed(12345 + b)
    idx <- sample(1:n, n, replace = TRUE)
    dat_b <- dat[idx, ]

    tryCatch({
      ebmr_b <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat_b, W_fn)
      ps_mat <- do.call(cbind, lapply(ebmr_b$ps_fit.list, function(pf) pf$fitted.values))
      h_x <- cbind(1, h_nu_fn(dat_b))
      n_b <- nrow(dat_b)

      if (strat$name == "C: Corner multistart") {
        # Try all corners + uniform, pick lowest objective
        inits <- list(c(1/3,1/3,1/3), c(1,0.001,0.001), c(0.001,1,0.001), c(0.001,0.001,1))
        best_obj <- Inf; best_nu <- NULL
        for (ini in inits) {
          r <- solve_nu(ps_mat, h_x, dat_b$r, n_b, init = ini)
          if (r$obj < best_obj) { best_obj <- r$obj; best_nu <- r }
        }
        ps_nu <- as.vector(ps_mat %*% best_nu$nu)
        mus[b] <- mean(dat_b$r * dat_b$teacher_report / ps_nu)
      } else if (strat$name == "D: Package default") {
        res_b <- ebmr_b$EBMR_IPW(h_nu = h_nu_fn, true_ps = NULL, se.fit = FALSE)
        mus[b] <- res_b$mu_ipw
      } else {
        r <- solve_nu(ps_mat, h_x, dat_b$r, n_b, init = strat$init)
        ps_nu <- as.vector(ps_mat %*% r$nu)
        mus[b] <- mean(dat_b$r * dat_b$teacher_report / ps_nu)
      }
    }, error = function(e) {})
  }

  boot_se <- sd(mus, na.rm = TRUE)
  cat(sprintf("  boot_se=%.4f, ratio=%.3f\n\n", boot_se, anal_se / boot_se))
}
