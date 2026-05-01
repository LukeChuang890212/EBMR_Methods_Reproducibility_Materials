## Check: Is (1,0,0) the global minimum of the CUE objective on original data?
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
h_nu <- function(data) cbind(health=data$health, father=data$father, parent_report=data$parent_report,
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

ebmr <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat, W_fn)
ps_mat <- do.call(cbind, lapply(ebmr$ps_fit.list, function(pf) pf$fitted.values))
h_x <- cbind(1, h_nu(dat))
r_vec <- dat$r

## CUE objective: G(nu)' W(nu) G(nu)
cue_obj <- function(nu) {
  ps_nu <- as.vector(ps_mat %*% nu)
  g_mat <- as.vector(r_vec / ps_nu - 1) * h_x
  G <- colMeans(g_mat)
  W_hat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(ncol(h_x)))
  as.numeric(t(G) %*% W_hat %*% G)
}

## W=I objective: G(nu)' G(nu)
wi_obj <- function(nu) {
  ps_nu <- as.vector(ps_mat %*% nu)
  g_mat <- as.vector(r_vec / ps_nu - 1) * h_x
  G <- colMeans(g_mat)
  sum(G^2)
}

cat("=== CUE and W=I objectives on original data ===\n\n")

# Package converged solution
res <- ebmr$EBMR_IPW(h_nu = h_nu, true_ps = NULL, se.fit = FALSE)
cat(sprintf("Package converged: nu=(%s), w=(%s)\n",
    paste(round(res$nu.hat, 4), collapse=","), paste(round(res$w.hat, 4), collapse=",")))
cat(sprintf("  CUE obj = %.6e\n", cue_obj(res$nu.hat)))
cat(sprintf("  W=I obj = %.6e\n\n", wi_obj(res$nu.hat)))

# Grid search over nu space
cat("Grid search:\n")
cat(sprintf("%-30s %12s %12s\n", "nu (approx w)", "CUE obj", "W=I obj"))
cat(paste(rep("-", 56), collapse=""), "\n")

test_points <- list(
  list(nu = c(1, 0.001, 0.001), desc = "M1 only"),
  list(nu = c(0.001, 1, 0.001), desc = "M2 only"),
  list(nu = c(0.001, 0.001, 1), desc = "M3 only"),
  list(nu = c(1/3, 1/3, 1/3), desc = "Equal"),
  list(nu = c(0.7, 0.3, 0.001), desc = "M1+M2"),
  list(nu = c(0.7, 0.001, 0.3), desc = "M1+M3"),
  list(nu = c(0.001, 0.7, 0.3), desc = "M2+M3"),
  list(nu = c(0.5, 0.3, 0.2), desc = "M1 heavy"),
  list(nu = c(0.9, 0.05, 0.05), desc = "M1 very heavy"),
  list(nu = res$nu.hat, desc = "Converged")
)

for (tp in test_points) {
  w <- tp$nu^2 / sum(tp$nu^2)
  cue <- tryCatch(cue_obj(tp$nu), error = function(e) NA)
  wi <- tryCatch(wi_obj(tp$nu), error = function(e) NA)
  cat(sprintf("%-30s %12.6e %12.6e\n",
      sprintf("%s w=(%s)", tp$desc, paste(round(w,3), collapse=",")), cue, wi))
}

# Fine grid around M1-dominated region
cat("\nFine grid around M1-dominated:\n")
for (nu2 in c(0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3)) {
  for (nu3 in c(0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3)) {
    nu1 <- 1
    nu <- c(nu1, nu2, nu3)
    w <- nu^2 / sum(nu^2)
    cue <- tryCatch(cue_obj(nu), error = function(e) NA)
    if (!is.na(cue)) {
      cat(sprintf("  nu=(1, %.3f, %.3f) w=(%.4f,%.4f,%.4f) CUE=%.6e\n",
          nu2, nu3, w[1], w[2], w[3], cue))
    }
  }
}

# Run nloptr from many starting points to find true global min
cat("\nGlobal search (100 random starts):\n")
set.seed(42)
best_cue <- Inf; best_nu <- NULL
for (i in 1:100) {
  init <- abs(rnorm(3))
  tryCatch({
    # CUE iterative
    est <- init / sum(abs(init))
    for (t in 1:50) {
      ps_nu <- as.vector(ps_mat %*% est)
      g_mat <- as.vector(r_vec / ps_nu - 1) * h_x
      W_hat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(ncol(h_x)))

      obj_W <- function(nu) {
        ps_nu <- as.vector(ps_mat %*% nu)
        g <- as.vector(r_vec / ps_nu - 1) * h_x
        G <- colMeans(g)
        as.numeric(t(G) %*% W_hat %*% G)
      }
      grad_W <- function(nu) {
        ps_nu <- as.vector(ps_mat %*% nu)
        g <- as.vector(r_vec / ps_nu - 1) * h_x
        G <- matrix(colMeans(g), ncol(h_x), 1)
        Gamma <- crossprod(h_x * (-r_vec/(ps_nu^2)), ps_mat) / n
        2 * as.vector(crossprod(Gamma, W_hat %*% G))
      }

      gn <- max(abs(grad_W(est)))
      if (gn < 1e-8) break
      res_i <- nloptr(x0=est, eval_f=obj_W, eval_grad_f=grad_W,
                      opts=list(algorithm="NLOPT_LD_LBFGS", maxeval=5000, xtol_rel=1e-14))
      est <- res_i$solution
    }

    cue_val <- cue_obj(est)
    if (cue_val < best_cue) {
      best_cue <- cue_val
      best_nu <- est
      w <- est^2 / sum(est^2)
      cat(sprintf("  Start %3d: NEW BEST CUE=%.6e, w=(%s)\n", i,
          cue_val, paste(round(w, 4), collapse=",")))
    }
  }, error = function(e) {})
}

cat(sprintf("\nGlobal best: CUE=%.6e, nu=(%s), w=(%s)\n",
    best_cue, paste(round(best_nu, 4), collapse=","),
    paste(round(best_nu^2/sum(best_nu^2), 4), collapse=",")))
