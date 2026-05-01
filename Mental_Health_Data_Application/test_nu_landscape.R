## Deep dive: nu objective landscape for rep 17 where 110/101 pick M1 but 111 switches
## M1 is NOT degenerate in this rep
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

W <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
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

# Use rep 17 (110/101 pick M1 w1>0.99 but 111 gives w=(0.0001, 0.8911, 0.1088))
for (rep_i in c(17, 18, 35)) {
  cat(sprintf("========== Rep %d ==========\n\n", rep_i))
  set.seed(12345 + rep_i)
  idx <- sample(1:n, n, replace = TRUE)
  dat_b <- dat[idx, ]
  n_b <- nrow(dat_b)

  ebmr_b <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat_b, W)

  # Individual model diagnostics
  for (i in 1:3) {
    pf <- ebmr_b$ps_fit.list[[i]]
    ps <- pf$fitted.values
    cat(sprintf("M%d: obj=%.2e, PS=[%.4f,%.4f], >0.99:%d, <0.01:%d\n",
        i, pf$gmm_fit$opt$objective, min(ps), max(ps), sum(ps>0.99), sum(ps<0.01)))
  }

  # Run 110, 101, 111
  r110 <- ebmr_b$EBMR_IPW(h_nu = h_nu_fn, model_indices = c(1,2), true_ps = NULL, se.fit = FALSE)
  r101 <- ebmr_b$EBMR_IPW(h_nu = h_nu_fn, model_indices = c(1,3), true_ps = NULL, se.fit = FALSE)
  r111 <- ebmr_b$EBMR_IPW(h_nu = h_nu_fn, true_ps = NULL, se.fit = FALSE)

  cat(sprintf("\n110: nu=(%s), w=(%s), mu=%.4f\n",
      paste(round(r110$nu.hat, 4), collapse=","), paste(round(r110$w.hat, 4), collapse=","), r110$mu_ipw))
  cat(sprintf("101: nu=(%s), w=(%s), mu=%.4f\n",
      paste(round(r101$nu.hat, 4), collapse=","), paste(round(r101$w.hat, 4), collapse=","), r101$mu_ipw))
  cat(sprintf("111: nu=(%s), w=(%s), mu=%.4f\n\n",
      paste(round(r111$nu.hat, 4), collapse=","), paste(round(r111$w.hat, 4), collapse=","), r111$mu_ipw))

  # Compute nu objective at different weight allocations for 111
  ps_mat <- do.call(cbind, lapply(ebmr_b$ps_fit.list, function(pf) pf$fitted.values))
  h_x <- cbind(1, h_nu_fn(dat_b))
  r_b <- dat_b$r

  compute_nu_obj <- function(nu) {
    ps_nu <- as.vector(ps_mat %*% nu)
    g_mat <- as.vector(r_b / ps_nu - 1) * h_x
    G_vec <- colMeans(g_mat)
    W_hat <- tryCatch(solve(crossprod(g_mat) / n_b), error = function(e) diag(ncol(h_x)))
    as.numeric(t(G_vec) %*% W_hat %*% G_vec)
  }

  # Objective at key points
  cat("Nu objective at different weight allocations:\n")
  test_nus <- list(
    list(nu = c(1, 0.001, 0.001), desc = "M1 only"),
    list(nu = c(0.001, 1, 0.001), desc = "M2 only"),
    list(nu = c(0.001, 0.001, 1), desc = "M3 only"),
    list(nu = c(1/3, 1/3, 1/3), desc = "Equal"),
    list(nu = r111$nu.hat, desc = "Converged 111"),
    list(nu = c(r110$nu.hat[1], r110$nu.hat[2], 0.001), desc = "110 solution extended"),
    list(nu = c(r101$nu.hat[1], 0.001, r101$nu.hat[2]), desc = "101 solution extended")
  )

  for (tn in test_nus) {
    obj <- tryCatch(compute_nu_obj(tn$nu), error = function(e) NA)
    w <- tn$nu^2 / sum(tn$nu^2)
    cat(sprintf("  %-25s nu=(%s) w=(%s) obj=%.4e\n",
        tn$desc, paste(round(tn$nu, 4), collapse=","),
        paste(round(w, 4), collapse=","), obj))
  }

  # Also compute 110 and 101 objectives (J=2)
  cat("\n110 nu objective (J=2):\n")
  ps_12 <- ps_mat[, 1:2]
  h_x_12 <- h_x
  compute_nu_obj_12 <- function(nu) {
    ps_nu <- as.vector(ps_12 %*% nu)
    g_mat <- as.vector(r_b / ps_nu - 1) * h_x_12
    G_vec <- colMeans(g_mat)
    W_hat <- tryCatch(solve(crossprod(g_mat) / n_b), error = function(e) diag(ncol(h_x_12)))
    as.numeric(t(G_vec) %*% W_hat %*% G_vec)
  }
  cat(sprintf("  At converged: obj=%.4e\n", compute_nu_obj_12(r110$nu.hat)))
  cat(sprintf("  At M1 only:   obj=%.4e\n", compute_nu_obj_12(c(1, 0.001))))
  cat(sprintf("  At M2 only:   obj=%.4e\n", compute_nu_obj_12(c(0.001, 1))))

  cat("\n101 nu objective (J=2):\n")
  ps_13 <- ps_mat[, c(1,3)]
  compute_nu_obj_13 <- function(nu) {
    ps_nu <- as.vector(ps_13 %*% nu)
    g_mat <- as.vector(r_b / ps_nu - 1) * h_x
    G_vec <- colMeans(g_mat)
    W_hat <- tryCatch(solve(crossprod(g_mat) / n_b), error = function(e) diag(ncol(h_x)))
    as.numeric(t(G_vec) %*% W_hat %*% G_vec)
  }
  cat(sprintf("  At converged: obj=%.4e\n", compute_nu_obj_13(r101$nu.hat)))
  cat(sprintf("  At M1 only:   obj=%.4e\n", compute_nu_obj_13(c(1, 0.001))))
  cat(sprintf("  At M3 only:   obj=%.4e\n\n", compute_nu_obj_13(c(0.001, 1))))
}
