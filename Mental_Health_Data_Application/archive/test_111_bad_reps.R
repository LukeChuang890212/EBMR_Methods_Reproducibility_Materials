## Find bootstrap reps where M2 or M3 gets heavy weight, examine nu convergence
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

ps_specifications <- list(
  formula.list = list(
    r ~ teacher_report + health + father + teacher_report:health + teacher_report:father + father:health,
    r ~ teacher_report + parent_report + health + teacher_report:parent_report + teacher_report:health + parent_report:health,
    r ~ teacher_report + parent_report + father + teacher_report:parent_report + teacher_report:father + parent_report:father
  ),
  h_alpha.list = list(
    c("health", "father", "parent_report", "fp", "fh", "hp", "fhp"),
    c("health", "father", "parent_report", "fp", "fh", "hp", "fhp"),
    c("health", "father", "parent_report", "fp", "fh", "hp", "fhp")
  ),
  outcome = "teacher_report",
  inv_link = function(eta) 1 / (1 + exp(eta)),
  optimizer = "L-BFGS-B"
)

h_nu <- function(data) cbind(health=data$health, father=data$father, parent_report=data$parent_report,
                              fp=data$fp, fh=data$fh, hp=data$hp, fhp=data$fhp)

# First scan to find bad reps
cat("=== Scanning for non-M1-dominated reps ===\n\n")
bad_reps <- c()
for (b in 1:200) {
  set.seed(12345 + b)
  idx <- sample(1:n, n, replace = TRUE)
  dat_b <- dat[idx, ]
  tryCatch({
    ebmr_b <- EBMRAlgorithmFast4$new("teacher_report", ps_specifications, dat_b, W)
    res_b <- ebmr_b$EBMR_IPW(h_nu = h_nu, true_ps = NULL, se.fit = FALSE)
    if (res_b$w.hat[1] < 0.5) {
      bad_reps <- c(bad_reps, b)
      cat(sprintf("  Rep %d: w=(%s), mu=%.4f\n", b,
          paste(round(res_b$w.hat, 4), collapse=", "), res_b$mu_ipw))
    }
  }, error = function(e) {})
}
cat(sprintf("\nFound %d non-M1-dominated reps out of 200\n\n", length(bad_reps)))

# Deep dive into first few bad reps
cat("=== Deep dive into bad reps ===\n\n")
for (b in head(bad_reps, 5)) {
  set.seed(12345 + b)
  idx <- sample(1:n, n, replace = TRUE)
  dat_b <- dat[idx, ]
  n_b <- nrow(dat_b)

  ebmr_b <- EBMRAlgorithmFast4$new("teacher_report", ps_specifications, dat_b, W)

  # Individual model info
  for (i in 1:3) {
    pf <- ebmr_b$ps_fit.list[[i]]
    ps <- pf$fitted.values
    cat(sprintf("Rep %d, M%d: obj=%.2e, grad=%.2e, PS=[%.4f,%.4f], >0.99:%d, <0.01:%d\n",
        b, i, pf$gmm_fit$opt$objective, pf$gmm_fit$opt$final_grad_norm,
        min(ps), max(ps), sum(ps > 0.99), sum(ps < 0.01)))
  }

  # Ensemble
  res_b <- ebmr_b$EBMR_IPW(h_nu = h_nu, true_ps = NULL, se.fit = FALSE)
  cat(sprintf("Rep %d, Ensemble: mu=%.4f, nu=(%s), w=(%s)\n",
      b, res_b$mu_ipw,
      paste(round(res_b$nu.hat, 4), collapse=", "),
      paste(round(res_b$w.hat, 4), collapse=", ")))

  # Nu objective at converged solution
  ps_mat <- do.call(cbind, lapply(ebmr_b$ps_fit.list, function(pf) pf$fitted.values))
  h_x <- cbind(1, h_nu(dat_b))
  ps_nu <- as.vector(ps_mat %*% res_b$nu.hat)
  g_mat <- as.vector(dat_b$r / ps_nu - 1) * h_x
  G_vec <- colMeans(g_mat)
  W_hat <- tryCatch(solve(crossprod(g_mat) / n_b), error = function(e) diag(ncol(h_x)))
  nu_obj <- as.numeric(t(G_vec) %*% W_hat %*% G_vec)
  neg_r_ps2 <- -dat_b$r / (ps_nu^2)
  Gamma_nu <- crossprod(h_x * neg_r_ps2, ps_mat) / n_b
  nu_grad <- max(abs(2 * as.vector(crossprod(Gamma_nu, W_hat %*% G_vec))))
  cat(sprintf("Rep %d, Nu GMM: obj=%.2e, grad=%.2e\n\n", b, nu_obj, nu_grad))
}

# Also check a few M1-dominated reps for comparison
cat("=== Comparison: M1-dominated reps ===\n\n")
good_count <- 0
for (b in 1:200) {
  if (good_count >= 3) break
  set.seed(12345 + b)
  idx <- sample(1:n, n, replace = TRUE)
  dat_b <- dat[idx, ]
  tryCatch({
    ebmr_b <- EBMRAlgorithmFast4$new("teacher_report", ps_specifications, dat_b, W)
    res_b <- ebmr_b$EBMR_IPW(h_nu = h_nu, true_ps = NULL, se.fit = FALSE)
    if (res_b$w.hat[1] > 0.99) {
      good_count <- good_count + 1
      ps_mat <- do.call(cbind, lapply(ebmr_b$ps_fit.list, function(pf) pf$fitted.values))
      h_x <- cbind(1, h_nu(dat_b))
      n_b <- nrow(dat_b)
      ps_nu <- as.vector(ps_mat %*% res_b$nu.hat)
      g_mat <- as.vector(dat_b$r / ps_nu - 1) * h_x
      G_vec <- colMeans(g_mat)
      W_hat <- tryCatch(solve(crossprod(g_mat) / n_b), error = function(e) diag(ncol(h_x)))
      nu_obj <- as.numeric(t(G_vec) %*% W_hat %*% G_vec)
      neg_r_ps2 <- -dat_b$r / (ps_nu^2)
      Gamma_nu <- crossprod(h_x * neg_r_ps2, ps_mat) / n_b
      nu_grad <- max(abs(2 * as.vector(crossprod(Gamma_nu, W_hat %*% G_vec))))
      cat(sprintf("Rep %d: w=(%s), mu=%.4f, nu_obj=%.2e, nu_grad=%.2e\n",
          b, paste(round(res_b$w.hat, 4), collapse=", "), res_b$mu_ipw, nu_obj, nu_grad))
    }
  }, error = function(e) {})
}
