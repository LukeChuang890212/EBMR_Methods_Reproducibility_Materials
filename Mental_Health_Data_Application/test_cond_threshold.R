## Test: Lower cond_threshold to reject degenerate solutions
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Mental_Health_Data_Application")
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)
library(Matrix)
library(numDeriv)
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
    r ~ teacher_report + health + father + teacher_report:health + teacher_report:father,
    r ~ teacher_report + health + parent_report + teacher_report:health + teacher_report:parent_report,
    r ~ teacher_report + parent_report + father + teacher_report:parent_report + teacher_report:father
  ),
  h_alpha.list = list(
    c("health", "father", "parent_report", "fp", "fh", "hp"),
    c("health", "father", "parent_report", "fp", "fh", "hp"),
    c("health", "father", "parent_report", "fp", "fh", "hp")
  ),
  outcome = "teacher_report",
  inv_link = function(eta) 1 / (1 + exp(eta))
)

h_nu <- function(data) {
  cbind(health = data$health, father = data$father, parent_report = data$parent_report,
        fp = data$fp, fh = data$fh, hp = data$hp, fhp = data$fhp)
}

# Test different cond_threshold values
thresholds <- c(1e8, 1e6, 1e4, 1e3, 1e2)

cat("=== Testing different cond_threshold values ===\n\n")

for (thresh in thresholds) {
  cat(sprintf("--- cond_threshold = %.0e ---\n", thresh))

  ps_spec_t <- ps_specifications
  ps_spec_t$cond_threshold <- thresh

  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("teacher_report", ps_spec_t, dat, W)

    for (i in 1:3) {
      pf <- ebmr$ps_fit.list[[i]]
      ps <- pf$fitted.values
      cat(sprintf("  Model %d: alpha=(%s)\n", i,
          paste(round(pf$coefficients, 3), collapse=", ")))
      cat(sprintf("           PS range=[%.4f, %.4f], PS>0.99: %d, PS<0.01: %d\n",
          min(ps), max(ps), sum(ps > 0.99), sum(ps < 0.01)))
      cat(sprintf("           GMM obj=%.2e, grad=%.2e, cond=%.2e\n",
          pf$gmm_fit$opt$objective, pf$gmm_fit$opt$final_grad_norm,
          tryCatch({
            Gamma <- pf$gmm_fit$Gamma.hat
            W_hat <- pf$gmm_fit$W.hat
            H <- crossprod(Gamma, W_hat %*% Gamma)
            ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
            max(ev) / max(min(ev), 1e-15)
          }, error = function(e) NA)))
    }

    result <- ebmr$EBMR_IPW(h_nu = h_nu, true_ps = NULL)
    cat(sprintf("  Ensemble: mu=%.4f, se=%.4f, w=(%s)\n",
        result$mu_ipw, result$se_ipw,
        paste(round(result$w.hat, 4), collapse=", ")))

    # Single model results
    for (idx in 1:3) {
      res_s <- ebmr$EBMR_IPW(h_nu = h_nu, model_indices = idx, true_ps = NULL)
      cat(sprintf("  M%d only:  mu=%.4f, se=%.4f\n", idx, res_s$mu_ipw, res_s$se_ipw))
    }

  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
  })
  cat("\n")
}
