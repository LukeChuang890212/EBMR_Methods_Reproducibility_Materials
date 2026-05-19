## Quick diagnostic: compare Fast4 alpha estimates with Fast2 reference
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

ebmr <- EBMRAlgorithmFast4$new("teacher_report", ps_specifications, dat, W)

cat("=== Model diagnostics ===\n\n")
for (i in 1:3) {
  pf <- ebmr$ps_fit.list[[i]]
  cat(sprintf("Model %d:\n", i))
  cat(sprintf("  Formula: %s\n", deparse(ps_specifications$formula.list[[i]])))
  cat(sprintf("  Design matrix cols: %s\n", paste(colnames(pf$design_matrix), collapse=", ")))
  cat(sprintf("  k=%d, h_dim=%d\n", ncol(pf$design_matrix), ncol(pf$h_x)))
  cat(sprintf("  alpha = (%s)\n", paste(round(pf$coefficients, 4), collapse=", ")))
  cat(sprintf("  SE    = (%s)\n", paste(round(pf$se, 4), collapse=", ")))
  cat(sprintf("  PS range: [%.6f, %.6f]\n", min(pf$fitted.values), max(pf$fitted.values)))
  cat(sprintf("  GMM converged: %s, obj=%.2e, grad=%.2e\n",
      pf$gmm_fit$opt$converged, pf$gmm_fit$opt$objective, pf$gmm_fit$opt$final_grad_norm))

  # Check for extreme PS
  ps <- pf$fitted.values
  cat(sprintf("  PS < 0.01: %d, PS > 0.99: %d\n", sum(ps < 0.01), sum(ps > 0.99)))

  # IPW mean with this single model
  mu_single <- mean(dat$r * dat$teacher_report / ps)
  cat(sprintf("  Single-model IPW mean: %.4f\n", mu_single))
  cat("\n")
}

# Check Fast2 reference if available
cat("=== Reference (Fast2 saved results) ===\n")
ref <- readRDS("MHD_results/popmean_summary.RDS")
print(ref)

# Check ensemble weights for all combinations
cat("\n=== Ensemble weights ===\n")
model_sets <- list(c(1,2), c(1,3), c(2,3), c(1,2,3))
labels <- c("12", "13", "23", "123")
for (j in seq_along(model_sets)) {
  ms <- model_sets[[j]]
  res <- ebmr$EBMR_IPW(h_nu = h_nu, model_indices = ms, true_ps = NULL)
  cat(sprintf("  %s: nu=(%s), w=(%s), mu=%.4f\n",
      labels[j],
      paste(round(res$nu.hat, 4), collapse=", "),
      paste(round(res$w.hat, 4), collapse=", "),
      res$mu_ipw))
}
