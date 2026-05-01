## Test new PS model specifications for MHD overall mean analysis
## Model 1: r ~ teacher_report + health + father + father:health
## Model 2: r ~ teacher_report + health + parent_report + health:parent_report
## Model 3: r ~ teacher_report + parent_report + father + parent_report:father
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
    r ~ teacher_report + health + father + father:health,
    r ~ teacher_report + health + parent_report + health:parent_report,
    r ~ teacher_report + parent_report + father + parent_report:father
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

cat("=== New PS model specifications ===\n\n")

ebmr <- EBMRAlgorithmFast4$new("teacher_report", ps_specifications, dat, W)

for (i in 1:3) {
  pf <- ebmr$ps_fit.list[[i]]
  ps <- pf$fitted.values
  cat(sprintf("Model %d: %s\n", i, deparse(ps_specifications$formula.list[[i]])))
  cat(sprintf("  Design cols: %s\n", paste(colnames(pf$design_matrix), collapse=", ")))
  cat(sprintf("  k=%d, h_dim=%d, overid=%d\n", ncol(pf$design_matrix), ncol(pf$h_x),
      ncol(pf$h_x) - ncol(pf$design_matrix)))
  cat(sprintf("  alpha = (%s)\n", paste(round(pf$coefficients, 4), collapse=", ")))
  cat(sprintf("  SE    = (%s)\n", paste(round(pf$se, 4), collapse=", ")))
  cat(sprintf("  PS range: [%.6f, %.6f], >0.99: %d, <0.01: %d\n",
      min(ps), max(ps), sum(ps > 0.99), sum(ps < 0.01)))
  cat(sprintf("  GMM obj=%.2e, grad=%.2e\n", pf$gmm_fit$opt$objective, pf$gmm_fit$opt$final_grad_norm))
  mu_single <- mean(dat$r * dat$teacher_report / ps)
  cat(sprintf("  Single-model IPW mean: %.4f\n\n", mu_single))
}

# All model combinations
mu_cc <- mean(dat$teacher_report[dat$r == 1])
se_cc <- sd(dat$teacher_report[dat$r == 1]) / sqrt(sum(dat$r))

cat("=== All model combinations ===\n\n")
cat(sprintf("  %-20s %10s %10s %20s\n", "Models", "Estimate", "SE", "w.hat"))
cat("  ", paste(rep("-", 65), collapse = ""), "\n")
cat(sprintf("  %-20s %10.4f %10.4f\n", "Complete Case", mu_cc, se_cc))

all_model_sets <- list(c(1), c(2), c(3), c(1,2), c(1,3), c(2,3), c(1,2,3))
model_names <- c("M1", "M2", "M3", "M1+M2", "M1+M3", "M2+M3", "M1+M2+M3")

for (j in seq_along(all_model_sets)) {
  ms <- all_model_sets[[j]]
  res <- ebmr$EBMR_IPW(h_nu = h_nu, model_indices = ms, true_ps = NULL)
  w_str <- paste(round(res$w.hat, 4), collapse=", ")
  cat(sprintf("  %-20s %10.4f %10.4f     (%s)\n",
      model_names[j], res$mu_ipw, res$se_ipw, w_str))
}
cat("  ", paste(rep("-", 65), collapse = ""), "\n")
