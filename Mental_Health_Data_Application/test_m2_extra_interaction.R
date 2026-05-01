## Test M2 with added health:parent_report interaction
## M2: r ~ teacher_report + health + parent_report + teacher_report:health + teacher_report:parent_report + health:parent_report
## k=7, h_dim=7 (with intercept), overid=0
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
h_alpha <- c("health", "father", "parent_report", "fp", "fh", "hp")
h_nu <- function(data) cbind(health=data$health, father=data$father, parent_report=data$parent_report,
                              fp=data$fp, fh=data$fh, hp=data$hp, fhp=data$fhp)

# M2 with health:parent_report added
ps_spec_m2 <- list(
  formula.list = list(
    r ~ teacher_report + health + parent_report + teacher_report:health + teacher_report:parent_report + health:parent_report
  ),
  h_alpha.list = list(h_alpha),
  outcome = "teacher_report",
  inv_link = function(eta) 1 / (1 + exp(eta))
)

cat("=== M2 with health:parent_report added ===\n\n")

# Point estimate
ebmr <- EBMRAlgorithmFast4$new("teacher_report", ps_spec_m2, dat, W)
pf <- ebmr$ps_fit.list[[1]]
ps <- pf$fitted.values

cat(sprintf("Design cols: %s\n", paste(colnames(pf$design_matrix), collapse=", ")))
cat(sprintf("k=%d, h_dim=%d, overid=%d\n", ncol(pf$design_matrix), ncol(pf$h_x),
    ncol(pf$h_x) - ncol(pf$design_matrix)))
cat(sprintf("alpha = (%s)\n", paste(round(pf$coefficients, 4), collapse=", ")))
cat(sprintf("SE    = (%s)\n", paste(round(pf$se, 4), collapse=", ")))
cat(sprintf("PS range: [%.6f, %.6f], >0.99:%d, <0.01:%d\n", min(ps), max(ps), sum(ps>0.99), sum(ps<0.01)))
cat(sprintf("GMM obj=%.2e, grad=%.2e\n", pf$gmm_fit$opt$objective, pf$gmm_fit$opt$final_grad_norm))

res <- ebmr$EBMR_IPW(h_nu = h_nu, true_ps = NULL)
cat(sprintf("mu_ipw = %.4f, se_ipw = %.4f\n\n", res$mu_ipw, res$se_ipw))

# Perturbation bootstrap
cat("=== Perturbation bootstrap (B=1000) ===\n\n")
B <- 1000
mus <- rep(NA, B)

for (b in 1:B) {
  if (b %% 100 == 0) cat(sprintf("  rep %d...\n", b))
  set.seed(12345 + b)
  wt <- rexp(n, rate = 1)
  tryCatch({
    ebmr_b <- EBMRAlgorithmFast4$new("teacher_report", ps_spec_m2, dat, W, wt = wt)
    res_b <- ebmr_b$EBMR_IPW(h_nu = h_nu, true_ps = NULL, se.fit = FALSE, wt = wt)
    mus[b] <- res_b$mu_ipw
  }, error = function(e) {})
}

valid <- !is.na(mus)
ptb_se <- sd(mus[valid])
cat(sprintf("\nValid: %d/%d\n", sum(valid), B))
cat(sprintf("Analytical SE: %.4f\n", res$se_ipw))
cat(sprintf("Perturbation SE: %.4f\n", ptb_se))
cat(sprintf("Ratio (Analytical/Perturb): %.3f\n", res$se_ipw / ptb_se))
