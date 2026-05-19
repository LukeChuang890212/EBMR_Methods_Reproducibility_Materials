## Investigate: Does M1 with k=6 (no father:health) stay stable across bootstrap?
## If so, does 111 ensemble improve?
## Compare k=7 vs k=6 for M1 stability
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
h_nu <- function(data) cbind(health=data$health, father=data$father, parent_report=data$parent_report,
                              fp=data$fp, fh=data$fh, hp=data$hp, fhp=data$fhp)

# k=7 formulas (current — M1 degenerates in bootstrap)
ps_k7 <- list(
  formula.list = list(
    r ~ teacher_report + health + father + teacher_report:health + teacher_report:father + father:health,
    r ~ teacher_report + parent_report + health + teacher_report:parent_report + teacher_report:health + parent_report:health,
    r ~ teacher_report + parent_report + father + teacher_report:parent_report + teacher_report:father + parent_report:father
  ),
  h_alpha.list = list(h_alpha, h_alpha, h_alpha),
  outcome = "teacher_report", inv_link = function(eta) 1/(1+exp(eta)), optimizer = "L-BFGS-B"
)

# k=6 formulas (original — M1 was stable)
ps_k6 <- list(
  formula.list = list(
    r ~ teacher_report + health + father + teacher_report:health + teacher_report:father,
    r ~ teacher_report + health + parent_report + teacher_report:health + teacher_report:parent_report,
    r ~ teacher_report + parent_report + father + teacher_report:parent_report + teacher_report:father
  ),
  h_alpha.list = list(h_alpha, h_alpha, h_alpha),
  outcome = "teacher_report", inv_link = function(eta) 1/(1+exp(eta)), optimizer = "L-BFGS-B"
)

B <- 300
cat("=== Comparing M1 stability: k=7 vs k=6 (B=300) ===\n\n")

for (spec_name in c("k=7", "k=6")) {
  ps_spec <- if (spec_name == "k=7") ps_k7 else ps_k6
  cat(sprintf("--- %s formulas ---\n", spec_name))

  # Point estimate
  ebmr <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat, W)
  res <- ebmr$EBMR_IPW(h_nu = h_nu, true_ps = NULL)
  m1_ps <- ebmr$ps_fit.list[[1]]$fitted.values
  cat(sprintf("  Point: mu=%.4f, se=%.4f, w=(%s)\n",
      res$mu_ipw, res$se_ipw, paste(round(res$w.hat, 4), collapse=",")))
  cat(sprintf("  M1 PS: [%.4f, %.4f], >0.99:%d\n\n",
      min(m1_ps), max(m1_ps), sum(m1_ps > 0.99)))

  # Bootstrap
  mus <- rep(NA, B)
  w_mat <- matrix(NA, B, 3)
  m1_degen <- rep(NA, B)

  for (b in 1:B) {
    if (b %% 100 == 0) cat(sprintf("  rep %d...\n", b))
    set.seed(12345 + b)
    idx <- sample(1:n, n, replace = TRUE)
    dat_b <- dat[idx, ]
    tryCatch({
      ebmr_b <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat_b, W)
      res_b <- ebmr_b$EBMR_IPW(h_nu = h_nu, true_ps = NULL, se.fit = FALSE)
      mus[b] <- res_b$mu_ipw
      w_mat[b, ] <- res_b$w.hat
      m1_ps_b <- ebmr_b$ps_fit.list[[1]]$fitted.values
      m1_degen[b] <- sum(m1_ps_b > 0.99) > 0 || sum(m1_ps_b < 0.01) > 0
    }, error = function(e) {})
  }

  valid <- !is.na(mus)
  cat(sprintf("  Bootstrap SE: %.4f (analytical: %.4f, ratio: %.3f)\n",
      sd(mus, na.rm=TRUE), res$se_ipw, res$se_ipw / sd(mus, na.rm=TRUE)))
  cat(sprintf("  M1 degenerate in bootstrap: %d/%d (%.1f%%)\n",
      sum(m1_degen, na.rm=TRUE), sum(valid), 100*mean(m1_degen, na.rm=TRUE)))
  cat(sprintf("  w1 > 0.9: %d/%d (%.1f%%)\n",
      sum(w_mat[valid,1] > 0.9), sum(valid), 100*mean(w_mat[valid,1] > 0.9)))
  cat(sprintf("  w2 > 0.9: %d/%d (%.1f%%)\n",
      sum(w_mat[valid,2] > 0.9), sum(valid), 100*mean(w_mat[valid,2] > 0.9)))
  cat(sprintf("  w3 > 0.9: %d/%d (%.1f%%)\n\n",
      sum(w_mat[valid,3] > 0.9), sum(valid), 100*mean(w_mat[valid,3] > 0.9)))
}
