## Perturbation bootstrap for 111 with k=6 formulas
## Perturbation preserves data composition — M1 should stay stable
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

# k=6 formulas
ps_spec <- list(
  formula.list = list(
    r ~ teacher_report + health + father + teacher_report:health + teacher_report:father,
    r ~ teacher_report + health + parent_report + teacher_report:health + teacher_report:parent_report,
    r ~ teacher_report + parent_report + father + teacher_report:parent_report + teacher_report:father
  ),
  h_alpha.list = list(h_alpha, h_alpha, h_alpha),
  outcome = "teacher_report", inv_link = function(eta) 1/(1+exp(eta)), optimizer = "L-BFGS-B"
)

# Point estimate
ebmr <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat, W)
res <- ebmr$EBMR_IPW(h_nu = h_nu, true_ps = NULL)
cat(sprintf("Point: mu=%.4f, se=%.4f, w=(%s)\n\n",
    res$mu_ipw, res$se_ipw, paste(round(res$w.hat, 4), collapse=",")))

B <- 500
cat(sprintf("=== Perturbation bootstrap (B=%d) ===\n\n", B))

mus <- rep(NA, B)
w_mat <- matrix(NA, B, 3)
m1_degen <- rep(NA, B)

for (b in 1:B) {
  if (b %% 100 == 0) cat(sprintf("  rep %d...\n", b))
  set.seed(12345 + b)
  wt <- rexp(n, rate = 1)
  tryCatch({
    ebmr_b <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat, W, wt = wt)
    res_b <- ebmr_b$EBMR_IPW(h_nu = h_nu, true_ps = NULL, se.fit = FALSE, wt = wt)
    mus[b] <- res_b$mu_ipw
    w_mat[b, ] <- res_b$w.hat
    m1_ps_b <- ebmr_b$ps_fit.list[[1]]$fitted.values
    m1_degen[b] <- sum(m1_ps_b > 0.99) > 0 || sum(m1_ps_b < 0.01) > 0
  }, error = function(e) {})
}

valid <- !is.na(mus)
ptb_se <- sd(mus[valid])
cat(sprintf("\nValid: %d/%d\n", sum(valid), B))
cat(sprintf("Analytical SE: %.4f\n", res$se_ipw))
cat(sprintf("Perturbation SE: %.4f\n", ptb_se))
cat(sprintf("Ratio: %.3f\n\n", res$se_ipw / ptb_se))

cat(sprintf("M1 degenerate: %d/%d (%.1f%%)\n",
    sum(m1_degen, na.rm=TRUE), sum(valid), 100*mean(m1_degen, na.rm=TRUE)))
cat(sprintf("w1 > 0.9: %d/%d (%.1f%%)\n",
    sum(w_mat[valid,1] > 0.9), sum(valid), 100*mean(w_mat[valid,1] > 0.9)))
cat(sprintf("w2 > 0.9: %d/%d (%.1f%%)\n",
    sum(w_mat[valid,2] > 0.9), sum(valid), 100*mean(w_mat[valid,2] > 0.9)))
cat(sprintf("w3 > 0.9: %d/%d (%.1f%%)\n",
    sum(w_mat[valid,3] > 0.9), sum(valid), 100*mean(w_mat[valid,3] > 0.9)))
