## Perturbation bootstrap for 111 ensemble
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
                              fp=data$fp, fh=data$fh, hp=data$hp)

# Point estimate
ebmr <- EBMRAlgorithmFast4$new("teacher_report", ps_specifications, dat, W)
res <- ebmr$EBMR_IPW(h_nu = h_nu, true_ps = NULL)
cat(sprintf("Point estimate: mu=%.4f, se=%.4f\n\n", res$mu_ipw, res$se_ipw))

# Perturbation bootstrap
B <- 1000
cat(sprintf("=== Perturbation bootstrap (B=%d) ===\n\n", B))

mus <- rep(NA, B)
w_mat <- matrix(NA, B, 3)

for (b in 1:B) {
  if (b %% 200 == 0) cat(sprintf("  rep %d...\n", b))
  set.seed(12345 + b)
  wt <- rexp(n, rate = 1)

  tryCatch({
    ebmr_b <- EBMRAlgorithmFast4$new("teacher_report", ps_specifications, dat, W, wt = wt)
    res_b <- ebmr_b$EBMR_IPW(h_nu = h_nu, true_ps = NULL, se.fit = FALSE, wt = wt)
    mus[b] <- res_b$mu_ipw
    w_mat[b, ] <- res_b$w.hat
  }, error = function(e) {})
}

valid <- !is.na(mus)
ptb_se <- sd(mus[valid])

cat(sprintf("\nValid: %d/%d\n", sum(valid), B))
cat(sprintf("Analytical SE: %.4f\n", res$se_ipw))
cat(sprintf("Perturbation SE: %.4f\n", ptb_se))
cat(sprintf("Ratio: %.3f\n\n", res$se_ipw / ptb_se))

# Weight distribution
cat("Weight distribution:\n")
for (j in 1:3) {
  wv <- w_mat[valid, j]
  cat(sprintf("  w%d: mean=%.4f, sd=%.4f, >0.9: %d/%d (%.1f%%)\n",
      j, mean(wv), sd(wv), sum(wv > 0.9), sum(valid), 100*mean(wv > 0.9)))
}

cat(sprintf("\n  cor(w1, mu): %.3f\n", cor(w_mat[valid,1], mus[valid])))
cat(sprintf("  cor(w2, mu): %.3f\n", cor(w_mat[valid,2], mus[valid])))
cat(sprintf("  cor(w3, mu): %.3f\n", cor(w_mat[valid,3], mus[valid])))
