## Examine nu estimation for 111 ensemble across bootstrap reps
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

B <- 500
cat("=== 111 ensemble nu diagnostics (B=500) ===\n\n")

mus <- rep(NA, B)
nu_mat <- matrix(NA, B, 3)
w_mat <- matrix(NA, B, 3)

for (b in 1:B) {
  if (b %% 100 == 0) cat(sprintf("  rep %d...\n", b))
  set.seed(12345 + b)
  idx <- sample(1:n, n, replace = TRUE)
  dat_b <- dat[idx, ]

  tryCatch({
    ebmr_b <- EBMRAlgorithmFast4$new("teacher_report", ps_specifications, dat_b, W)
    res_b <- ebmr_b$EBMR_IPW(h_nu = h_nu, true_ps = NULL, se.fit = FALSE)
    mus[b] <- res_b$mu_ipw
    nu_mat[b, ] <- res_b$nu.hat
    w_mat[b, ] <- res_b$w.hat
  }, error = function(e) {})
}

valid <- !is.na(mus)
cat(sprintf("\nValid: %d/%d\n", sum(valid), B))
cat(sprintf("mu: mean=%.4f, sd=%.4f\n\n", mean(mus, na.rm=TRUE), sd(mus, na.rm=TRUE)))

# Nu distribution
cat("=== Nu distribution ===\n")
for (j in 1:3) {
  nv <- nu_mat[valid, j]
  cat(sprintf("  nu%d: mean=%.4f, sd=%.4f, min=%.4f, Q25=%.4f, med=%.4f, Q75=%.4f, max=%.4f\n",
      j, mean(nv), sd(nv), min(nv), quantile(nv, 0.25), median(nv), quantile(nv, 0.75), max(nv)))
}

# Weight distribution
cat("\n=== Weight distribution ===\n")
for (j in 1:3) {
  wv <- w_mat[valid, j]
  cat(sprintf("  w%d: mean=%.4f, sd=%.4f, min=%.4f, Q25=%.4f, med=%.4f, Q75=%.4f, max=%.4f\n",
      j, mean(wv), sd(wv), min(wv), quantile(wv, 0.25), median(wv), quantile(wv, 0.75), max(wv)))
}
cat(sprintf("\n  w1 > 0.9: %d/%d (%.1f%%)\n", sum(w_mat[valid,1] > 0.9), sum(valid), 100*mean(w_mat[valid,1] > 0.9)))
cat(sprintf("  w2 > 0.9: %d/%d (%.1f%%)\n", sum(w_mat[valid,2] > 0.9), sum(valid), 100*mean(w_mat[valid,2] > 0.9)))
cat(sprintf("  w3 > 0.9: %d/%d (%.1f%%)\n", sum(w_mat[valid,3] > 0.9), sum(valid), 100*mean(w_mat[valid,3] > 0.9)))

# Correlation with mu
cat(sprintf("\n  cor(w1, mu): %.3f\n", cor(w_mat[valid,1], mus[valid])))
cat(sprintf("  cor(w2, mu): %.3f\n", cor(w_mat[valid,2], mus[valid])))
cat(sprintf("  cor(w3, mu): %.3f\n", cor(w_mat[valid,3], mus[valid])))

# What happens when w1 < 0.9 (non-M1-dominated)?
non_m1 <- valid & w_mat[,1] < 0.9
cat(sprintf("\n=== Reps where M1 weight < 0.9: %d/%d ===\n", sum(non_m1), sum(valid)))
if (sum(non_m1) > 0) {
  cat(sprintf("  mu in these reps: mean=%.4f, sd=%.4f\n", mean(mus[non_m1]), sd(mus[non_m1])))
  cat(sprintf("  mu in M1-dominated reps: mean=%.4f, sd=%.4f\n",
      mean(mus[valid & w_mat[,1] >= 0.9]), sd(mus[valid & w_mat[,1] >= 0.9])))
  # Show a few examples
  non_m1_idx <- which(non_m1)
  cat("\n  Examples:\n")
  for (i in head(non_m1_idx, 10)) {
    cat(sprintf("    rep %d: mu=%.4f, w=(%.4f, %.4f, %.4f)\n",
        i, mus[i], w_mat[i,1], w_mat[i,2], w_mat[i,3]))
  }
}
