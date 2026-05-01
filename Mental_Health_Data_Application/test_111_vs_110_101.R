## Why does 111 have severe switching but 110 and 101 don't?
## Same bootstrap samples, compare weights and mu across estimators
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

ps_spec <- list(
  formula.list = list(
    r ~ teacher_report + health + father + teacher_report:health + teacher_report:father,
    r ~ teacher_report + health + parent_report + teacher_report:health + teacher_report:parent_report,
    r ~ teacher_report + parent_report + father + teacher_report:parent_report + teacher_report:father
  ),
  h_alpha.list = list(h_alpha, h_alpha, h_alpha),
  outcome = "teacher_report", inv_link = function(eta) 1/(1+exp(eta)), optimizer = "L-BFGS-B"
)

B <- 300
cat("=== Comparing 110, 101, 111 on same bootstrap samples (B=300) ===\n\n")

mu_110 <- mu_101 <- mu_111 <- rep(NA, B)
w_110 <- matrix(NA, B, 2)  # w for M1, M2
w_101 <- matrix(NA, B, 2)  # w for M1, M3
w_111 <- matrix(NA, B, 3)  # w for M1, M2, M3
m1_degen <- rep(NA, B)

for (b in 1:B) {
  if (b %% 100 == 0) cat(sprintf("  rep %d...\n", b))
  set.seed(12345 + b)
  idx <- sample(1:n, n, replace = TRUE)
  dat_b <- dat[idx, ]

  tryCatch({
    ebmr_b <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat_b, W)

    # Check M1 degeneracy
    m1_ps <- ebmr_b$ps_fit.list[[1]]$fitted.values
    m1_degen[b] <- sum(m1_ps > 0.99) > 0 || sum(m1_ps < 0.01) > 0

    # 110 (M1+M2)
    r110 <- ebmr_b$EBMR_IPW(h_nu = h_nu, model_indices = c(1,2), true_ps = NULL, se.fit = FALSE)
    mu_110[b] <- r110$mu_ipw
    w_110[b, ] <- r110$w.hat

    # 101 (M1+M3)
    r101 <- ebmr_b$EBMR_IPW(h_nu = h_nu, model_indices = c(1,3), true_ps = NULL, se.fit = FALSE)
    mu_101[b] <- r101$mu_ipw
    w_101[b, ] <- r101$w.hat

    # 111 (M1+M2+M3)
    r111 <- ebmr_b$EBMR_IPW(h_nu = h_nu, true_ps = NULL, se.fit = FALSE)
    mu_111[b] <- r111$mu_ipw
    w_111[b, ] <- r111$w.hat
  }, error = function(e) {})
}

valid <- !is.na(mu_111)
cat(sprintf("\nValid: %d/%d\n", sum(valid), B))
cat(sprintf("M1 degenerate: %d/%d (%.1f%%)\n\n", sum(m1_degen, na.rm=TRUE), sum(valid), 100*mean(m1_degen, na.rm=TRUE)))

cat("=== Bootstrap SE ===\n")
cat(sprintf("  110: %.4f\n", sd(mu_110, na.rm=TRUE)))
cat(sprintf("  101: %.4f\n", sd(mu_101, na.rm=TRUE)))
cat(sprintf("  111: %.4f\n\n", sd(mu_111, na.rm=TRUE)))

cat("=== Weight stability ===\n")
cat(sprintf("  110: w1>0.9: %d/%d (%.1f%%), w1 sd=%.4f\n",
    sum(w_110[valid,1]>0.9), sum(valid), 100*mean(w_110[valid,1]>0.9), sd(w_110[valid,1])))
cat(sprintf("  101: w1>0.9: %d/%d (%.1f%%), w1 sd=%.4f\n",
    sum(w_101[valid,1]>0.9), sum(valid), 100*mean(w_101[valid,1]>0.9), sd(w_101[valid,1])))
cat(sprintf("  111: w1>0.9: %d/%d (%.1f%%), w1 sd=%.4f\n\n",
    sum(w_111[valid,1]>0.9), sum(valid), 100*mean(w_111[valid,1]>0.9), sd(w_111[valid,1])))

# When M1 degenerates, what happens to each estimator?
degen <- valid & m1_degen
non_degen <- valid & !m1_degen
cat("=== When M1 degenerates ===\n")
cat(sprintf("  110 mu sd: %.4f (degen) vs %.4f (non-degen)\n",
    sd(mu_110[degen], na.rm=TRUE), sd(mu_110[non_degen], na.rm=TRUE)))
cat(sprintf("  101 mu sd: %.4f (degen) vs %.4f (non-degen)\n",
    sd(mu_101[degen], na.rm=TRUE), sd(mu_101[non_degen], na.rm=TRUE)))
cat(sprintf("  111 mu sd: %.4f (degen) vs %.4f (non-degen)\n\n",
    sd(mu_111[degen], na.rm=TRUE), sd(mu_111[non_degen], na.rm=TRUE)))

cat("=== When M1 degenerates, weight distributions ===\n")
cat(sprintf("  110 w1: mean=%.4f, sd=%.4f\n", mean(w_110[degen,1], na.rm=TRUE), sd(w_110[degen,1], na.rm=TRUE)))
cat(sprintf("  101 w1: mean=%.4f, sd=%.4f\n", mean(w_101[degen,1], na.rm=TRUE), sd(w_101[degen,1], na.rm=TRUE)))
cat(sprintf("  111 w1: mean=%.4f, sd=%.4f\n\n", mean(w_111[degen,1], na.rm=TRUE), sd(w_111[degen,1], na.rm=TRUE)))

# Do 110 and 101 still give M1 high weight when M1 degenerates?
cat(sprintf("  When M1 degen, 110 w1>0.9: %d/%d (%.1f%%)\n",
    sum(w_110[degen,1]>0.9, na.rm=TRUE), sum(degen), 100*mean(w_110[degen,1]>0.9, na.rm=TRUE)))
cat(sprintf("  When M1 degen, 101 w1>0.9: %d/%d (%.1f%%)\n",
    sum(w_101[degen,1]>0.9, na.rm=TRUE), sum(degen), 100*mean(w_101[degen,1]>0.9, na.rm=TRUE)))
cat(sprintf("  When M1 degen, 111 w1>0.9: %d/%d (%.1f%%)\n\n",
    sum(w_111[degen,1]>0.9, na.rm=TRUE), sum(degen), 100*mean(w_111[degen,1]>0.9, na.rm=TRUE)))

# Show a few reps where 110/101 keep M1 but 111 switches
cat("=== Examples: 110/101 keep M1 but 111 switches ===\n")
count <- 0
for (b in which(valid)) {
  if (count >= 10) break
  if (w_110[b,1] > 0.9 && w_101[b,1] > 0.9 && w_111[b,1] < 0.5) {
    count <- count + 1
    cat(sprintf("  Rep %d: 110 w1=%.4f mu=%.4f | 101 w1=%.4f mu=%.4f | 111 w=(%s) mu=%.4f | M1_degen=%s\n",
        b, w_110[b,1], mu_110[b], w_101[b,1], mu_101[b],
        paste(round(w_111[b,], 4), collapse=","), mu_111[b], m1_degen[b]))
  }
}
if (count == 0) cat("  None found\n")
