setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Mental_Health_Data_Application")
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)
library(Matrix); library(numDeriv)
source("MHD_functions.R")

original_data <- read.csv("data_application.csv")
percent <- original_data$Percentage
n_full <- 2486
class_n <- round(n_full * percent / 100)
dat_full <- gen_data(original_data, class_n, n_full)
dat_full$y <- dat_full$teacher_report
dat <- dat_full[dat_full$health == 1, ]
n <- nrow(dat)

W <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

ps_specifications <- list(
  formula.list = list(
    r ~ teacher_report + father,
    r ~ teacher_report + parent_report,
    r ~ father + parent_report
  ),
  h_alpha.list = list(
    c("father", "parent_report"),
    c("father", "parent_report"),
    c("father", "parent_report")
  ),
  outcome = "teacher_report",
  inv_link = function(eta) 1 / (1 + exp(eta))
)

h_nu <- function(data) {
  cbind(father = data$father, parent_report = data$parent_report,
        fp = data$father * data$parent_report)
}

# Test single perturbation rep
set.seed(12346)
wt <- rexp(n, rate = 1)

cat("=== Debug perturbation bootstrap ===\n\n")
cat(sprintf("n = %d, wt range = [%.3f, %.3f]\n\n", n, min(wt), max(wt)))

# Try without wt first
cat("--- Without wt (should work) ---\n")
tryCatch({
  ebmr0 <- EBMRAlgorithmFast4$new("teacher_report", ps_specifications, dat, W)
  res0 <- ebmr0$EBMR_IPW(h_nu = h_nu, true_ps = NULL, se.fit = FALSE)
  cat(sprintf("mu_ipw = %.4f\n", res0$mu_ipw))
}, error = function(e) cat(sprintf("ERROR: %s\n", e$message)))

# Try with wt
cat("\n--- With wt (perturbation) ---\n")
tryCatch({
  ebmr1 <- EBMRAlgorithmFast4$new("teacher_report", ps_specifications, dat, W, wt = wt)
  cat("Constructor succeeded.\n")
  for (i in 1:3) {
    pf <- ebmr1$ps_fit.list[[i]]
    cat(sprintf("  M%d: coef=(%s), PS=[%.4f,%.4f]\n", i,
        paste(round(pf$coefficients, 3), collapse=", "),
        min(pf$fitted.values), max(pf$fitted.values)))
  }
  res1 <- ebmr1$EBMR_IPW(h_nu = h_nu, true_ps = NULL, se.fit = FALSE, wt = wt)
  cat(sprintf("mu_ipw = %.4f\n", res1$mu_ipw))
}, error = function(e) {
  cat(sprintf("ERROR: %s\n", e$message))
  cat("Traceback:\n")
  traceback()
})

# Try each model separately with wt
cat("\n--- Each model separately with wt ---\n")
for (i in 1:3) {
  ps_spec_i <- list(
    formula.list = list(ps_specifications$formula.list[[i]]),
    h_alpha.list = list(ps_specifications$h_alpha.list[[i]]),
    outcome = ps_specifications$outcome,
    inv_link = ps_specifications$inv_link
  )
  tryCatch({
    ebmr_i <- EBMRAlgorithmFast4$new("teacher_report", ps_spec_i, dat, W, wt = wt)
    pf <- ebmr_i$ps_fit.list[[1]]
    cat(sprintf("M%d: coef=(%s), PS=[%.4f,%.4f]\n", i,
        paste(round(pf$coefficients, 3), collapse=", "),
        min(pf$fitted.values), max(pf$fitted.values)))
    res_i <- ebmr_i$EBMR_IPW(h_nu = h_nu, true_ps = NULL, se.fit = FALSE, wt = wt)
    cat(sprintf("  mu_ipw = %.4f\n", res_i$mu_ipw))
  }, error = function(e) {
    cat(sprintf("M%d ERROR: %s\n", i, e$message))
  })
}
