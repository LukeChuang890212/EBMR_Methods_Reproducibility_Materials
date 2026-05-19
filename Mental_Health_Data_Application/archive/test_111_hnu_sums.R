## Test 111 with h_nu = (father+health, pr+father, pr+health)
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

# h_nu: intercept (auto) + father+health, pr+father, pr+health
h_nu <- function(data) cbind(fh_sum = data$father + data$health,
                              pf_sum = data$parent_report + data$father,
                              ph_sum = data$parent_report + data$health)

ebmr <- EBMRAlgorithmFast4$new("teacher_report", ps_specifications, dat, W)
res <- ebmr$EBMR_IPW(h_nu = h_nu, true_ps = NULL)
cat(sprintf("Point: mu=%.4f, se=%.4f, w=(%s)\n\n",
    res$mu_ipw, res$se_ipw, paste(round(res$w.hat, 4), collapse=", ")))

B <- 500
mus <- rep(NA, B)
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
    w_mat[b, ] <- res_b$w.hat
  }, error = function(e) {})
}

boot_se <- sd(mus, na.rm = TRUE)
valid <- !is.na(mus)
cat(sprintf("\nAnalytical SE: %.4f\nBootstrap SE: %.4f\nRatio: %.3f\n\n",
    res$se_ipw, boot_se, res$se_ipw / boot_se))

cat("Weight distribution:\n")
for (j in 1:3) cat(sprintf("  w%d: mean=%.4f, sd=%.4f, >0.9: %d/%d\n",
    j, mean(w_mat[valid,j]), sd(w_mat[valid,j]), sum(w_mat[valid,j]>0.9), sum(valid)))
