## Verify: no change when perturbation is not used
## Compare results with and without the wt fixes by running without wt
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

# Expected values (from previous runs before the fix)
expected <- list(
  `100` = list(mu = 0.1434, se = 0.0301),
  `010` = list(mu = 0.4679, se = 0.1964),
  `001` = list(mu = 0.1985, se = 0.0156),
  `110` = list(mu = 0.1435, se = 0.0301),
  `101` = list(mu = 0.1434, se = 0.0300),
  `011` = list(mu = 0.1985, se = 0.0156),
  `111` = list(mu = 0.1435, se = 0.0296)
)

cat("=== Verify no change without perturbation ===\n\n")

ebmr <- EBMRAlgorithmFast4$new("teacher_report", ps_specifications, dat, W)

all_model_sets <- list(c(1), c(2), c(3), c(1,2), c(1,3), c(2,3), c(1,2,3))
labels <- c("100", "010", "001", "110", "101", "011", "111")

all_pass <- TRUE
cat(sprintf("  %-8s %10s %10s %10s %10s %8s\n",
            "Label", "mu_now", "mu_exp", "se_now", "se_exp", "Match"))
cat("  ", paste(rep("-", 62), collapse = ""), "\n")

for (j in 1:7) {
  res <- ebmr$EBMR_IPW(h_nu = h_nu, model_indices = all_model_sets[[j]], true_ps = NULL)
  exp <- expected[[labels[j]]]
  mu_match <- abs(res$mu_ipw - exp$mu) < 0.0005
  se_match <- abs(res$se_ipw - exp$se) < 0.0005
  pass <- mu_match & se_match
  if (!pass) all_pass <- FALSE
  cat(sprintf("  %-8s %10.4f %10.4f %10.4f %10.4f %8s\n",
              labels[j], res$mu_ipw, exp$mu, res$se_ipw, exp$se,
              if (pass) "OK" else "FAIL"))
}

cat("  ", paste(rep("-", 62), collapse = ""), "\n")
cat(sprintf("\n  Overall: %s\n", if (all_pass) "ALL PASSED" else "SOME FAILED"))
