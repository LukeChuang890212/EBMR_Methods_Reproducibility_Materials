## Systematic search for M2 formulas that are non-degenerate
## M2 must include teacher_report and at least one covariate NOT in M1
## M1 = r ~ teacher_report + health + father + teacher_report:health + teacher_report:father
## M3 = r ~ teacher_report + parent_report + father + teacher_report:parent_report + teacher_report:father
## M2 should be distinct from both
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

cat("=== Systematic M2 formula search ===\n")
cat("M1: r ~ y + health + father + y:health + y:father\n")
cat("M3: r ~ y + parent_report + father + y:parent_report + y:father\n")
cat("Looking for M2 that is non-degenerate and distinct from M1, M3\n\n")

# All candidate M2 formulas involving teacher_report
# Must involve parent_report (to distinguish from M1 which has health+father)
# or health+parent_report combination different from M1 and M3
candidates <- list(
  # parent_report only (no other y-interactions)
  list(f = r ~ teacher_report + parent_report, desc = "y + pr", k = 3),
  list(f = r ~ teacher_report + parent_report + health, desc = "y + pr + health", k = 4),
  list(f = r ~ teacher_report + parent_report + father, desc = "y + pr + father", k = 4),
  list(f = r ~ teacher_report + parent_report + health + father, desc = "y + pr + health + father", k = 5),

  # with teacher_report:parent_report
  list(f = r ~ teacher_report + parent_report + teacher_report:parent_report, desc = "y + pr + y:pr", k = 4),
  list(f = r ~ teacher_report + parent_report + health + teacher_report:parent_report, desc = "y + pr + health + y:pr", k = 5),
  list(f = r ~ teacher_report + parent_report + father + teacher_report:parent_report, desc = "y + pr + father + y:pr", k = 5),

  # with covariate interactions only
  list(f = r ~ teacher_report + parent_report + health + health:parent_report, desc = "y + pr + health + h:pr", k = 5),
  list(f = r ~ teacher_report + parent_report + father + father:parent_report, desc = "y + pr + father + f:pr", k = 5),

  # health-focused (with parent_report to distinguish from M1)
  list(f = r ~ teacher_report + health + parent_report + teacher_report:health, desc = "y + health + pr + y:health", k = 5),
  list(f = r ~ teacher_report + health + parent_report + health:parent_report, desc = "y + health + pr + h:pr", k = 5),

  # father + parent_report combos (to distinguish from M1 and M3)
  list(f = r ~ teacher_report + father + parent_report + father:parent_report, desc = "y + father + pr + f:pr", k = 5),
  list(f = r ~ teacher_report + father + parent_report + teacher_report:father, desc = "y + father + pr + y:father", k = 5),

  # Simple: just teacher_report + one covariate
  list(f = r ~ teacher_report + health, desc = "y + health", k = 3),
  list(f = r ~ teacher_report + father, desc = "y + father", k = 3),

  # With y:health (like M1 but without father)
  list(f = r ~ teacher_report + health + teacher_report:health, desc = "y + health + y:health", k = 4),
  list(f = r ~ teacher_report + father + teacher_report:father, desc = "y + father + y:father", k = 4),

  # health + father (no parent_report, but with different interactions from M1)
  list(f = r ~ teacher_report + health + father + health:father, desc = "y + health + father + h:f", k = 5),
  list(f = r ~ teacher_report + health + father, desc = "y + health + father", k = 4),

  # Three covariates + one covariate interaction
  list(f = r ~ teacher_report + health + father + parent_report + health:parent_report, desc = "y + h + f + pr + h:pr", k = 6),
  list(f = r ~ teacher_report + health + father + parent_report + father:parent_report, desc = "y + h + f + pr + f:pr", k = 6),
  list(f = r ~ teacher_report + health + father + parent_report + health:father, desc = "y + h + f + pr + h:f", k = 6)
)

cat(sprintf("%-35s %4s %15s %6s %6s %10s %8s\n",
    "Formula", "k", "PS range", ">0.99", "<0.01", "mu_ipw", "Status"))
cat(paste(rep("-", 95), collapse=""), "\n")

for (cand in candidates) {
  tryCatch({
    ps_spec <- list(formula.list = list(cand$f), h_alpha.list = list(h_alpha),
                    outcome = "teacher_report", inv_link = function(eta) 1/(1+exp(eta)))
    ebmr <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat, W)
    pf <- ebmr$ps_fit.list[[1]]
    ps <- pf$fitted.values
    mu <- mean(dat$r * dat$teacher_report / ps)
    degen <- sum(ps > 0.99) > 0 || sum(ps < 0.01) > 0
    # Check if same as M1 or M3
    same_m1 <- abs(mu - 0.2113) < 0.001
    same_m3 <- abs(mu - 0.2859) < 0.001 || abs(mu - 0.2865) < 0.001
    status <- if (degen) "DEGEN" else if (same_m1) "=M1" else if (same_m3) "=M3" else "OK"
    cat(sprintf("%-35s %4d [%.4f,%.4f] %6d %6d %10.4f %8s\n",
        cand$desc, ncol(pf$design_matrix), min(ps), max(ps),
        sum(ps > 0.99), sum(ps < 0.01), mu, status))
  }, error = function(e) {
    cat(sprintf("%-35s ERROR: %s\n", cand$desc, e$message))
  })
}
