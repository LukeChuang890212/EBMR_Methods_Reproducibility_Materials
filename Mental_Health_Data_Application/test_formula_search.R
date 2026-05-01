## Search for non-degenerate M2 and M3 formulas
## Constraints:
##   M1: must include teacher_report + health (current, works fine)
##   M2: must include teacher_report + father (to distinguish from M1)
##   M3: must include teacher_report + parent_report (to distinguish from M1)
## Try all possible formulas with these required terms + optional extras
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
h_nu_fn <- function(data) {
  cbind(health = data$health, father = data$father, parent_report = data$parent_report,
        fp = data$fp, fh = data$fh, hp = data$hp, fhp = data$fhp)
}

## M1 is fixed (works well)
m1_formula <- r ~ teacher_report + health + father + teacher_report:health + teacher_report:father

## Generate candidate M2 formulas (must include teacher_report + father)
## Optional additions: health, parent_report, and interactions with teacher_report
m2_candidates <- list(
  # No extra covariates
  list(f = r ~ teacher_report + father, desc = "y + father"),
  # Add one covariate
  list(f = r ~ teacher_report + father + health, desc = "y + father + health"),
  list(f = r ~ teacher_report + father + parent_report, desc = "y + father + pr"),
  # With teacher_report interactions
  list(f = r ~ teacher_report + father + teacher_report:father, desc = "y + father + y:father"),
  list(f = r ~ teacher_report + father + health + teacher_report:father, desc = "y + father + health + y:father"),
  list(f = r ~ teacher_report + father + parent_report + teacher_report:father, desc = "y + father + pr + y:father"),
  # With covariate interactions (no y interactions)
  list(f = r ~ teacher_report + father + health + father:health, desc = "y + father + health + f:h"),
  list(f = r ~ teacher_report + father + parent_report + father:parent_report, desc = "y + father + pr + f:pr"),
  # Two covariates
  list(f = r ~ teacher_report + father + health + parent_report, desc = "y + father + health + pr"),
  # Two covariates + covariate interaction
  list(f = r ~ teacher_report + father + health + parent_report + father:health, desc = "y + father + health + pr + f:h"),
  list(f = r ~ teacher_report + father + health + parent_report + father:parent_report, desc = "y + father + health + pr + f:pr")
)

## Generate candidate M3 formulas (must include teacher_report + parent_report)
m3_candidates <- list(
  list(f = r ~ teacher_report + parent_report, desc = "y + pr"),
  list(f = r ~ teacher_report + parent_report + health, desc = "y + pr + health"),
  list(f = r ~ teacher_report + parent_report + father, desc = "y + pr + father"),
  list(f = r ~ teacher_report + parent_report + teacher_report:parent_report, desc = "y + pr + y:pr"),
  list(f = r ~ teacher_report + parent_report + health + teacher_report:parent_report, desc = "y + pr + health + y:pr"),
  list(f = r ~ teacher_report + parent_report + father + teacher_report:parent_report, desc = "y + pr + father + y:pr"),
  list(f = r ~ teacher_report + parent_report + health + parent_report:health, desc = "y + pr + health + pr:h"),
  list(f = r ~ teacher_report + parent_report + father + parent_report:father, desc = "y + pr + father + pr:f"),
  list(f = r ~ teacher_report + parent_report + health + father, desc = "y + pr + health + father"),
  list(f = r ~ teacher_report + parent_report + health + father + parent_report:health, desc = "y + pr + health + father + pr:h"),
  list(f = r ~ teacher_report + parent_report + health + father + parent_report:father, desc = "y + pr + health + father + pr:f")
)

cat("=== Searching for non-degenerate M2 formulas ===\n\n")
cat(sprintf("%-45s %4s %10s %6s %6s %10s\n", "Formula", "k", "PS range", ">0.99", "<0.01", "mu_ipw"))
cat(paste(rep("-", 90), collapse=""), "\n")

for (cand in m2_candidates) {
  tryCatch({
    ps_spec <- list(
      formula.list = list(cand$f),
      h_alpha.list = list(h_alpha),
      outcome = "teacher_report",
      inv_link = function(eta) 1 / (1 + exp(eta))
    )
    ebmr <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat, W)
    pf <- ebmr$ps_fit.list[[1]]
    ps <- pf$fitted.values
    k <- ncol(pf$design_matrix)
    mu <- mean(dat$r * dat$teacher_report / ps)
    degen <- sum(ps > 0.99) > 0 || sum(ps < 0.01) > 0
    flag <- if(degen) " *DEGEN*" else " OK"
    cat(sprintf("%-45s %4d [%.4f,%.4f] %6d %6d %10.4f%s\n",
        cand$desc, k, min(ps), max(ps), sum(ps > 0.99), sum(ps < 0.01), mu, flag))
  }, error = function(e) {
    cat(sprintf("%-45s ERROR: %s\n", cand$desc, e$message))
  })
}

cat("\n\n=== Searching for non-degenerate M3 formulas ===\n\n")
cat(sprintf("%-45s %4s %10s %6s %6s %10s\n", "Formula", "k", "PS range", ">0.99", "<0.01", "mu_ipw"))
cat(paste(rep("-", 90), collapse=""), "\n")

for (cand in m3_candidates) {
  tryCatch({
    ps_spec <- list(
      formula.list = list(cand$f),
      h_alpha.list = list(h_alpha),
      outcome = "teacher_report",
      inv_link = function(eta) 1 / (1 + exp(eta))
    )
    ebmr <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat, W)
    pf <- ebmr$ps_fit.list[[1]]
    ps <- pf$fitted.values
    k <- ncol(pf$design_matrix)
    mu <- mean(dat$r * dat$teacher_report / ps)
    degen <- sum(ps > 0.99) > 0 || sum(ps < 0.01) > 0
    flag <- if(degen) " *DEGEN*" else " OK"
    cat(sprintf("%-45s %4d [%.4f,%.4f] %6d %6d %10.4f%s\n",
        cand$desc, k, min(ps), max(ps), sum(ps > 0.99), sum(ps < 0.01), mu, flag))
  }, error = function(e) {
    cat(sprintf("%-45s ERROR: %s\n", cand$desc, e$message))
  })
}
