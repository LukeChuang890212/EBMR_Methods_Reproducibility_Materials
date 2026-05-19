## Quick test: M2 = r ~ teacher_report + parent_report + health + teacher_report:health
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
h_nu <- function(data) cbind(health=data$health, father=data$father, parent_report=data$parent_report,
                              fp=data$fp, fh=data$fh, hp=data$hp, fhp=data$fhp)

ps_spec <- list(
  formula.list = list(r ~ teacher_report + parent_report + health + teacher_report:parent_report),
  h_alpha.list = list(h_alpha),
  outcome = "teacher_report",
  inv_link = function(eta) 1 / (1 + exp(eta))
)

ebmr <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat, W)
pf <- ebmr$ps_fit.list[[1]]
ps <- pf$fitted.values
res <- ebmr$EBMR_IPW(h_nu = h_nu, true_ps = NULL)

cat(sprintf("M2: r ~ teacher_report + parent_report + health + teacher_report:health\n"))
cat(sprintf("k=%d, h_dim=%d, overid=%d\n", ncol(pf$design_matrix), ncol(pf$h_x),
    ncol(pf$h_x) - ncol(pf$design_matrix)))
cat(sprintf("alpha = (%s)\n", paste(round(pf$coefficients, 4), collapse=", ")))
cat(sprintf("SE    = (%s)\n", paste(round(pf$se, 4), collapse=", ")))
cat(sprintf("PS range: [%.6f, %.6f], >0.99:%d, <0.01:%d\n",
    min(ps), max(ps), sum(ps > 0.99), sum(ps < 0.01)))
cat(sprintf("mu_ipw = %.4f, se_ipw = %.4f\n", res$mu_ipw, res$se_ipw))
