setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Mental_Health_Data_Application")
devtools::load_all("../EBMRalgorithmFast4", quiet=TRUE)
source("MHD_functions.R")
dat <- gen_data(read.csv("data_application.csv"), round(2486*read.csv("data_application.csv")$Percentage/100), 2486)
dat$y <- dat$teacher_report
W <- function(g) solve(t(g)%*%g/nrow(g))
ps <- list(formula.list=list(
  r~teacher_report+health+father+teacher_report:health+teacher_report:father+father:health,
  r~teacher_report+parent_report+health+teacher_report:parent_report+teacher_report:health+parent_report:health,
  r~teacher_report+parent_report+father+teacher_report:parent_report+teacher_report:father+parent_report:father),
  h_alpha.list=list(c("health","father","parent_report","fp","fh","hp","fhp"),
                     c("health","father","parent_report","fp","fh","hp","fhp"),
                     c("health","father","parent_report","fp","fh","hp","fhp")),
  outcome="teacher_report", inv_link=function(eta)1/(1+exp(eta)), optimizer="L-BFGS-B")
ebmr <- EBMRAlgorithmFast4$new("teacher_report", ps, dat, W)
ps1 <- ebmr$ps_fit.list[[1]]$fitted.values
ps2 <- ebmr$ps_fit.list[[2]]$fitted.values
ps3 <- ebmr$ps_fit.list[[3]]$fitted.values

cat("=== PS collinearity ===\n\n")
cat("Correlation matrix:\n")
print(round(cor(cbind(pi1=ps1, pi2=ps2, pi3=ps3)), 4))

cat("\nDistinct PS values per model:\n")
cat("  M1:", length(unique(round(ps1,4))), "\n")
cat("  M2:", length(unique(round(ps2,4))), "\n")
cat("  M3:", length(unique(round(ps3,4))), "\n")

cat("\nPS summary:\n")
cat("  M1:", round(summary(ps1), 4), "\n")
cat("  M2:", round(summary(ps2), 4), "\n")
cat("  M3:", round(summary(ps3), 4), "\n")

# Check if PS columns are linearly dependent
ps_mat <- cbind(1, ps1, ps2, ps3)
cat("\nRank of (1, pi1, pi2, pi3):", qr(ps_mat)$rank, "out of 4\n")
sv <- svd(scale(cbind(ps1, ps2, ps3)))$d
cat("Singular values:", round(sv, 4), "\n")
cat("Condition number:", round(max(sv)/min(sv), 2), "\n")
