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
for(i in 1:3) {
  pf <- ebmr$ps_fit.list[[i]]
  cat(sprintf("M%d: obj=%.4e, grad=%.4e, PS=[%.4f,%.4f], >0.99:%d, <0.01:%d\n",
      i, pf$gmm_fit$opt$objective, pf$gmm_fit$opt$final_grad_norm,
      min(pf$fitted.values), max(pf$fitted.values),
      sum(pf$fitted.values>0.99), sum(pf$fitted.values<0.01)))
}
