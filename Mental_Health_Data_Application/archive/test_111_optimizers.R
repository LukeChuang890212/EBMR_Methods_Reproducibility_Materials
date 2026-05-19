## Test different optimizers for 111 ensemble SE discrepancy
## k=6 formulas, h_alpha=h_nu=saturated
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Mental_Health_Data_Application")
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)
library(Matrix); library(numDeriv); library(nloptr)
source("MHD_functions.R")

original_data <- read.csv("data_application.csv")
percent <- original_data$Percentage
n <- 2486
class_n <- round(n * percent / 100)
dat <- gen_data(original_data, class_n, n)
dat$y <- dat$teacher_report

W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
compute_pi <- function(eta) plogis(-pmin(pmax(eta, -20), 20))

h_alpha_names <- c("health", "father", "parent_report", "fp", "fh", "hp", "fhp")
h_nu_fn <- function(data) cbind(health=data$health, father=data$father, parent_report=data$parent_report,
                                 fp=data$fp, fh=data$fh, hp=data$hp, fhp=data$fhp)

# Fast4 package spec
ps_spec <- list(
  formula.list = list(
    r ~ teacher_report + health + father + teacher_report:health + teacher_report:father,
    r ~ teacher_report + health + parent_report + teacher_report:health + teacher_report:parent_report,
    r ~ teacher_report + parent_report + father + teacher_report:parent_report + teacher_report:father
  ),
  h_alpha.list = list(h_alpha_names, h_alpha_names, h_alpha_names),
  outcome = "teacher_report", inv_link = function(eta) 1/(1+exp(eta))
)

B <- 300
cat("=== Testing optimizers for 111 ensemble (B=300) ===\n\n")

## Strategy 1: Fast4 L-BFGS-B (default)
cat("--- 1. Fast4 L-BFGS-B ---\n")
ps1 <- ps_spec; ps1$optimizer <- "L-BFGS-B"
ebmr <- EBMRAlgorithmFast4$new("teacher_report", ps1, dat, W_fn)
res <- ebmr$EBMR_IPW(h_nu = h_nu_fn, true_ps = NULL)
cat(sprintf("  Point: mu=%.4f, se=%.4f\n", res$mu_ipw, res$se_ipw))
mus1 <- rep(NA, B)
for (b in 1:B) {
  set.seed(12345+b); dat_b <- dat[sample(1:n,n,replace=TRUE),]
  tryCatch({
    eb <- EBMRAlgorithmFast4$new("teacher_report", ps1, dat_b, W_fn)
    r <- eb$EBMR_IPW(h_nu=h_nu_fn, true_ps=NULL, se.fit=FALSE)
    mus1[b] <- r$mu_ipw
  }, error=function(e){})
}
cat(sprintf("  Boot SE=%.4f, Ratio=%.3f\n\n", sd(mus1,na.rm=TRUE), res$se_ipw/sd(mus1,na.rm=TRUE)))

## Strategy 2: Fast4 constrained_nr
cat("--- 2. Fast4 constrained_nr ---\n")
ps2 <- ps_spec; ps2$optimizer <- "constrained_nr"
ebmr2 <- EBMRAlgorithmFast4$new("teacher_report", ps2, dat, W_fn)
res2 <- ebmr2$EBMR_IPW(h_nu = h_nu_fn, true_ps = NULL)
cat(sprintf("  Point: mu=%.4f, se=%.4f\n", res2$mu_ipw, res2$se_ipw))
mus2 <- rep(NA, B)
for (b in 1:B) {
  set.seed(12345+b); dat_b <- dat[sample(1:n,n,replace=TRUE),]
  tryCatch({
    eb <- EBMRAlgorithmFast4$new("teacher_report", ps2, dat_b, W_fn)
    r <- eb$EBMR_IPW(h_nu=h_nu_fn, true_ps=NULL, se.fit=FALSE)
    mus2[b] <- r$mu_ipw
  }, error=function(e){})
}
cat(sprintf("  Boot SE=%.4f, Ratio=%.3f\n\n", sd(mus2,na.rm=TRUE), res2$se_ipw/sd(mus2,na.rm=TRUE)))

## Strategy 3: Two-step GMM (only 2 W updates) via package with max_outer_override
## Can't easily control this in package, so use L-BFGS-B and check
## Actually, let's test with Fast2 package for comparison
cat("--- 3. Fast2 Gauss-Newton ---\n")
tryCatch({
  library(EBMRalgorithmFast2)
  ps_f2 <- list(
    formula.list = ps_spec$formula.list,
    h_x_names.list = list(h_alpha_names, h_alpha_names, h_alpha_names),
    outcome = "teacher_report", inv_link = function(eta) 1/(1+exp(eta))
  )
  h_nu_f2 <- function(data) cbind(health=data$health, father=data$father, parent_report=data$parent_report,
                                    fp=data$fp, fh=data$fh, hp=data$hp, fhp=data$fhp)
  ebmr3 <- EBMRalgorithmFast2$new("teacher_report", ps_f2, dat, W_fn)
  res3 <- ebmr3$EBMR_IPW(h_nu = h_nu_f2, true_ps = NULL)
  cat(sprintf("  Point: mu=%.4f, se=%.4f\n", res3$mu_ipw, res3$se_ipw))
  mus3 <- rep(NA, B)
  for (b in 1:B) {
    set.seed(12345+b); dat_b <- dat[sample(1:n,n,replace=TRUE),]
    tryCatch({
      eb <- EBMRalgorithmFast2$new("teacher_report", ps_f2, dat_b, W_fn)
      r <- eb$EBMR_IPW(h_nu=h_nu_f2, true_ps=NULL, se.fit=FALSE)
      mus3[b] <- r$mu_ipw
    }, error=function(e){})
  }
  cat(sprintf("  Boot SE=%.4f, Ratio=%.3f\n\n", sd(mus3,na.rm=TRUE), res3$se_ipw/sd(mus3,na.rm=TRUE)))
}, error=function(e) cat(sprintf("  ERROR: %s\n\n", e$message)))

## Strategy 4: nloptr L-BFGS inner loop (found better M2 solution earlier)
cat("--- 4. nloptr L-BFGS (inline GMM) ---\n")
build_dm <- function(d, m) {
  if (m==1) cbind(1, d$y, d$health, d$father, d$y*d$health, d$y*d$father)
  else if (m==2) cbind(1, d$y, d$health, d$parent_report, d$y*d$health, d$y*d$parent_report)
  else cbind(1, d$y, d$parent_report, d$father, d$y*d$parent_report, d$y*d$father)
}
build_hx <- function(d) {
  h1 <- d[h_alpha_names]; for(j in 1:ncol(h1)) h1[,j] <- as.factor(h1[,j])
  model.matrix(lm(rep(1,nrow(d))~., data=h1))
}

run_nloptr_gmm <- function(dm, hx, rv, nn) {
  k<-ncol(dm); hd<-ncol(hx)
  gf<-function(a)(rv/compute_pi(as.vector(dm%*%a))-1)*hx
  Gf<-function(a){g<-gf(a);matrix(.colMeans(g,nn,hd),hd,1)}
  Gammaf<-function(a){pv<-compute_pi(as.vector(dm%*%a));crossprod(hx*(rv*(1-pv)/pv),dm)/nn}
  of<-function(a,W){G<-Gf(a);as.numeric(crossprod(G,W%*%G))}
  grf<-function(a,W){2*as.vector(crossprod(Gammaf(a),W%*%Gf(a)))}
  W<-diag(hd)
  r<-nloptr(x0=rep(0,k),eval_f=function(x)of(x,W),eval_grad_f=function(x)grf(x,W),
            opts=list(algorithm="NLOPT_LD_LBFGS",maxeval=5000,xtol_rel=1e-12))
  est<-r$solution
  for(t in 1:200){
    gm<-gf(est);W<-tryCatch(W_fn(gm),error=function(e)diag(hd))
    if(max(abs(grf(est,W)))<1e-7) break
    r<-nloptr(x0=est,eval_f=function(x)of(x,W),eval_grad_f=function(x)grf(x,W),
              opts=list(algorithm="NLOPT_LD_LBFGS",maxeval=5000,xtol_rel=1e-12))
    est<-r$solution
  }
  compute_pi(as.vector(dm%*%est))
}

# Point estimate with nloptr
ps_all <- list()
for(m in 1:3) ps_all[[m]] <- run_nloptr_gmm(build_dm(dat,m), build_hx(dat), dat$r, n)
ps_mat <- do.call(cbind, ps_all)
# Use package for nu estimation
ebmr4 <- EBMRAlgorithmFast4$new("teacher_report", ps1, dat, W_fn)
# Override fitted values
for(m in 1:3) ebmr4$ps_fit.list[[m]]$fitted.values <- ps_all[[m]]
# Can't easily override — just report the inline result
mu4 <- mean(dat$r * dat$teacher_report / ps_all[[1]])  # M1 dominated
cat(sprintf("  M1 mu=%.4f\n", mu4))

mus4 <- rep(NA, B)
for(b in 1:B) {
  if(b%%100==0) cat(sprintf("  rep %d...\n",b))
  set.seed(12345+b); dat_b <- dat[sample(1:n,n,replace=TRUE),]
  tryCatch({
    hx_b <- build_hx(dat_b)
    ps_b <- list()
    for(m in 1:3) ps_b[[m]] <- run_nloptr_gmm(build_dm(dat_b,m), hx_b, dat_b$r, nrow(dat_b))
    # Simple: use M1 (assuming M1 dominates)
    mus4[b] <- mean(dat_b$r * dat_b$teacher_report / ps_b[[1]])
  }, error=function(e){})
}
cat(sprintf("  M1-only Boot SE=%.4f\n\n", sd(mus4,na.rm=TRUE)))

cat("=== Summary ===\n")
cat(sprintf("  %-25s %10s %10s %10s\n", "Optimizer", "Anal SE", "Boot SE", "Ratio"))
cat(paste(rep("-",60),collapse=""),"\n")
cat(sprintf("  %-25s %10.4f %10.4f %10.3f\n", "Fast4 L-BFGS-B", res$se_ipw, sd(mus1,na.rm=TRUE), res$se_ipw/sd(mus1,na.rm=TRUE)))
cat(sprintf("  %-25s %10.4f %10.4f %10.3f\n", "Fast4 constrained_nr", res2$se_ipw, sd(mus2,na.rm=TRUE), res2$se_ipw/sd(mus2,na.rm=TRUE)))
