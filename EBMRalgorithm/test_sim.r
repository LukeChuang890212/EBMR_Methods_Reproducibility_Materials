###################################
#  ---------- Test ----------------
###################################

devtools::install_github("LukeChuang890212/EBMR_Methods_Reproducibility_Materials/EBMRalgorithm")
library(EBMRalgorithm)

library(numDeriv)
library(tidyverse)
library(Matrix)

Chuang2023.1.1 = function(n){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2)  1/(1+exp(-0.2+0.1*y+0.1*u1+0.1*u2))

  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

data = Chuang2023.1.1(1000)
r = data$r
y = data$y
u1 = data$u1
u2 = data$u2
z1 = data$z1
z2 = data$z2

J = 3
dat = data
alpha.true = c(-0.2, 0.1, 0.1, 0.1)
true.pi.model = function(y, u1, u2, r, n, alpha.true) 1/(1+exp(cbind(rep(1, 1000), y, u1, u2)%*%alpha.true))
true.pi = true.pi.model(y, u1, u2, r, 1000, alpha.true)

source("G:\\Other computers\\我的電腦\\MNAR-Simulation\\MNAR_2023\\ChuangChao2023_SM.r")
source("G:\\Other computers\\我的電腦\\MNAR-Simulation\\MNAR_2023\\Wang2014_12.r")
n = sum(r)
N = 1000
h.list = list(function(u1, u2, z1, z2) list(cbind(as.factor(u1), as.factor(z1)), cbind(u2, z2)),
              function(u1, u2, z1, z2) list(cbind(as.factor(u1), as.factor(z1)), cbind(z2)),
              function(u1, u2, z1, z2) list(cbind(as.factor(z1)), cbind(u2, z2)))

propensity.list = list(list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = c("u1"),
                            model.x2.names =c("u2")),
                       list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = c("u1"),
                            model.x2.names = NULL),
                       list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = NULL,
                            model.x2.names = c("u2")))

start = Sys.time()
pi.fit.list = list()
for(j in 1:J){
  pi.fit.list[[j]] = Wang2014.1(auxilliary = h.list[[j]](u1, u2, z1, z2),
                                model.y = propensity.list[[j]]$model.y(y),
                                model.x1.names = propensity.list[[j]]$model.x1.names,
                                model.x2.names = propensity.list[[j]]$model.x2.names,
                                w = propensity.list[[j]]$w, w.prime = propensity.list[[j]]$w.prime)
}
pi.fit.list[[1]]$theta.hat
pi.fit.list[[1]]$se
auxilliary = list(cbind(as.factor(u1)), cbind(u2), y)
est.res1 = ChuangChao2023(pi.fit.list, auxilliary, family, ortho = TRUE, true.pi)
est1 =  unlist(est.res1[1:4])
est1
Sys.time() - start

# source("C:\\Users\\stat-pc\\Desktop\\NTHU_Research\\EBMR_Methods_Reproducibility_Materials\\EBMRalgorithm\\R\\Methods.r")
ps_specifications = list(
  formula.list = list(
    r ~ o(y) + u1 + u2,
    r ~ o(y) + u1,
    r ~ o(y) + u2
  ),
  h_x_names.list = list(
    c("u1", "u2", "z1", "z2"),
    c("u1", "z1", "z2"),
    c("u2", "z1", "z2")
  ),
  inv_link = function(eta) 1/(1+exp(eta))
)

W = function(g.matrix){
  return(solve(t(g.matrix)%*%g.matrix/nrow(g.matrix)))
}

start = Sys.time()
ebmr <- EBMRAlgorithm$new("y", ps_specifications, data, W)
ebmr$ps_fit.list[[1]]$coefficients
ebmr$ps_fit.list[[1]]$se
result = ebmr$EBMR_IPW(h_x_names = c("u1", "u2"), true_ps = true.pi)
result = unlist(result[1:4])
result
Sys.time() - start

# ps_fit.list = list()
# for(j in 1:J){
#   formula = ps_specifications$formula.list[[j]]
#   h_x_names = ps_specifications$h_x_names.list[[j]]
#   inv_link = ps_specifications$inv_link
#   if(is.null(wt)){
#     ps_fit.list[[j]] = WangShaoKim2014(formula, h_x_names, inv_link)
#   }else{
#     ps_fit.list[[j]] = WangShaoKim2014(formula, h_x_names, inv_link, wt, se.fit = F)
#   }
# }
