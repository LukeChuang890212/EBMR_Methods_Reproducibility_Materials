############################################
# Preparation for computing the integration
############################################
H = function(k, x){
  summands = list()
  summands[[1]] = c(-2)
  summands[[2]] = c(-2, 4)
  for(i in 3:k){
    last.coef = summands[[i-1]]
    if(i%%2 == 0){
      last.order = seq(1, i-1, 2)
      current.coef = rep(0, i%/%2+1)
      current.coef[1] = last.coef[1]
      for(j in 2:(i%/%2)){
        current.coef[j] = last.coef[j-1]*(-2)+last.coef[j]*last.order[j]
      }
      current.coef[length(current.coef)] = last.coef[length(last.coef)]*(-2)
      summands[[i]] = current.coef
    }else{
      last.order = seq(0, i-1, 2)
      current.coef = rep(0, i%/%2+1)
      for(j in 1:(i%/%2)){
        current.coef[j] = last.coef[j]*(-2)+last.coef[j+1]*last.order[j+1]
      }
      current.coef[length(current.coef)] = last.coef[length(last.coef)]*(-2)
      summands[[i]] = current.coef
    }
  }
  coef = summands[[k]]
  current.order = NA
  if(k%%2 == 0){
    current.order = seq(0, k, 2)
  }else{
    current.order = seq(1, k, 2)
  }

  return(((-1)^k)*sum(coef*(x^current.order)))
}

# Hk = Vectorize(function(x) H(k = 15, x))
# t = uniroot.all(Hk, lower = -5, upper = 5, n = 19)
#
# Hk = Vectorize(function(x) H(k = 14, x))
# w.i = ((2^(15-1))*factorial(15)*sqrt(pi))/((15^2)*(Hk(t))^2)

############################################################
# All functions needed for the case: pi(y, u)
############################################################

pi.s = function(y, u) exp(1+y+u)/(1+exp(1+y+u))

p = function(x, beta) exp(x%*%beta)/(1+exp(x%*%beta))

phi = function(y, x, beta) (p(x, beta)^y)*(1-p(x, beta))^(1-y)

# The pseudo-1st-order derivative for the density of y given x with respect to beta
phi.prime = function(y, x, beta) ((-1)^as.numeric(y == 0))*p(x, beta)*(1-p(x, beta))

# The pseudo-1st-order derivative for the density of y given x with respect to beta
phi.2 = function(y, x, beta) ((-1)^as.numeric(y == 0))*p(x, beta)*(1-p(x, beta))*(1-2*p(x, beta))

f.beta = function(y, x, beta){
  phi.prime.i = phi.prime(y, x, beta)
  apply(x, 2, function(xj) xj%*%phi.prime.i)
}

S.eff.g = function(beta, y, z, u, r, K){
  N = length(r); n = sum(r)
  if(is.null(u)) u = rep(0, length(y))
  x = cbind(1, z, u)
  x.dim = ncol(x)

  w = phi.prime(y, x, beta[1:x.dim])/phi(y, x, beta[1:x.dim])

  B.i = matrix(NA, n, x.dim)
  complete.index = which(r == 1)
  for(i in 1:n){
    w.num = phi.prime(y[complete.index[i]], x, beta[1:x.dim])*K(u-u[complete.index[i]])/N
    num = apply(x, 2, function(xj) xj%*%w.num)
    denom = mean(phi(y[complete.index[i]], x, beta[1:x.dim])*K(u-u[complete.index[i]]))
    B.i[i,] = num/denom
  }
  return(apply(x, 2, function(xj) xj%*%(r*w))-apply(B.i, 2, sum))
}

# b = function(beta, y.i, z, u, u.i, K){
#   x = cbind(1, z, u)
#
#   pseudoInt_f.beta_pi.s.i = phi.prime(0, x, beta)*pi.s(0, u[i])+phi.prime(1, x, beta)*pi.s(1, u[i])
#   int_f.Y_X_pi.s.i = phi(0, x, beta)*pi.s(0, u[i])+phi(1, x, beta)*pi.s(1, u[i])
#
#   phi.i = phi(y.i, x, beta)
#   phi.prime.i = phi.prime(y.i, x, beta)
#   K.i = K(u-u.i)
#
#   w.num = (phi.prime.i+pseudoInt_f.beta_pi.s.i/(1-int_f.Y_X_pi.s.i)*phi.i)*K.i/N
#   num = apply(x, 2, function(xj) xj%*%w.num)
#   denom = mean(phi.i*K.i)
#
#   g.i = num/denom
#
#   num = (mean(phi(0, x, beta)/(1-int_f.Y_X_pi.s.i)*phi.i*K.i)*b(beta, 0, z, u, u.i, K)
#                      + mean(phi(1, x, beta)/(1-int_f.Y_X_pi.s.i)*phi.i*K.i)*b(beta, 1, z, u, u.i, K))
#   int_k_b.i = num/denom
#
#   return(g.i+int_k_b.i)
# }

b = function(beta, supp.yu, z, u, K, b.m){
  if(is.vector(z)) z = as.matrix(z)
  if(is.vector(u)) u = as.matrix(u)
  supp.len = nrow(supp.yu)
  x.dim = ifelse(is.null(u), 1+ncol(z), 1+ncol(z)+ncol(u))

  g.i = int_k_b.i = matrix(NA, supp.len, x.dim)
  for(i in 1:supp.len){
    y.i = supp.yu[i, 1]; u.i = supp.yu[i, 2]
    x = cbind(1, z, u.i)

    pseudoInt_f.beta_pi.s.i = phi.prime(0, x, beta)*pi.s(0, u.i)+phi.prime(1, x, beta)*pi.s(1, u.i)
    int_f.Y_X_pi.s.i = phi(0, x, beta)*pi.s(0, u.i)+phi(1, x, beta)*pi.s(1, u.i)

    phi.i = phi(y.i, x, beta)
    phi.prime.i = phi.prime(y.i, x, beta)
    K.i = K(u-u.i)

    w.num = (phi.prime.i+pseudoInt_f.beta_pi.s.i/(1-int_f.Y_X_pi.s.i)*phi.i)*K.i/N
    num = apply(x, 2, function(xj) xj%*%w.num)
    denom = mean(phi.i*K.i)

    g.i[i,] = num/denom

    num = (mean(phi(0, x, beta)/(1-int_f.Y_X_pi.s.i)*phi.i*K.i)*b.m[which(supp.yu[, 1] == 0 & supp.yu[, 2] == u.i),]
           + mean(phi(1, x, beta)/(1-int_f.Y_X_pi.s.i)*phi.i*K.i)*b.m[which(supp.yu[, 1] == 1 & supp.yu[, 2] == u.i),])

    int_k_b.i[i,] = num/denom
  }

  return(g.i-int_k_b.i)
}

# supp.yu = expand.grid(sort(unique(y[r == 1])), sort(unique(u[r == 1])))
#
# conv.err = 10^8
# t = 1
#
# b.arr = array(0, dim = c(nrow(supp.yu), beta.dim, 1))
# while(conv.err > 10^(-4)){
#   b.hat = b(beta, supp.yu, z, u, K, b.arr[,,t])
#
#   b.arr = abind(b.arr, b.hat)
#   t = t + 1
#
#   conv.err = max(abs(b.arr[,,t]-b.arr[,,t-1]))
#   print(c(t, conv.err))
# }

S.eff.g.prime = function(beta, y, z, u, r, K){
  if(is.null(u)) u = rep(0, length(y))
  N = length(r); n = sum(r)
  x = cbind(1, z, u)
  x.dim = ncol(x)

  phi.i = phi(y, x, beta[1:x.dim])
  phi.prime.i = phi.prime(y, x, beta[1:x.dim])
  phi.2.i = phi.2(y, x, beta[1:x.dim])
  w = (phi.2.i/phi.i)-(phi.prime.i/phi.i)^2

  xx.i = array(NA, dim = c(x.dim, x.dim, N))
  for(i in 1:N){
    xx.i[,,i] = x[i,]%*%t(x[i,])
  }

  # A_beta: the partial derivative of A_i with respect to beta
  A_beta = apply(xx.i, c(1, 2), function(xij) xij%*%(w*r))

  B.i_beta = array(NA, dim = c(x.dim, x.dim, n))
  complete.index = which(r == 1)
  for(i in 1:length(complete.index)){
    phi.i = phi(y[complete.index[i]], x, beta[1:x.dim])
    phi.prime.i = phi.prime(y[complete.index[i]], x, beta[1:x.dim])
    phi.2.i = phi.2(y[complete.index[i]], x, beta[1:x.dim])

    K.i = 1
    if(any(u != 0)) K.i = K(u-u[complete.index[i]])

    weighted.xx.i = apply(xx.i, c(1, 2), function(xij) xij%*%(phi.2.i*K.i)/N)
    weighted.x.i = apply(x, 2, function(xj) xj%*%(phi.prime.i*K.i)/N)

    B.i_beta[,,i] = weighted.xx.i/mean(phi.i*K.i)-(weighted.x.i%*%t(weighted.x.i))/(mean(phi.i*K.i)^2)
  }

  # B_beta: the partial derivative of B with respect to beta
  B_beta = apply(B.i_beta, c(1, 2), sum)

  return(A_beta-B_beta)
  # return(-A_beta/sigma)
}

##################################################################
# The main function for estimation: pi(y, u)
##################################################################

ZhaoMa2022 = function(y, z, u, r, init, K){
  if(is.vector(z)) z = as.matrix(z)
  if(is.vector(u)) u = as.matrix(u)
  x.dim = ifelse(is.null(u), 1+ncol(z), 1+ncol(z)+ncol(u))
  beta.dim = length(init)

  # start = Sys.time()
  conv.err = 10^8
  t = 1

  beta.m = matrix(init, nrow = beta.dim)
  while(conv.err > 10^(-3)){
    beta.hat = rep(0, beta.dim)

    first.deriv_beta = S.eff.g(beta.m[,t], y, z, u, r, K)[1:x.dim]
    second.deriv_beta = S.eff.g.prime(beta.m[,t], y, z, u, r, K)[1:x.dim, 1:x.dim]

    beta.hat[1:x.dim] = beta.m[1:x.dim, t]-solve(second.deriv_beta)%*%first.deriv_beta

    beta.m = cbind(beta.m, beta.hat)
    t = t + 1

    conv.err = max(abs(beta.m[,t]-beta.m[,t-1]))
    print(c(t, conv.err))
  }
  beta.hat = beta.m[,t]

  # print(Sys.time()-start)
  return(list(sol.path = beta.m, beta.hat = beta.hat))
}

# N = 500
# z = rnorm(N, 0.5, 0.5)
# beta = c(0.25, -0.5)
# y = cbind(1, z)%*%beta+rnorm(N)
# r = rbinom(N, size = 1, prob = exp(1+y)/(1+exp(1+y)))
# dat = data.frame(y = y, z = z, r = r)

##############################
# Test: coordinate descent
##############################

# obj = function(theta, y, z, u, r){
#   n = length(y)
#   value = sum(S.eff(theta, y, z, u, r)^2)
#   if(is.infinite(value)) print(value)
#   return(ifelse(is.infinite(value), 10^8, value))
# }
#
# cc.fit = lm(y~z, dat, subset = r == 1)
# theta.init = c(cc.fit$coef, 0, summary(cc.fit)$sigma)
# # theta.init = c(20, 10, 0, 4)
#
# start = Sys.time()
# conv.err = 10^8
# t = 1
#
# theta.m = matrix(theta.init, nrow = length(theta.init))
# while(conv.err > 10^(-3)){
#   for(i in 1:nrow(theta.m)){
#     obj.i = function(theta.i, y, z, u, r){
#       theta = theta.m[, t]
#       theta[i] = theta.i
#       value = obj(theta, y, z, u, r)
#       return(ifelse(is.infinite(value), 10^8, value))
#     }
#     opt = optim(theta.m[i, t],
#                 fn = obj.i,
#                 y = y, z = z, u = u, r = r,
#                 method = "L-BFGS-B")
#     theta.hat = theta.m[, t]
#     theta.hat[i] = opt$par
#     theta.m = cbind(theta.m, theta.hat)
#     t = t + 1
#   }
#   conv.err = max(abs(theta.m[,t]-theta.m[,t-1]))
#   print(c(t, conv.err))
# }
# print(Sys.time()-start)
# theta.m[,t]



