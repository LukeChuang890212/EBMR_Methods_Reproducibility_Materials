Wang2014.0 = function(auxilliary, model.y, model.x1.names, model.x2.names, w, w.prime){
  propensity = function(y, x, theta, L) 1/w(theta, y, x, L)

  r = dat$r
  x1 = auxilliary[[1]]; x2 = auxilliary[[2]]; auxilliary.y = auxilliary[[3]];
  if(!is.null(x1)) x1 = as.matrix(x1)
  if(!is.null(x2)) x2 = as.matrix(x2)
  model.x1 = NULL; model.x2 = NULL;
  if(!is.null(model.x1.names)) model.x1 = as.matrix(dat[model.x1.names])
  if(!is.null(model.x2.names)) model.x2 = as.matrix(dat[model.x2.names])

  N = length(r)
  discrete.dim = 0
  if(!is.null(x1)) for(j in 1:ncol(x1)) discrete.dim = discrete.dim + length(unique(x1[,j]))
  discrete.dim = ifelse(discrete.dim == 0, 0, discrete.dim - ncol(x1) + 1) # to avoid collinearity
  continous.dim = ifelse(!is.null(x2), ncol(x2), 0)
  auxilliary.dim = discrete.dim + continous.dim + 1

  mnar = ifelse(is.null(model.y), 0, 1)
  model.discrete.dim = ifelse(is.null(model.x1), 0, ncol(model.x1))
  model.continous.dim = ifelse(is.null(model.x2), 0, ncol(model.x2))
  param.num = c(model.discrete.dim, model.continous.dim)

  M = 1+mnar+sum(param.num)+1
  init = rep(0.5, M)

  d = NULL
  if(!is.null(x1)){
    x1 =  as.data.frame(x1)
    for(j in 1:ncol(x1)) x1[, j] = as.factor(x1[, j])
    d = model.matrix(lm(rep(1, N)~., data =  x1))
  }


  g = function(theta){
    g.m = matrix(NA, auxilliary.dim, N)
    rw = r*w(theta[1:(M-1)], model.y, cbind(model.x1, model.x2), N)
    if(discrete.dim > 0){
      for(l in 1:discrete.dim){g.m[l,] = t(d)[l,]*(rw-1)}
    }
    if(continous.dim > 0){
      for(l in (discrete.dim+1):(discrete.dim+continous.dim)) g.m[l,] = t(as.matrix(x2))[l-discrete.dim,]*(rw-1)
    }
    g.m[auxilliary.dim,] = y*rw-theta[M]
    return(g.m)
  }

  G = function(g.m){
    return(matrix(apply(g.m, 1, mean), auxilliary.dim, 1))
  }

  W = function(g.m){
    return(solve(g.m%*%t(g.m)/N))
  }

  Gamma = function(theta){
    gamma = array(NA, dim = c(auxilliary.dim, M, N))
    r_w.prime = r*w.prime(theta[1:(M-1)], model.y, cbind(model.x1, model.x2), N)
    if(discrete.dim > 0){
      for(l in 1:discrete.dim){
        gamma[l,1,] = t(d)[l,]*r_w.prime
        if(mnar == 1) gamma[l,2,] = t(d)[l,]*r_w.prime*model.y
        if(model.discrete.dim > 0){
          for(c in (1+mnar+1):(1+mnar+param.num[1])){
            gamma[l,c,] = t(d)[l,]*r_w.prime*t(as.matrix(model.x1))[c-(1+mnar),]
          }
        }
        if(model.continous.dim > 0){
          for(c in (1+mnar+param.num[1]+1):(1+mnar+sum(param.num))){
            gamma[l,c,] = t(d)[l,]*r_w.prime*t(as.matrix(model.x2))[c-(1+mnar+param.num[1]),]
          }
        }
      }
    }
    if(continous.dim > 0){
      for(l in (discrete.dim+1):(discrete.dim+continous.dim)){
        gamma[l,1,] = t(as.matrix(x2))[l-discrete.dim,]*r_w.prime
        if(mnar == 1) gamma[l,2,] = t(as.matrix(x2))[l-discrete.dim,]*r_w.prime*model.y
        if(model.discrete.dim >0){
          for(c in (1+mnar+1):(1+mnar+param.num[1])){
            gamma[l,c,] = t(as.matrix(x2))[l-discrete.dim,]*r_w.prime*t(as.matrix(model.x1))[c-(1+mnar),]
          }
        }
        if(model.continous.dim > 0){
          for(c in (1+mnar+param.num[1]+1):(1+mnar+sum(param.num))){
            gamma[l,c,] = t(as.matrix(x2))[l-discrete.dim,]*r_w.prime*t(as.matrix(model.x2))[c-(1+mnar+param.num[1]),]
          }
        }
      }
    }
    gamma[auxilliary.dim, 1:(1+mnar), ] = t(cbind(rep(1, length(auxilliary.y)), model.y)*as.vector(r_w.prime*auxilliary.y))
    if(length(model.x1.names)+length(model.x2.names) > 0){
      gamma[auxilliary.dim, (1+mnar+1):(M-1), ] = t(cbind(model.x1, model.x2)*as.vector(r_w.prime*auxilliary.y))
    }
    gamma[-auxilliary.dim, M, ] = 0
    gamma[auxilliary.dim, M, ] = -1
    return(apply(gamma, c(1, 2), mean))
  }

  # Gamma(theta)
  #
  # Gamma(theta){
  #   Gamma.hat = array(NA, dim = c(auxilliary.dim, N, M))
  #   for(l in 1:auxilliary.dim){
  #     Gamma.hat[l,,] = jacobian(function(theta.v) g(theta.v)[l,], theta)
  #   }
  #   return(apply(Gamma.hat, c(1, 3), mean))
  # }

  obj = function(theta){
    g.m = g(theta)
    G.hat = G(g.m)
    value = t(G.hat)%*%G.hat
    return(ifelse(is.infinite(value) || is.na(value), 10^8, value))
  }

  theta.m = matrix(init, M)
  conv.err = 10^8
  t = 1

  while (conv.err > 10^(-8) & t < 1000){
    opt = optim(theta.m[,t], obj, method = "L-BFGS-B")
    theta.m = cbind(theta.m, opt$par)
    g.m = g(theta.m[,t+1]); W.hat = W(g.m);
    obj = function(theta){
      g.m = g(theta); G.hat = G(g.m);
      value = t(G.hat)%*%W.hat%*%G.hat
      return(ifelse(is.infinite(value) || is.na(value), 10^8, value))
    }
    conv.err = max(abs(theta.m[,t+1]-theta.m[,t]))
    t = t + 1
  }

  theta.hat = theta.m[, t]
  Gamma.hat = Gamma(theta.hat)
  g.m = g(theta.hat); W.hat = W(g.m);
  cov.hat = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)/N
  se = sqrt(diag(cov.hat))

  return(list(sol.path = theta.m, theta.hat = theta.hat, cov.hat = cov.hat,
              se = se, lower = theta.hat-qnorm(0.975)*se, upper = theta.hat+qnorm(0.975)*se,
              g.m = g.m, GammaW = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)%*%t(Gamma.hat)%*%W.hat,
              model = propensity, w.prime = w.prime, model.x.names = c(model.x1.names, model.x2.names), h = t(cbind(d, x2))))
}

Morikawa2021 = function(auxilliary, model.y, model.x1.names, model.x2.names, w, w.prime){
  propensity = function(y, x, theta, L) 1/w(theta, y, x, L)

  r = dat$r
  x1 = auxilliary[[1]]; x2 = auxilliary[[2]]; auxilliary.y = auxilliary[[3]];
  if(!is.null(x1)) x1 = as.matrix(x1)
  if(!is.null(x2)) x2 = as.matrix(x2)

  weighted.m = x2[, ncol(x2)]

  model.x1 = NULL; model.x2 = NULL;
  if(!is.null(model.x1.names)) model.x1 = as.matrix(dat[model.x1.names])
  if(!is.null(model.x2.names)) model.x2 = as.matrix(dat[model.x2.names])

  N = length(r)
  discrete.dim = 0
  if(!is.null(x1)) for(j in 1:ncol(x1)) discrete.dim = discrete.dim + length(unique(x1[,j]))
  discrete.dim = ifelse(discrete.dim == 0, 0, discrete.dim - ncol(x1) + 1) # to avoid collinearity
  continous.dim = ifelse(!is.null(x2), ncol(x2), 0)
  auxilliary.dim = discrete.dim + continous.dim + 1

  mnar = ifelse(is.null(model.y), 0, 1)
  model.discrete.dim = ifelse(is.null(model.x1), 0, ncol(model.x1))
  model.continous.dim = ifelse(is.null(model.x2), 0, ncol(model.x2))
  param.num = c(model.discrete.dim, model.continous.dim)

  M = 1+mnar+sum(param.num)+1
  init = rep(0.5, M)

  d = NULL
  if(!is.null(x1)){
    x1 =  as.data.frame(x1)
    for(j in 1:ncol(x1)) x1[, j] = as.factor(x1[, j])
    d = model.matrix(lm(rep(1, N)~., data =  x1))
  }

  g = function(theta){
    g.m = matrix(NA, auxilliary.dim, N)
    rw = r*w(theta[1:(M-1)], model.y, cbind(model.x1, model.x2), N)
    if(discrete.dim > 0){
      for(l in 1:discrete.dim){g.m[l,] = t(d)[l,]*(rw-1)}
    }
    if(continous.dim > 0){
      for(l in (discrete.dim+1):(discrete.dim+continous.dim)) g.m[l,] = t(as.matrix(x2))[l-discrete.dim,]*(rw-1)
    }
    g.m[auxilliary.dim,] = rw*(y-theta[M])-(rw-1)*weighted.m
    return(g.m)
  }

  G = function(g.m){
    return(matrix(apply(g.m, 1, mean), auxilliary.dim, 1))
  }

  W = function(g.m){
    return(solve(g.m%*%t(g.m)/N))
  }

  # Gamma = function(theta){
  #   gamma = array(NA, dim = c(auxilliary.dim, M, N))
  #   r_w.prime = r*w.prime(theta[1:(M-1)], model.y, cbind(model.x1, model.x2), N)
  #   if(discrete.dim > 0){
  #     for(l in 1:discrete.dim){
  #       gamma[l,1,] = t(d)[l,]*r_w.prime
  #       if(mnar == 1) gamma[l,2,] = t(d)[l,]*r_w.prime*model.y
  #       if(model.discrete.dim > 0){
  #         for(c in (1+mnar+1):(1+mnar+param.num[1])){
  #           gamma[l,c,] = t(d)[l,]*r_w.prime*t(as.matrix(model.x1))[c-(1+mnar),]
  #         }
  #       }
  #       if(model.continous.dim > 0){
  #         for(c in (1+mnar+param.num[1]+1):(1+mnar+sum(param.num))){
  #           gamma[l,c,] = t(d)[l,]*r_w.prime*t(as.matrix(model.x2))[c-(1+mnar+param.num[1]),]
  #         }
  #       }
  #     }
  #   }
  #   if(continous.dim > 0){
  #     for(l in (discrete.dim+1):(discrete.dim+continous.dim)){
  #       gamma[l,1,] = t(as.matrix(x2))[l-discrete.dim,]*r_w.prime
  #       if(mnar == 1) gamma[l,2,] = t(as.matrix(x2))[l-discrete.dim,]*r_w.prime*model.y
  #       if(model.discrete.dim >0){
  #         for(c in (1+mnar+1):(1+mnar+param.num[1])){
  #           gamma[l,c,] = t(as.matrix(x2))[l-discrete.dim,]*r_w.prime*t(as.matrix(model.x1))[c-(1+mnar),]
  #         }
  #       }
  #       if(model.continous.dim > 0){
  #         for(c in (1+mnar+param.num[1]+1):(1+mnar+sum(param.num))){
  #           gamma[l,c,] = t(as.matrix(x2))[l-discrete.dim,]*r_w.prime*t(as.matrix(model.x2))[c-(1+mnar+param.num[1]),]
  #         }
  #       }
  #     }
  #   }
  #   gamma[auxilliary.dim, 1:(1+mnar), ] = t(cbind(rep(1, length(auxilliary.y)), model.y)*as.vector(r_w.prime*auxilliary.y))
  #   if(length(model.x1.names)+length(model.x2.names) > 0){
  #     gamma[auxilliary.dim, (1+mnar+1):(M-1), ] = t(cbind(model.x1, model.x2)*as.vector(r_w.prime*auxilliary.y))
  #   }
  #   gamma[-auxilliary.dim, M, ] = 0
  #   gamma[auxilliary.dim, M, ] = -1
  #   return(apply(gamma, c(1, 2), mean))
  # }

  Gamma = function(theta){
    Gamma.hat = array(NA, dim = c(auxilliary.dim, N, M))
    for(l in 1:auxilliary.dim){
      Gamma.hat[l,,] = jacobian(function(theta.v) g(theta.v)[l,], theta)
    }
    return(apply(Gamma.hat, c(1, 3), mean))
  }

  obj = function(theta){
    g.m = g(theta)
    G.hat = G(g.m)
    value = t(G.hat)%*%G.hat
    return(ifelse(is.infinite(value) || is.na(value), 10^8, value))
  }

  theta.m = matrix(init, M)
  conv.err = 10^8
  t = 1

  while (conv.err > 10^(-8) & t < 1000){
    opt = optim(theta.m[,t], obj, method = "L-BFGS-B")
    theta.m = cbind(theta.m, opt$par)
    g.m = g(theta.m[,t+1]); W.hat = W(g.m);
    obj = function(theta){
      g.m = g(theta); G.hat = G(g.m);
      value = t(G.hat)%*%W.hat%*%G.hat
      return(ifelse(is.infinite(value) || is.na(value), 10^8, value))
    }
    conv.err = max(abs(theta.m[,t+1]-theta.m[,t]))
    t = t + 1
  }

  theta.hat = theta.m[, t]
  Gamma.hat = Gamma(theta.hat)
  g.m = g(theta.hat); W.hat = W(g.m);
  cov.hat = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)/N
  se = sqrt(diag(cov.hat))

  return(list(sol.path = theta.m, theta.hat = theta.hat, cov.hat = cov.hat,
              se = se, lower = theta.hat-qnorm(0.975)*se, upper = theta.hat+qnorm(0.975)*se,
              g.m = g.m, GammaW = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)%*%t(Gamma.hat)%*%W.hat,
              model = propensity, w.prime = w.prime, model.x.names = c(model.x1.names, model.x2.names), h = t(cbind(d, x2))))
}

Wang2014.1 = function(auxilliary = list(x1, x2), model.y, model.x1.names, model.x2.names, w, w.prime, init = NULL){
  propensity = function(y, x, theta, L) 1/w(theta, y, x, L)
  r = dat$r

  x1 = auxilliary[[1]]; x2 = auxilliary[[2]];
  if(!is.null(x1)) x1 = as.matrix(x1)
  if(!is.null(x2)) x2 = as.matrix(x2)

  model.x1 = NULL; model.x2 = NULL;
  if(!is.null(model.x1.names)) model.x1 = as.matrix(dat[model.x1.names])
  if(!is.null(model.x2.names)) model.x2 = as.matrix(dat[model.x2.names])

  N = length(r)
  discrete.dim = 0
  if(!is.null(x1)) for(j in 1:ncol(x1)) discrete.dim = discrete.dim + length(unique(x1[,j]))
  discrete.dim = ifelse(discrete.dim == 0, 0, discrete.dim - ncol(x1) + 1) # to avoid collinearity
  continous.dim = ifelse(!is.null(x2), ncol(x2), 0)
  auxilliary.dim = discrete.dim + continous.dim

  mnar = ifelse(is.null(model.y), 0, 1)
  model.discrete.dim = ifelse(is.null(model.x1), 0, ncol(model.x1))
  model.continous.dim = ifelse(is.null(model.x2), 0, ncol(model.x2))
  param.num = c(model.discrete.dim, model.continous.dim)

  M = 1+mnar+length(model.x1.names) + length(model.x2.names)
  if(is.null(init)) init = rep(0, M)

  d = NULL
  if(!is.null(x1)){
    x1 =  as.data.frame(x1)
    for(j in 1:ncol(x1)) x1[, j] = as.factor(x1[, j])
    d = model.matrix(lm(rep(1, N)~., data =  x1))
  }
  # x1 = as.numeric(as.matrix(x1, ncol(x1)))

  g = function(theta){
    g.m = matrix(NA, auxilliary.dim, N)
    rw = r*w(theta, model.y, cbind(model.x1, model.x2), N)
    if(discrete.dim > 0){
      for(l in 1:discrete.dim){g.m[l,] = t(d)[l,]*(rw-1)}
    }
    if(continous.dim > 0){
      for(l in (discrete.dim+1):(discrete.dim+continous.dim)) g.m[l,] = t(as.matrix(x2))[l-discrete.dim,]*(rw-1)
    }
    return(g.m)
  }

  G = function(g.m){
    return(matrix(apply(g.m, 1, mean), auxilliary.dim, 1))
  }

  W = function(g.m){
    return(solve(g.m%*%t(g.m)/N))
  }

  Gamma = function(theta){
    gamma = array(NA, dim = c(auxilliary.dim, M, N))
    r_w.prime = r*w.prime(theta, model.y, cbind(model.x1, model.x2), N)
    if(discrete.dim > 0){
      for(l in 1:discrete.dim){
        gamma[l,1,] = t(d)[l,]*r_w.prime
        if(mnar == 1) gamma[l,2,] = t(d)[l,]*r_w.prime*model.y
        if(model.discrete.dim > 0){
          for(c in (1+mnar+1):(1+mnar+param.num[1])){
            gamma[l,c,] = t(d)[l,]*r_w.prime*t(as.matrix(model.x1))[c-(1+mnar),]
          }
        }
        if(model.continous.dim > 0){
          for(c in (1+mnar+param.num[1]+1):(1+mnar+sum(param.num))){
            gamma[l,c,] = t(d)[l,]*r_w.prime*t(as.matrix(model.x2))[c-(1+mnar+param.num[1]),]
          }
        }
      }
    }
    if(continous.dim > 0){
      for(l in (discrete.dim+1):(discrete.dim+continous.dim)){
        gamma[l,1,] = t(as.matrix(x2))[l-discrete.dim,]*r_w.prime
        if(mnar == 1) gamma[l,2,] = t(as.matrix(x2))[l-discrete.dim,]*r_w.prime*model.y
        if(model.discrete.dim >0){
          for(c in (1+mnar+1):(1+mnar+param.num[1])){
            gamma[l,c,] = t(as.matrix(x2))[l-discrete.dim,]*r_w.prime*t(as.matrix(model.x1))[c-(1+mnar),]
          }
        }
        if(model.continous.dim > 0){
          for(c in (1+mnar+param.num[1]+1):(1+mnar+sum(param.num))){
            gamma[l,c,] = t(as.matrix(x2))[l-discrete.dim,]*r_w.prime*t(as.matrix(model.x2))[c-(1+mnar+param.num[1]),]
          }
        }
      }
    }
    return(apply(gamma, c(1, 2), mean))
  }

  obj = function(theta){
    g.m = g(theta)
    G.hat = G(g.m)
    value = t(G.hat)%*%G.hat
    return(ifelse(is.infinite(value) || is.na(value), 10^8, value))
  }

  # theta.m = matrix(rep(0, (1+1+p)), (1+1+p))
  theta.m = matrix(init, M)
  # theta.m = matrix(cc.miss(y, u, r, w, (1+1+p)), (1+1+p))
  conv.err = 10^8
  t = 1

  while (conv.err > 10^(-8) & t < 1000){
    opt = optim(theta.m[,t], obj, method = "L-BFGS-B")
    theta.m = cbind(theta.m, opt$par)
    g.m = g(theta.m[,t+1]); W.hat = W(g.m);
    obj = function(theta){
      g.m = g(theta); G.hat = G(g.m);
      value = t(G.hat)%*%W.hat%*%G.hat
      return(ifelse(is.infinite(value) || is.na(value), 10^8, value))
    }
    conv.err = max(abs(theta.m[,t+1]-theta.m[,t]))
    t = t + 1
  }

  theta.hat = theta.m[, t]
  Gamma.hat = Gamma(theta.hat)
  g.m = g(theta.hat); W.hat = W(g.m);
  S = var(t(g.m))
  # cov.hat = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)/N
  cov.hat = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)%*%t(Gamma.hat)%*%W.hat%*%S%*%W.hat%*%Gamma.hat%*%solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)/N
  se = sqrt(diag(cov.hat))

  return(list(sol.path = theta.m, theta.hat = theta.hat, cov.hat = cov.hat,
              se = se, lower = theta.hat-qnorm(0.975)*se, upper = theta.hat+qnorm(0.975)*se,
              g.m = g.m, GammaW = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)%*%t(Gamma.hat)%*%W.hat,
              model = propensity, w.prime = w.prime, model.x.names = c(model.x1.names, model.x2.names), h = t(cbind(d, x2))))
}

Wang2014.2 = function(auxilliary, model.y, model.x1.names, model.x2.names, w, w.prime){
  propensity = function(y, x, theta, L) 1/w(theta, y, x, L)

  r = dat$r
  x1 = auxilliary[[1]]; x2 = auxilliary[[2]]; auxilliary.y = auxilliary[[3]]; auxilliary.m0 = auxilliary[[4]];
  auxilliary.y[r == 0] = auxilliary.m0[r == 0]
  if(!is.null(x1)) x1 = as.matrix(x1)
  if(!is.null(x2)) x2 = as.matrix(x2)
  model.x1 = NULL; model.x2 = NULL;
  if(!is.null(model.x1.names)) model.x1 = as.matrix(dat[model.x1.names])
  if(!is.null(model.x2.names)) model.x2 = as.matrix(dat[model.x2.names])

  N = length(r)
  discrete.dim = 0
  if(!is.null(x1)) for(j in 1:ncol(x1)) discrete.dim = discrete.dim + length(unique(x1[,j]))
  discrete.dim = ifelse(discrete.dim == 0, 0, discrete.dim - ncol(x1) + 1) # to avoid collinearity
  continous.dim = ifelse(!is.null(x2), ncol(x2), 0)
  auxilliary.dim = discrete.dim + continous.dim + 1

  mnar = ifelse(is.null(model.y), 0, 1)
  model.discrete.dim = ifelse(is.null(model.x1), 0, ncol(model.x1))
  model.continous.dim = ifelse(is.null(model.x2), 0, ncol(model.x2))
  param.num = c(model.discrete.dim, model.continous.dim)

  M = 1+mnar+sum(param.num)
  init = rep(0.5, M)

  d = NULL
  if(!is.null(x1)){
    x1 =  as.data.frame(x1)
    for(j in 1:ncol(x1)) x1[, j] = as.factor(x1[, j])
    d = model.matrix(lm(rep(1, N)~., data =  x1))
  }


  g = function(theta){
    g.m = matrix(NA, auxilliary.dim, N)
    rw = r*w(theta, model.y, cbind(model.x1, model.x2), N)
    if(discrete.dim > 0){
      for(l in 1:discrete.dim){g.m[l,] = t(d)[l,]*(rw-1)}
    }
    if(continous.dim > 0){
      for(l in (discrete.dim+1):(discrete.dim+continous.dim)) g.m[l,] = t(as.matrix(x2))[l-discrete.dim,]*(rw-1)
    }
    g.m[auxilliary.dim,] = auxilliary.y*(rw-1)
    return(g.m)
  }

  G = function(g.m){
    return(matrix(apply(g.m, 1, mean), auxilliary.dim, 1))
  }

  W = function(g.m){
    return(solve(g.m%*%t(g.m)/N))
  }

  # Gamma = function(theta){
  #   gamma = array(NA, dim = c(auxilliary.dim, M, N))
  #   r_w.prime = r*w.prime(theta[1:(M-1)], model.y, cbind(model.x1, model.x2), N)
  #   if(discrete.dim > 0){
  #     for(l in 1:discrete.dim){
  #       gamma[l,1,] = t(d)[l,]*r_w.prime
  #       if(mnar == 1) gamma[l,2,] = t(d)[l,]*r_w.prime*model.y
  #       if(model.discrete.dim > 0){
  #         for(c in (1+mnar+1):(1+mnar+param.num[1])){
  #           gamma[l,c,] = t(d)[l,]*r_w.prime*t(as.matrix(model.x1))[c-(1+mnar),]
  #         }
  #       }
  #       if(model.continous.dim > 0){
  #         for(c in (1+mnar+param.num[1]+1):(1+mnar+sum(param.num))){
  #           gamma[l,c,] = t(d)[l,]*r_w.prime*t(as.matrix(model.x2))[c-(1+mnar+param.num[1]),]
  #         }
  #       }
  #     }
  #   }
  #   if(continous.dim > 0){
  #     for(l in (discrete.dim+1):(discrete.dim+continous.dim)){
  #       gamma[l,1,] = t(as.matrix(x2))[l-discrete.dim,]*r_w.prime
  #       if(mnar == 1) gamma[l,2,] = t(as.matrix(x2))[l-discrete.dim,]*r_w.prime*model.y
  #       if(model.discrete.dim >0){
  #         for(c in (1+mnar+1):(1+mnar+param.num[1])){
  #           gamma[l,c,] = t(as.matrix(x2))[l-discrete.dim,]*r_w.prime*t(as.matrix(model.x1))[c-(1+mnar),]
  #         }
  #       }
  #       if(model.continous.dim > 0){
  #         for(c in (1+mnar+param.num[1]+1):(1+mnar+sum(param.num))){
  #           gamma[l,c,] = t(as.matrix(x2))[l-discrete.dim,]*r_w.prime*t(as.matrix(model.x2))[c-(1+mnar+param.num[1]),]
  #         }
  #       }
  #     }
  #   }
  #   gamma[auxilliary.dim, 1:(1+mnar), ] = t(cbind(1, model.y)*as.vector(r_w.prime*auxilliary.y))
  #   if(length(model.x1.names)+length(model.x2.names) > 0){
  #     gamma[auxilliary.dim, (1+mnar+1):(M-1), ] = t(cbind(model.x1, model.x2)*as.vector(r_w.prime*auxilliary.y))
  #   }
  #   gamma[-auxilliary.dim, M, ] = 0
  #   gamma[auxilliary.dim, M, ] = -1
  #   return(apply(gamma, c(1, 2), mean))
  # }

  # Gamma(theta)
  #
  Gamma = function(theta){
    Gamma.hat = array(NA, dim = c(auxilliary.dim, N, M))
    for(l in 1:auxilliary.dim){
      Gamma.hat[l,,] = jacobian(function(theta.v) g(theta.v)[l,], theta)
    }
    return(apply(Gamma.hat, c(1, 3), mean))
  }

  obj = function(theta){
    g.m = g(theta)
    G.hat = G(g.m)
    value = t(G.hat)%*%G.hat
    return(ifelse(is.infinite(value) || is.na(value), 10^8, value))
  }

  theta.m = matrix(init, M)
  conv.err = 10^8
  t = 1

  while (conv.err > 10^(-8) & t < 1000){
    opt = optim(theta.m[,t], obj, method = "L-BFGS-B")
    theta.m = cbind(theta.m, opt$par)
    g.m = g(theta.m[,t+1]); W.hat = W(g.m);
    obj = function(theta){
      g.m = g(theta); G.hat = G(g.m);
      value = t(G.hat)%*%W.hat%*%G.hat
      return(ifelse(is.infinite(value) || is.na(value), 10^8, value))
    }
    conv.err = max(abs(theta.m[,t+1]-theta.m[,t]))
    t = t + 1
  }

  theta.hat = theta.m[, t]
  Gamma.hat = Gamma(theta.hat)
  g.m = g(theta.hat); W.hat = W(g.m);
  cov.hat = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)/N
  se = sqrt(diag(cov.hat))

  return(list(sol.path = theta.m, theta.hat = theta.hat, cov.hat = cov.hat,
              se = se, lower = theta.hat-qnorm(0.975)*se, upper = theta.hat+qnorm(0.975)*se,
              g.m = g.m, GammaW = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)%*%t(Gamma.hat)%*%W.hat,
              model = propensity, w.prime = w.prime, model.x.names = c(model.x1.names, model.x2.names), h = t(cbind(d, x2))))
}

# Wang2014.2 = function(h.opt, model.x = list(y, u1, u2), r, w, w.prime, init){
#   propensity = function(y, u1, u2, theta, L) 1/w(theta, y, u1, u2, L)
#
#   y = model.x[[1]]; model.u1 = model.x[[2]]; model.u2 = model.x[[3]];
#
#   n = length(r)
#   M = length(init)
#
#   mnar = ifelse(is.null(y), 0, 1)
#   model.p1 = ifelse(is.null(model.u1), 0, ncol(as.matrix(model.u1))); model.p2 = ifelse(is.null(model.u2), 0, ncol(as.matrix(model.u2)));
#   param.num = c(model.p1, model.p2)
#
#   g = function(theta){
#     propensity.fit = as.vector(propensity(y, model.u1, model.u2, theta[1:M], n))
#     g.m = t((r/propensity.fit-1)*h.opt)
#     return(g.m)
#   }
#
#   G = function(g.m){
#     return(matrix(apply(g.m, 1, mean), M, 1))
#   }
#
#   W = function(g.m){
#     return(solve(g.m%*%t(g.m)/n))
#   }
#
#   Gamma = function(theta){
#     gamma = array(NA, dim = c(M, M, n))
#     for(i in 1:M){
#       gamma[i,,] = t(jacobian(function(theta.v) g(theta.v)[i,], theta))
#     }
#     return(apply(gamma, c(1, 2), mean))
#   }
#
#   obj = function(theta){
#     g.m = g(theta)
#     G.hat = G(g.m)
#     value = t(G.hat)%*%G.hat
#     return(ifelse(is.infinite(value) || is.na(value), 10^8, value))
#   }
#
#   # theta.m = matrix(rep(0, (1+1+p)), (1+1+p))
#   theta.m = matrix(init, (1+mnar+sum(param.num)))
#   # theta.m = matrix(cc.miss(y, u, r, w, (1+1+p)), (1+1+p))
#   conv.err = 10^8
#   t = 1
#
#   while (conv.err > 10^(-8) & t < 1000){
#     opt = optim(theta.m[,t], obj, method = "Nelder-Mead")
#     theta.m = cbind(theta.m, opt$par)
#     g.m = g(theta.m[,t+1]); W.hat = W(g.m);
#     obj = function(theta){
#       g.m = g(theta); G.hat = G(g.m);
#       value = t(G.hat)%*%W.hat%*%G.hat
#       return(ifelse(is.infinite(value) || is.na(value), 10^8, value))
#     }
#     conv.err = max(abs(theta.m[,t+1]-theta.m[,t]))
#     t = t + 1
#   }
#
#   theta.hat = theta.m[, t]
#   Gamma.hat = Gamma(theta.hat)
#   g.m = g(theta.hat); W.hat = W(g.m);
#   cov.hat = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)/n
#   se = sqrt(diag(cov.hat))
#
#   return(list(sol.path = theta.m, theta.hat = theta.hat, cov.hat = cov.hat,
#               se = se, lower = theta.hat-qnorm(0.975)*se, upper = theta.hat+qnorm(0.975)*se,
#               g.m = g.m, GammaW = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)%*%t(Gamma.hat)%*%W.hat,
#               model = propensity, w.prime = w.prime, model.x = model.x))
# }


# Wang2014.2 = function(model.x = list(y, u1, u2), r, w, w.prime, init){
#   propensity = function(y, u1, u2, theta, L) 1/w(theta, y, u1, u2, L)
#
#   y = model.x[[1]]; model.u1 = model.x[[2]]; model.u2 = model.x[[3]];
#
#   n = length(r)
#   M = length(init)
#
#   mnar = ifelse(is.null(y), 0, 1)
#   model.p1 = ifelse(is.null(model.u1), 0, ncol(as.matrix(model.u1))); model.p2 = ifelse(is.null(model.u2), 0, ncol(as.matrix(model.u2)));
#   param.num = c(model.p1, model.p2)
#
#   g = function(theta){
#     propensity.fit = as.vector(propensity(y, model.u1, model.u2, theta[1:M], n))
#     propensity.prime = jacobian(function(theta.v) propensity(y, model.u1, model.u2, theta.v, n), theta[1:M])
#
#     g.m = t((r/propensity.fit-1)*(propensity.prime/(1-propensity.fit)))
#     return(g.m)
#   }
#
#   G = function(g.m){
#     return(matrix(apply(g.m, 1, mean), M, 1))
#   }
#
#   W = function(g.m){
#     return(solve(g.m%*%t(g.m)/n))
#   }
#
#   Gamma = function(theta){
#     gamma = array(NA, dim = c(M, M, n))
#     for(i in 1:M){
#       gamma[i,,] = t(jacobian(function(theta.v) g(theta.v)[i,], theta))
#     }
#     return(apply(gamma, c(1, 2), mean))
#   }
#
#   obj = function(theta){
#     g.m = g(theta)
#     G.hat = G(g.m)
#     value = t(G.hat)%*%G.hat
#     return(ifelse(is.infinite(value) || is.na(value), 10^8, value))
#   }
#
#   # theta.m = matrix(rep(0, (1+1+p)), (1+1+p))
#   theta.m = matrix(init, (1+mnar+sum(param.num)))
#   # theta.m = matrix(cc.miss(y, u, r, w, (1+1+p)), (1+1+p))
#   conv.err = 10^8
#   t = 1
#
#   while (conv.err > 10^(-8) & t < 1000){
#     opt = optim(theta.m[,t], obj, method = "Nelder-Mead")
#     theta.m = cbind(theta.m, opt$par)
#     g.m = g(theta.m[,t+1]); W.hat = W(g.m);
#     obj = function(theta){
#       g.m = g(theta); G.hat = G(g.m);
#       value = t(G.hat)%*%W.hat%*%G.hat
#       return(ifelse(is.infinite(value) || is.na(value), 10^8, value))
#     }
#     conv.err = max(abs(theta.m[,t+1]-theta.m[,t]))
#     t = t + 1
#   }
#
#   theta.hat = theta.m[, t]
#   Gamma.hat = Gamma(theta.hat)
#   g.m = g(theta.hat); W.hat = W(g.m);
#   cov.hat = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)/n
#   se = sqrt(diag(cov.hat))
#
#   return(list(sol.path = theta.m, theta.hat = theta.hat, cov.hat = cov.hat,
#               se = se, lower = theta.hat-qnorm(0.975)*se, upper = theta.hat+qnorm(0.975)*se,
#               g.m = g.m, GammaW = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)%*%t(Gamma.hat)%*%W.hat, model = propensity))
# }

cc.miss = function(y, u, r, w, d){
  l = function(theta) -sum(log((1/w(theta, y, u))^r))
  theta.cc = nlminb(rep(0, d), l)$par
  return(theta.cc)
}
