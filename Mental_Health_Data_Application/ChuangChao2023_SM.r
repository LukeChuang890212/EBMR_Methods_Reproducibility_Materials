inv.link.v = c("identity" = function(x) x,
             "binomial" = function(x) exp(x)/(1+exp(x)))

rm.extreme = function(v){
  iqr <- IQR(v)

  # Identify the values outside the lower and upper bounds
  lower_bound <- quantile(v, 0.25) - 3 * iqr
  upper_bound <- quantile(v, 0.75) + 3 * iqr

  # Exclude the extreme values
  return(v[v >= lower_bound & v <= upper_bound])
}

impute.m1 = function(outcome.list, beta.v, dim.beta, u1, u2, z1, z2, L){
  m1.m = matrix(NA, L, K)
  X = list()
  for(k in 1:K){
    X[[k]] = model.matrix(lm(paste(outcome.list[[k]]$call)[2], data = data.frame(y = rep(1, L), u1 = u1, u2 = u2, z1 = z1, z2 = z2)))
    m1.m[, k] = inv.link(X[[k]]%*%beta.v[(sum(dim.beta[0:(k-1)])+1):sum(dim.beta[1:k])])
  }
  return(list(m1.m = m1.m, X = X))
}

# m1k = function(beta.v, outcome.fit, u1.value, u2.value, z1.value, z2.value, L){
#   X = model.matrix(lm(paste(outcome.fit$call)[2], data = data.frame(y = rep(1, L), u1 = u1.value, u2 = u2.value, z1 = z1.value, z2 = z2.value)))
#   return(inv.link(X%*%beta.v))
# }
#
# m1 = function(beta.v, w.m1.v, u1.value, u2.value, z1.value, z2.value, L){
#   m1.m = matrix(NA, L, K)
#   for(k in 1:K){
#     m1.m[k] = m1k(beta.v[(sum(dim.beta[0:(k-1)])+1):sum(dim.beta[1:k])], outcome.list[[k]], u1.value, u2.value, z1.value, z2.value, L)
#   }
#   m1 = m1.m%*%(w.m1.v^2)/sum(w.m1.v^2)
#   return(m1)
# }

impute.m1 = function(outcome.list, compressed.m1, pi.fit.list, pi.list, alpha.list, w.m1.fit, w.pi.fit, beta.pi.prime, situation, family = "gaussian", respondent = FALSE){
  norm.adj = function(w){
    if(length(w) > 1){
      return((diag(2*w)*sum(w^2)-2*(w^2)%*%t(w))/(sum(w^2)^2))
    }else{
      return(0)
    }
  }

  w.m1 = (w.m1.fit$coef^2)/sum(w.m1.fit$coef^2)
  w.pi = (w.pi.fit$theta.hat^2)/sum(w.pi.fit$theta.hat^2)

  J = length(pi.list)
  K = length(outcome.list)
  dim.alpha = unlist(lapply(alpha.list, length))
  # beta.list = lapply(outcome.list, coefficients)
  # dim.beta = unlist(lapply(beta.list, length))

  L = length(r)
  sigma = NA
  if(family == "gaussian") sigma = sd(y[r == 1]-compressed.m1[r == 1])

  m1 = rep(NA, sum(r == ifelse(respondent, 1, 0)))
  id = which(r == ifelse(respondent, 1, 0))
  G.x = matrix(0, length(id), sum(dim.alpha))
  B.x = matrix(0, length(id), J)
  S.x = matrix(0, length(id), sum(dim.alpha))
  # h.opt.x = matrix(0, length(id), sum(dim.alpha))

  start = Sys.time()
  if(family == "nonparametric"){
  }else{
    # For parametric method
    rgen = switch(family,
                  gaussian  = function(m1.x) rnorm(L, mean = m1.x, sd = sigma),
                  binomial = function(m1.x) rbinom(L, size = 1, prob = m1.x))
    inv.link = switch(family,
                      gaussian = function(x) x, binomial = function(x) exp(x)/(1+exp(x)))

    for(i in 1:length(id)){
      dat.value = dat[id[i],]
      # z1.value = as.numeric(z1[id[i]])
      # z2.value = ifelse(is.null(z2[id[i]]), NA, z2[id[i]])
      # u1.value = ifelse(is.null(u1[id[i]]), NA, u1[id[i]])
      # u2.value = ifelse(is.null(u2[id[i]]), NA, u2[id[i]])

      m1.x.v = rep(NA, K)
      for(k in 1:K){
        # m1.x.v[k] = predict(outcome.list[[k]],
        #                     data.frame(z1 = z1.value, z2 = z2.value, u1 = u1.value, u2 = u2.value), type = "response")
        m1.x.v[k] = predict(outcome.list[[k]], newdata = dat.value, type = "response")
      }
      m1.x = as.numeric(m1.x.v%*%w.m1)
      # m1.x = as.numeric(impute.m1(outcome.list, unlist(beta.list), dim.beta, z1 = z1.value, z2 = z2.value, u1 = u1.value, u2 = u2.value, 1)$m1.m%*%w.m1)

      y.s = rgen(m1.x)

      pi.x.m = matrix(NA, L, J)
      pi.x1.alpha.m = matrix(NA, L, sum(dim.alpha))
      for(j in 1:J){
        # pi.x.m[, j] = pi.list[[j]](y.s, u1.value, u2.value, alpha.list[[j]], L)
        # pi.x0.alpha.m[, (sum(dim.alpha[0:(j-1)])+1):sum(dim.alpha[1:j])] = jacobian(function(alpha.v) pi.list[[j]](y.s, u1.value, u2.value, alpha.v, L), alpha.list[[j]])
        pi.x.m[, j] = pi.list[[j]](y.s, rep(1, L)%*%as.matrix(dat.value[pi.fit.list[[j]]$model.x.names]), alpha.list[[j]], L)
        pi.x1.alpha.m[, (sum(dim.alpha[0:(j-1)])+1):sum(dim.alpha[1:j])] = jacobian(function(alpha.v) pi.list[[j]](y.s, rep(1, L)%*%as.matrix(dat.value[pi.fit.list[[j]]$model.x.names]), alpha.v, L), alpha.list[[j]])
      }
      pi.x = as.vector(pi.x.m%*%w.pi)

      m1[i] = mean(y.s)

      # S.x[i,] = (pi.x1.alpha.m/(1-pi.x))[1,]

      # h.opt.x[i,] = apply((pi.x0.alpha.m/((pi.x)^2))/mean((1-pi.x)/((pi.x)^2)), 2, mean)
    }
  }

  Sys.time()-start

  return(list(m1 = m1))
}

impute.m0 = function(outcome.list, compressed.m1, pi.fit.list, pi.list, alpha.list, w.m1.fit, w.pi.fit, beta.pi.prime, situation, family = "gaussian", respondent = FALSE){
  norm.adj = function(w){
    if(length(w) > 1){
      return((diag(2*w)*sum(w^2)-2*(w^2)%*%t(w))/(sum(w^2)^2))
    }else{
      return(0)
    }
  }

  w.m1 = (w.m1.fit$coef^2)/sum(w.m1.fit$coef^2)
  w.pi = (w.pi.fit$theta.hat^2)/sum(w.pi.fit$theta.hat^2)

  J = length(pi.list)
  K = length(outcome.list)
  dim.alpha = unlist(lapply(alpha.list, length))
  # beta.list = lapply(outcome.list, coefficients)
  # dim.beta = unlist(lapply(beta.list, length))

  L = length(r)
  sigma = NA
  if(family == "gaussian") sigma = sd(y[r == 1]-compressed.m1[r == 1])

  a.v = m0 = rep(NA, sum(r == ifelse(respondent, 1, 0)))
  id = which(r == ifelse(respondent, 1, 0))
  G.x = matrix(0, length(id), sum(dim.alpha))
  B.x = matrix(0, length(id), J)
  # S.x = matrix(0, length(id), sum(dim.alpha))
  # h.opt.x = matrix(0, length(id), sum(dim.alpha))

  start = Sys.time()
  if(family == "nonparametric"){
    # For nonparametric method
    impute.y = function(beta.v, w.m1.v){
      m1.m = m1.x.m = matrix(NA, L, K)
      m1.m = impute.m1(outcome.list, beta.v, dim.beta, u1[r == 1], u2[r == 1], z1[r == 1], z2[r == 1], sum(r))$m1.m
      m1.x.m = impute.m1(outcome.list, beta.v, dim.beta, u1.value, u2.value, z1.value, z2.value, L)$m1.m
      return((m1.x.m-m1.m)%*%w.m1.v+y[r == 1])
    }

    for(i in 1:length(id)){
      z1.value = as.numeric(z1[id[i]])
      z2.value = ifelse(is.null(z2[id[i]]), NA, z2[id[i]])
      u1.value = ifelse(is.null(u1[id[i]]), NA, u1[id[i]])
      u2.value = ifelse(is.null(u2[id[i]]), NA, u2[id[i]])

      y.s = as.vector(impute.y(unlist(beta.list), w.m1))

      pi.x.m = matrix(NA, L, J)
      pi.x0.alpha.m = matrix(NA, L, sum(dim.alpha))
      for(j in 1:J){
        pi.x.m[, j] = pi.list[[j]](y.s, u1.value, u2.value, alpha.list[[j]], L)
        pi.x0.alpha.m[, (sum(dim.alpha[0:(j-1)])+1):sum(dim.alpha[1:j])] = jacobian(function(alpha.v) pi.list[[j]](y.s, u1.value, u2.value, alpha.v, L), alpha.list[[j]])
      }
      pi.x = pi.x.m%*%w.pi

      y.s.beta.m = jacobian(function(beta.v) impute.y(beta.v, w.m1), unlist(beta.list))
      y.s.w.m1.m = jacobian(function(w.m1.v) impute.y(unlist(beta.list), w.m1.v), w.m1)

      pi.x.beta.m = array(NA, dim = c(L, sum(dim.beta), J))
      pi.x.w.m1.m = array(NA, dim = c(L, K, J))
      for(j in 1:J){
        pi.x.beta.m[,, j] = jacobian(function(beta.v) pi.list[[j]](impute.y(beta.v, w.m1), u1.value, u2.value, alpha.list[[j]], L), unlist(beta.list))
        pi.x.w.m1.m[,, j] = jacobian(function(w.m1.v) pi.list[[j]](impute.y(unlist(beta.list), w.m1.v), u1.value, u2.value, alpha.list[[j]], L), w.m1)
      }

      a = a.v[i] = mean((1-pi.x)/pi.x)
      m0[i] = mean(((1-pi.x)/pi.x)/a*(y.s))
      var.hat[i] = var(((1-pi.x)/pi.x)/a*(y.s))
      gamma.xa[i] = mean(-((1-pi.x)/pi.x)/(a^2)*(y.s))
      odds.var.hat[i] = var((1-pi.x)/pi.x)

      G.x[i, ] = cbind(as.vector((y.s/(a*(pi.x^2)))*(m0[i]/y.s-1))*t(t(pi.x0.alpha.m)*rep(w.pi, dim.alpha)),
                       as.vector((1/pi.x-1)/a)*y.s.beta.m-as.vector((y.s/((pi.x^2)*a))*(m0[i]/y.s-1))*apply(pi.x.beta.m, c(1, 2), weighted.mean, w.pi),
                       as.vector(y.s/(a*(pi.x^2))*(m0[i]/y.s-1))*pi.x.m,
                       as.vector((1/pi.x-1)/a)*y.s.w.m1.m-as.vector(y.s/(a*(pi.x^2))*(m0[i]/y.s-1))*apply(pi.x.w.m1.m, c(1, 2), weighted.mean, w.pi)) %>%
        apply(., 2, function(G.xj) mean(G.xj))
    }
  }else{
    # For parametric method
    rgen = switch(family,
                  gaussian  = function(m1.x) rnorm(L, mean = m1.x, sd = sigma),
                  binomial = function(m1.x) rbinom(L, size = 1, prob = m1.x))
    inv.link = switch(family,
                      gaussian = function(x) x, binomial = function(x) exp(x)/(1+exp(x)))

    for(i in 1:length(id)){
      dat.value = dat[id[i],]
      # z1.value = as.numeric(z1[id[i]])
      # z2.value = ifelse(is.null(z2[id[i]]), NA, z2[id[i]])
      # u1.value = ifelse(is.null(u1[id[i]]), NA, u1[id[i]])
      # u2.value = ifelse(is.null(u2[id[i]]), NA, u2[id[i]])

      m1.x.v = rep(NA, K)
      for(k in 1:K){
        # m1.x.v[k] = predict(outcome.list[[k]],
        #                     data.frame(z1 = z1.value, z2 = z2.value, u1 = u1.value, u2 = u2.value), type = "response")
        m1.x.v[k] = predict(outcome.list[[k]], newdata = dat.value, type = "response")
      }
      m1.x = as.numeric(m1.x.v%*%w.m1)
      # m1.x = as.numeric(impute.m1(outcome.list, unlist(beta.list), dim.beta, z1 = z1.value, z2 = z2.value, u1 = u1.value, u2 = u2.value, 1)$m1.m%*%w.m1)

      y.s = rgen(m1.x)

      pi.x.m = matrix(NA, L, J)
      pi.x0.alpha.m = matrix(NA, L, sum(dim.alpha))
      for(j in 1:J){
        # pi.x.m[, j] = pi.list[[j]](y.s, u1.value, u2.value, alpha.list[[j]], L)
        # pi.x0.alpha.m[, (sum(dim.alpha[0:(j-1)])+1):sum(dim.alpha[1:j])] = jacobian(function(alpha.v) pi.list[[j]](y.s, u1.value, u2.value, alpha.v, L), alpha.list[[j]])
        pi.x.m[, j] = pi.list[[j]](y.s, rep(1, L)%*%as.matrix(dat.value[pi.fit.list[[j]]$model.x.names]), alpha.list[[j]], L)
        pi.x0.alpha.m[, (sum(dim.alpha[0:(j-1)])+1):sum(dim.alpha[1:j])] = jacobian(function(alpha.v) pi.list[[j]](y.s, rep(1, L)%*%as.matrix(dat.value[pi.fit.list[[j]]$model.x.names]), alpha.v, L), alpha.list[[j]])
      }
      pi.x = as.vector(pi.x.m%*%w.pi)

      a = a.v[i] = mean((1-pi.x)/pi.x)
      m0[i] = mean(((1-pi.x)/pi.x)/a*(y.s))

      m0.alpha = apply(-as.vector((y.s/a)*((pi.x)^(-2)))*(t(t(pi.x0.alpha.m)*rep(w.pi, dim.alpha))), 2, mean)
      m0.a = apply(-as.matrix((1/pi.x-1)/(a^2)*y.s), 2, mean)
      m0.beta.pi = apply(-as.vector((y.s/a)*((pi.x)^(-2)))*(pi.x.m%*%norm.adj(w.pi.fit$theta.hat)), 2, mean)
      a.alpha = apply(-(pi.x)^(-2)*(t(t(pi.x0.alpha.m)*rep(w.pi, dim.alpha))), 2, mean)
      a.beta.pi = apply(-(pi.x)^(-2)*(pi.x.m%*%(norm.adj(w.pi.fit$theta.hat))), 2, mean)

      G.x[i,] = m0.alpha+m0.a%*%a.alpha+(m0.beta.pi+m0.a%*%a.beta.pi)%*%beta.pi.prime
      B.x[i,] = m0.beta.pi+m0.a%*%a.beta.pi
      # S.x[i,] = (pi.x0.alpha.m/(1-pi.x))[1,]

      # h.opt.x[i,] = apply((pi.x0.alpha.m/((pi.x)^2))/mean((1-pi.x)/((pi.x)^2)), 2, mean)
    }
  }

  Sys.time()-start

  return(list(m0 = m0, a.v = a.v, G.x = G.x, B.x = B.x))
}

impute.m = function(outcome.list, compressed.m1, pi.fit.list, pi.list, alpha.list, w.m1.fit, w.pi.fit, beta.pi.prime, situation, family = "gaussian", respondent = FALSE){
  norm.adj = function(w){
    if(length(w) > 1){
      return((diag(2*w)*sum(w^2)-2*(w^2)%*%t(w))/(sum(w^2)^2))
    }else{
      return(0)
    }
  }

  w.m1 = (w.m1.fit$coef^2)/sum(w.m1.fit$coef^2)
  w.pi = (w.pi.fit$theta.hat^2)/sum(w.pi.fit$theta.hat^2)

  J = length(pi.list)
  K = length(outcome.list)
  dim.alpha = unlist(lapply(alpha.list, length))

  L = length(r)
  sigma = NA
  if(family == "gaussian") sigma = sd(y[r == 1]-compressed.m1[r == 1])

  a.v = m = rep(NA, sum(r == ifelse(respondent, 1, 0)))
  id = which(r == ifelse(respondent, 1, 0))
  G.x = matrix(0, length(id), sum(dim.alpha))
  B.x = matrix(0, length(id), J)
  # h.opt.x = matrix(0, length(id), sum(dim.alpha))

  start = Sys.time()
  if(family == "nonparametric"){
  }else{
    # For parametric method
    rgen = switch(family,
                  gaussian  = function(m1.x) rnorm(L, mean = m1.x, sd = sigma),
                  binomial = function(m1.x) rbinom(L, size = 1, prob = m1.x))
    inv.link = switch(family,
                      gaussian = function(x) x, binomial = function(x) exp(x)/(1+exp(x)))

    for(i in 1:length(id)){
      dat.value = dat[id[i],]

      m1.x.v = rep(NA, K)
      for(k in 1:K){
        m1.x.v[k] = predict(outcome.list[[k]], newdata = dat.value, type = "response")
      }
      m1.x = as.numeric(m1.x.v%*%w.m1)
      # m1.x = as.numeric(impute.m1(outcome.list, unlist(beta.list), dim.beta, z1 = z1.value, z2 = z2.value, u1 = u1.value, u2 = u2.value, 1)$m1.m%*%w.m1)

      y.s = rgen(m1.x)

      pi.x.m = matrix(NA, L, J)
      pi.x0.alpha.m = matrix(NA, L, sum(dim.alpha))
      for(j in 1:J){
        pi.x.m[, j] = pi.list[[j]](y.s, rep(1, L)%*%as.matrix(dat.value[pi.fit.list[[j]]$model.x.names]), alpha.list[[j]], L)
        pi.x0.alpha.m[, (sum(dim.alpha[0:(j-1)])+1):sum(dim.alpha[1:j])] = jacobian(function(alpha.v) pi.list[[j]](y.s, rep(1, L)%*%as.matrix(dat.value[pi.fit.list[[j]]$model.x.names]), alpha.v, L), alpha.list[[j]])
      }
      pi.x = as.vector(pi.x.m%*%w.pi)

      a = a.v[i] = mean(1/pi.x)
      m[i] = mean((1/pi.x)/a*(y.s))

      m.alpha = apply(-as.vector((y.s/a)*((pi.x)^(-2)))*(t(t(pi.x0.alpha.m)*rep(w.pi, dim.alpha))), 2, mean)
      m.a = apply(-as.matrix((1/pi.x)/(a^2)*y.s), 2, mean)
      m.beta.pi = apply(-as.vector((y.s/a)*((pi.x)^(-2)))*(pi.x.m%*%norm.adj(w.pi.fit$theta.hat)), 2, mean)
      a.alpha = apply(-(pi.x)^(-2)*(t(t(pi.x0.alpha.m)*rep(w.pi, dim.alpha))), 2, mean)
      a.beta.pi = apply(-(pi.x)^(-2)*(pi.x.m%*%(norm.adj(w.pi.fit$theta.hat))), 2, mean)

      G.x[i,] = m.alpha+m.a%*%a.alpha+(m.beta.pi+m.a%*%a.beta.pi)%*%beta.pi.prime
      B.x[i,] = m.beta.pi+m.a%*%a.beta.pi
      # S.x[i,] = (pi.x0.alpha.m/(1-pi.x))[1,]
      # h.opt.x[i,] = apply((pi.x0.alpha.m/((pi.x)^2))/mean(1/((pi.x)^2)), 2, mean)
    }
  }
  Sys.time()-start

  return(list(m = m, a.v = a.v, G.x = G.x, B.x = B.x))
}

impute.weighted.m = function(outcome.list, compressed.m1, pi.fit.list, pi.list, alpha.list, w.m1.fit, w.pi.fit, beta.pi.prime, situation, family = "gaussian"){
  norm.adj = function(w){
    if(length(w) > 1){
      return((diag(2*w)*sum(w^2)-2*(w^2)%*%t(w))/(sum(w^2)^2))
    }else{
      return(0)
    }
  }

  w.m1 = (w.m1.fit$coef^2)/sum(w.m1.fit$coef^2)
  w.pi = (w.pi.fit$theta.hat^2)/sum(w.pi.fit$theta.hat^2)

  J = length(pi.list)
  K = length(outcome.list)
  dim.alpha = unlist(lapply(alpha.list, length))
  dim.L = nrow(pi.fit.list[[1]]$h)

  L = N = length(r)
  sigma = NA
  if(family == "gaussian") sigma = sd(y[r == 1]-compressed.m1[r == 1])

  a.v = weighted.m = rep(NA, N)
  G.x = matrix(0, N, sum(dim.alpha))
  B.x = matrix(0, N, J)
  weighted.S.x = matrix(0, N, sum(dim.alpha))

  start = Sys.time()
  if(family == "nonparametric"){
  }else{
    # For parametric method
    rgen = switch(family,
                  gaussian  = function(m1.x) rnorm(L, mean = m1.x, sd = sigma),
                  binomial = function(m1.x) rbinom(L, size = 1, prob = m1.x))
    inv.link = switch(family,
                      gaussian = function(x) x, binomial = function(x) exp(x)/(1+exp(x)))

    for(i in 1:N){
      dat.value = dat[i,]

      m1.x.v = rep(NA, K)
      for(k in 1:K){
        m1.x.v[k] = predict(outcome.list[[k]], newdata = dat.value, type = "response")
      }
      m1.x = as.numeric(m1.x.v%*%w.m1)
      # m1.x = as.numeric(impute.m1(outcome.list, unlist(beta.list), dim.beta, z1 = z1.value, z2 = z2.value, u1 = u1.value, u2 = u2.value, 1)$m1.m%*%w.m1)

      y.s = rgen(m1.x)

      pi.x.m = matrix(NA, L, J)
      pi.x0.alpha.m = matrix(NA, L, sum(dim.alpha))
      for(j in 1:J){
        pi.x.m[, j] = pi.list[[j]](y.s, rep(1, L)%*%as.matrix(dat.value[pi.fit.list[[j]]$model.x.names]), alpha.list[[j]], L)
        pi.x0.alpha.m[, (sum(dim.alpha[0:(j-1)])+1):sum(dim.alpha[1:j])] = jacobian(function(alpha.v) pi.list[[j]](y.s, rep(1, L)%*%as.matrix(dat.value[pi.fit.list[[j]]$model.x.names]), alpha.v, L), alpha.list[[j]])
      }
      pi.x = as.vector(pi.x.m%*%w.pi)

      a = a.v[i] = mean((1-pi.x)/((pi.x)^2))
      weighted.m[i] = mean((1-pi.x)/((pi.x)^2)/a*(y.s))

      weighted.m.alpha = apply(as.vector((y.s/a)*((pi.x)^(-2)-2*(pi.x)^(-3)))*(t(t(pi.x0.alpha.m)*rep(w.pi, dim.alpha))), 2, mean)
      weighted.m.a = apply(-as.matrix((1/pi.x-1)/(a^2)*y.s), 2, mean)
      weighted.m.beta.pi = apply(as.vector((y.s/a)*((pi.x)^(-2)-2*(pi.x)^(-3)))*(pi.x.m%*%norm.adj(w.pi.fit$theta.hat)), 2, mean)
      a.alpha = apply(((pi.x)^(-2)-2*(pi.x)^(-3))*(t(t(pi.x0.alpha.m)*rep(w.pi, dim.alpha))), 2, mean)
      a.beta.pi = apply(((pi.x)^(-2)-2*(pi.x)^(-3))*(pi.x.m%*%(norm.adj(w.pi.fit$theta.hat))), 2, mean)

      G.x[i,] = weighted.m.alpha+weighted.m.a%*%a.alpha+(weighted.m.beta.pi+weighted.m.a%*%a.beta.pi)%*%beta.pi.prime
      B.x[i,] = weighted.m.beta.pi+weighted.m.a%*%a.beta.pi

      weighted.S.x[i,] = apply((1-pi.x)/((pi.x)^2)/a*(pi.x0.alpha.m/(1-pi.x)), 2, mean)
    }
  }
  Sys.time()-start

  return(list(weighted.m = weighted.m, a.v = a.v, G.x = G.x, B.x = B.x, weighted.S.x = weighted.S.x))
}

impute.weighted.m.v2 = function(outcome.list, compressed.m1, pi.fit.list, pi.list, alpha.list, w.m1.fit, w.pi.fit, beta.pi.prime, situation, family = "gaussian"){
  rf = randomForest(x = t(pi.fit.list[[1]]$h), y = (r/compressed.pi)/(r/compressed.pi-1)*y)
  return(list(weighted.m = predict(rf)))
}

calculate.h.opt = function(outcome.list, compressed.m1, pi.fit.list, pi.list, alpha.list, w.m1.fit, w.pi.fit, beta.pi.prime, situation, family = "gaussian"){
  norm.adj = function(w){
    if(length(w) > 1){
      return((diag(2*w)*sum(w^2)-2*(w^2)%*%t(w))/(sum(w^2)^2))
    }else{
      return(0)
    }
  }

  w.m1 = (w.m1.fit$coef^2)/sum(w.m1.fit$coef^2)
  w.pi = (w.pi.fit$theta.hat^2)/sum(w.pi.fit$theta.hat^2)

  J = length(pi.list)
  K = length(outcome.list)
  dim.alpha = unlist(lapply(alpha.list, length))
  dim.L = nrow(pi.fit.list[[1]]$h)

  L = N
  sigma = NA
  if(family == "gaussian") sigma = sd(y[r == 1]-compressed.m1[r == 1])

  E.yg.x = matrix(NA, dim.L, N)
  E.gg.x = array(NA, dim = c(dim.L, dim.L, N))

  start = Sys.time()
  if(family == "nonparametric"){
    # For nonparametric method
    impute.y = function(beta.v, w.m1.v){
      m1.m = m1.x.m = matrix(NA, L, K)
      m1.m = impute.m1(outcome.list, beta.v, dim.beta, u1[r == 1], u2[r == 1], z1[r == 1], z2[r == 1], sum(r))$m1.m
      m1.x.m = impute.m1(outcome.list, beta.v, dim.beta, u1.value, u2.value, z1.value, z2.value, L)$m1.m
      return((m1.x.m-m1.m)%*%w.m1.v+y[r == 1])
    }

    for(i in 1:length(id)){
      z1.value = as.numeric(z1[id[i]])
      z2.value = ifelse(is.null(z2[id[i]]), NA, z2[id[i]])
      u1.value = ifelse(is.null(u1[id[i]]), NA, u1[id[i]])
      u2.value = ifelse(is.null(u2[id[i]]), NA, u2[id[i]])

      y.s = as.vector(impute.y(unlist(beta.list), w.m1))

      pi.x.m = matrix(NA, L, J)
      pi.x0.alpha.m = matrix(NA, L, sum(dim.alpha))
      for(j in 1:J){
        pi.x.m[, j] = pi.list[[j]](y.s, u1.value, u2.value, alpha.list[[j]], L)
        pi.x0.alpha.m[, (sum(dim.alpha[0:(j-1)])+1):sum(dim.alpha[1:j])] = jacobian(function(alpha.v) pi.list[[j]](y.s, u1.value, u2.value, alpha.v, L), alpha.list[[j]])
      }
      pi.x = pi.x.m%*%w.pi

      y.s.beta.m = jacobian(function(beta.v) impute.y(beta.v, w.m1), unlist(beta.list))
      y.s.w.m1.m = jacobian(function(w.m1.v) impute.y(unlist(beta.list), w.m1.v), w.m1)

      pi.x.beta.m = array(NA, dim = c(L, sum(dim.beta), J))
      pi.x.w.m1.m = array(NA, dim = c(L, K, J))
      for(j in 1:J){
        pi.x.beta.m[,, j] = jacobian(function(beta.v) pi.list[[j]](impute.y(beta.v, w.m1), u1.value, u2.value, alpha.list[[j]], L), unlist(beta.list))
        pi.x.w.m1.m[,, j] = jacobian(function(w.m1.v) pi.list[[j]](impute.y(unlist(beta.list), w.m1.v), u1.value, u2.value, alpha.list[[j]], L), w.m1)
      }

      a = a.v[i] = mean((1-pi.x)/pi.x)
      m[i] = mean(((1-pi.x)/pi.x)/a*(y.s))
      var.hat[i] = var(((1-pi.x)/pi.x)/a*(y.s))
      gamma.xa[i] = mean(-((1-pi.x)/pi.x)/(a^2)*(y.s))
      odds.var.hat[i] = var((1-pi.x)/pi.x)

      G.x[i, ] = cbind(as.vector((y.s/(a*(pi.x^2)))*(m[i]/y.s-1))*t(t(pi.x0.alpha.m)*rep(w.pi, dim.alpha)),
                       as.vector((1/pi.x-1)/a)*y.s.beta.m-as.vector((y.s/((pi.x^2)*a))*(m[i]/y.s-1))*apply(pi.x.beta.m, c(1, 2), weighted.mean, w.pi),
                       as.vector(y.s/(a*(pi.x^2))*(m[i]/y.s-1))*pi.x.m,
                       as.vector((1/pi.x-1)/a)*y.s.w.m1.m-as.vector(y.s/(a*(pi.x^2))*(m[i]/y.s-1))*apply(pi.x.w.m1.m, c(1, 2), weighted.mean, w.pi)) %>%
        apply(., 2, function(G.xj) mean(G.xj))
    }
  }else{
    # For parametric method
    rgen = switch(family,
                  gaussian  = function(m1.x) rnorm(L, mean = m1.x, sd = sigma),
                  binomial = function(m1.x) rbinom(L, size = 1, prob = m1.x))
    inv.link = switch(family,
                      gaussian = function(x) x, binomial = function(x) exp(x)/(1+exp(x)))

    for(i in 1:N){
      dat.value = dat[i,]

      m1.x.v = rep(NA, K)
      for(k in 1:K){
        m1.x.v[k] = predict(outcome.list[[k]], newdata = dat.value, type = "response")
      }
      m1.x = as.numeric(m1.x.v%*%w.m1)
      # m1.x = as.numeric(impute.m1(outcome.list, unlist(beta.list), dim.beta, z1 = z1.value, z2 = z2.value, u1 = u1.value, u2 = u2.value, 1)$m1.m%*%w.m1)

      y.s = rgen(m1.x)

      pi.x.m = matrix(NA, L, J)
      pi.x0.alpha.m = matrix(NA, L, sum(dim.alpha))
      for(j in 1:J){
        pi.x.m[, j] = pi.list[[j]](y.s, rep(1, L)%*%as.matrix(dat.value[pi.fit.list[[j]]$model.x.names]), alpha.list[[j]], L)
        pi.x0.alpha.m[, (sum(dim.alpha[0:(j-1)])+1):sum(dim.alpha[1:j])] = jacobian(function(alpha.v) pi.list[[j]](y.s, rep(1, L)%*%as.matrix(dat.value[pi.fit.list[[j]]$model.x.names]), alpha.v, L), alpha.list[[j]])
      }
      pi.x = as.vector(pi.x.m%*%w.pi)

      a = mean(1/pi.x)
      h.x = as.matrix(pi.fit.list[[1]]$h[, i])
      E.yg.x[, i] = mean((1/pi.x)/a*(1/pi.x-1)*y.s)*h.x
      E.gg.x[, ,i] = mean((1/pi.x)/a*(1/pi.x-1))*h.x%*%t(h.x)

      # m.alpha = apply(-as.vector((y.s/a)*((pi.x)^(-2)))*(t(t(pi.x0.alpha.m)*rep(w.pi, dim.alpha))), 2, mean)
      # m.a = apply(-as.matrix((1/pi.x)/(a^2)*y.s), 2, mean)
      # m.beta.pi = apply(-as.vector((y.s/a)*((pi.x)^(-2)))*(pi.x.m%*%norm.adj(w.pi.fit$theta.hat)), 2, mean)
      # a.alpha = apply(-(pi.x)^(-2)*(t(t(pi.x0.alpha.m)*rep(w.pi, dim.alpha))), 2, mean)
      # a.beta.pi = apply(-(pi.x)^(-2)*(pi.x.m%*%(norm.adj(w.pi.fit$theta.hat))), 2, mean)
      #
      # G.x[i,] = m.alpha+m.a%*%a.alpha+(m.beta.pi+m.a%*%a.beta.pi)%*%beta.pi.prime
      # B.x[i,] = m.beta.pi+m.a%*%a.beta.pi

      # h.opt.x[i,] = apply((pi.x0.alpha.m/((pi.x)^2))/mean(1/((pi.x)^2)), 2, mean)
    }
  }
  Sys.time()-start

  E.yg = apply(E.yg.x, 1, mean)
  E.gg = apply(E.gg.x, c(1, 2), mean)
  h.opt = as.vector(E.yg%*%solve(E.gg)%*%pi.fit.list[[1]]$h)

  return(h.opt)
}

find.w.pi = function(auxilliary, pi.m, r, init, perturb = rep(1, N), ortho = FALSE){
  eta = diag(rep(1, ncol(pi.m)))
  if(ortho){
     eta = eigen(cov(pi.m))$vector
  }

  x1 = auxilliary[[1]]; x2 = auxilliary[[2]];
  if(!is.null(x1)) x1 = as.matrix(x1)
  if(!is.null(x2)) x2 = as.matrix(x2)

  N = length(r)
  discrete.dim = 0
  if(!is.null(x1)) for(j in 1:ncol(x1)) discrete.dim = discrete.dim + length(unique(x1[,j]))
  discrete.dim = ifelse(discrete.dim == 0, 0, discrete.dim - ncol(x1) + 1) # to avoid collinearity
  continous.dim = ifelse(!is.null(x2), ncol(x2), 0)
  auxilliary.dim = discrete.dim + continous.dim

  M = ncol(pi.m)
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
    rw = rep(0, N)
    rw[r == 1] = 1/(pi.m%*%eta%*%theta)
    if(discrete.dim > 0){
      for(l in 1:discrete.dim){g.m[l,] = t(d)[l,]*(rw-1)}
    }
    if(continous.dim > 0){
      for(l in (discrete.dim+1):(discrete.dim+continous.dim)) g.m[l,] = t(x2)[l-discrete.dim,]*(rw-1)
    }
    g.m = t(perturb*t(g.m))
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
    for(i in 1:auxilliary.dim){
      gamma[i,,] = t(jacobian(function(theta.v) g(theta.v)[i,], theta))
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
    opt = optim(theta.m[, t], obj, method = "L-BFGS-B",
                lower = rep(-Inf, ncol(pi.m)), upper = rep(Inf, ncol(pi.m)))
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

  theta.hat = as.vector(eta%*%theta.hat)
  cov.hat = eta%*%cov.hat%*%t(eta)
  GammaW = eta%*%solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)%*%t(Gamma.hat)%*%W.hat
  se = sqrt(diag(cov.hat))

  return(list(sol.path = theta.m, theta.hat = theta.hat, cov.hat = cov.hat,
              se = se, lower = theta.hat-qnorm(0.975)*se, upper = theta.hat+qnorm(0.975)*se,
              g.m = g.m, GammaW = GammaW, h = t(cbind(d, x2)), eta = eta))
}

ChuangChao2023.GMM = function(auxilliary, pi.m, y, r, init){
  x1 = auxilliary[[1]]; x2 = auxilliary[[2]];
  if(!is.null(x1)) x1 = as.matrix(x1)
  if(!is.null(x2)) x2 = as.matrix(x2)

  N = length(r)
  discrete.dim = 0
  if(!is.null(x1)) for(j in 1:ncol(x1)) discrete.dim = discrete.dim + length(unique(x1[,j]))
  discrete.dim = ifelse(discrete.dim == 0, 0, discrete.dim - ncol(x1) + 1) # to avoid collinearity
  continous.dim = ifelse(!is.null(x2), ncol(x2), 0)
  auxilliary.dim = discrete.dim + continous.dim + 1
  M = length(init)

  d = NULL
  if(!is.null(x1)){
    x1 =  as.data.frame(x1)
    for(j in 1:ncol(x1)) x1[, j] = as.factor(x1[, j])
    d = model.matrix(lm(rep(1, N)~., data =  x1))
  }
  # x1 = as.numeric(as.matrix(x1, ncol(x1)))

  g = function(theta){
    g.m = matrix(NA, auxilliary.dim, N)
    rw = rep(0, N)
    rw[r == 1] = 1/(pi.m%*%theta[-M])
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
    for(i in 1:auxilliary.dim){
      gamma[i,,] = t(jacobian(function(theta.v) g(theta.v)[i,], theta))
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
    opt = optim(theta.m[, t], obj, method = "L-BFGS-B",
                lower = rep(-Inf, ncol(pi.m)), upper = rep(Inf, ncol(pi.m)))
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
  cov.hat = NULL; GammaW = matrix(0, M, auxilliary.dim);
  if(det(t(Gamma.hat)%*%W.hat%*%Gamma.hat) < 10^(-8)){
    cov.hat = ((W.hat[auxilliary.dim, auxilliary.dim])^(-1))/N
  }else{
    cov.hat = cov.hat = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)/N
    GammaW = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)%*%t(Gamma.hat)%*%W.hat
  }
  se = sqrt(diag(as.matrix(cov.hat)))

  return(list(sol.path = theta.m, theta.hat = theta.hat, cov.hat = cov.hat,
              se = se, lower = theta.hat-qnorm(0.975)*se, upper = theta.hat+qnorm(0.975)*se,
              g.m = g.m, GammaW = GammaW, h = t(cbind(d, x2, y))))
}

ChuangChao2023.cGMM = function(auxilliary, pi.fit.list){
  J = length(pi.fit.list)
  alpha.list = lapply(pi.fit.list, function(pi.fit) pi.fit$theta.hat[-length(pi.fit$theta.hat)])
  mu.hat.v = unlist(lapply(pi.fit.list, function(pi.fit) pi.fit$theta.hat[length(pi.fit$theta.hat)]))
  pi.m = matrix(NA, n, J)
  for(j in 1:J){
    pi.m[, j] = pi.fit.list[[j]]$model(y[r == 1], as.matrix(dat[r == 1, pi.fit.list[[j]]$model.x.names]), alpha.list[[j]], n)
  }

  w.pi.fit = find.w.pi(auxilliary, pi.m, r, init = rep(1/J, J))
  w.pi = w.pi.fit$theta.hat
  w.pi = as.matrix((w.pi^2)/as.numeric(t(w.pi)%*%w.pi))
  mu.hat = mu.hat.v%*%w.pi

  return(mu.hat)
}

collect.outcome = function(outcome.list){
  m1.m = NULL
  #------------------------------------------------------------------------------#
  # Collect the outcome models
  #------------------------------------------------------------------------------#
  K = length(outcome.list)
  m1.m = matrix(NA, N, K)
  for(k in 1:K){
    m1.m[, k] = predict(outcome.list[[k]], newdata = dat, type = "response")
  }
  #------------------------------------------------------------------------------#
  return(list(outcome.list = outcome.list, m1.m = m1.m))
}

collect.propensity = function(pi.list, pi.fit.list, alpha.list){
  n = sum(r)
  pi.m = NULL
  #------------------------------------------------------------------------------#
  # Collect the propensity score models
  #------------------------------------------------------------------------------#
  J = length(pi.list)
  pi.m = matrix(NA, n, J)
  for(j in 1:J){
    pi.m[, j] = pi.list[[j]](y[r == 1], as.matrix(dat[r == 1, pi.fit.list[[j]]$model.x.names]), alpha.list[[j]], n)
  }
  #------------------------------------------------------------------------------#
  return(list(pi.list = pi.list, pi.fit.list = pi.fit.list, alpha.list = alpha.list, pi.m = pi.m))
}

get.pi.m = function(pi.fit.list, pi.list, alpha.list, nonrespondent = TRUE){
  J = length(pi.fit.list)
  id = which(r == ifelse(nonrespondent, 0, 1))
  n = length(id)
  pi.m = matrix(NA, n, J)
  for(j in 1:J){
    pi.m[, j] = pi.list[[j]](y[id],  as.matrix(dat[id, pi.fit.list[[j]]$model.x.names]), alpha.list[[j]], n)
  }
  return(pi.m)
}

get.pi.m2 = function(pi.fit.list, pi.list, alpha.list, y, u1, u2, r, nonrespondent = TRUE){
  J = length(pi.fit.list)
  id = which(r == ifelse(nonrespondent, 0, 1))
  n = length(id)
  pi.m = matrix(NA, n, J)
  for(j in 1:J){
    pi.j = pi.list[[j]](y[id],  as.matrix(dat[id, pi.fit.list[[j]]$model.x.names]), alpha.list[[j]], n)
    pi.m[, j] = apply(apply(pi.j, 1, function(pi.j.i) rbinom(N, size = 1, prob = pi.j.i)), 2, mean)
  }
  return(pi.m)
}

est.AIPW.se = function(h){
  augment.alpha = apply(-r*h/(imputed.pi^2)*(t(t(pi.x.alpha.m)*rep(w.pi, dim.alpha))), 2, mean)+apply((r/imputed.pi-1)*G.x, 2, mean)
  augment.beta = apply(-r*h/(imputed.pi^2)*(imputed.pi.m%*%norm.w.pi.adj), 2, mean)+apply((r/imputed.pi-1)*B.x, 2, mean)
  E.Q = (ipw.alpha-augment.alpha)+(ipw.beta-augment.beta)%*%beta.pi.prime
  E.T = ipw.beta+augment.beta

  mu.AIPW.iid = as.vector(r/imputed.pi*y-(r/imputed.pi-1)*h+E.Q%*%alpha.iid+E.T%*%beta.pi.iid)
  se.AIPW = sqrt(var(mu.AIPW.iid)/N)

  mu.AIPW.true.iid = as.vector(r/imputed.pi*y-(r/imputed.pi-1)*h)
  se.AIPW.true = sqrt(var(mu.AIPW.true.iid)/N)

  return(list(se.AIPW = se.AIPW, se.AIPW.true = se.AIPW.true))
}

ChuangChao2023 = function(pi.fit.list, auxilliary, family = "gaussian", ortho = FALSE, true.pi = NULL, perturb = rep(1, N)){
  ################################################################################
  # Collect the propensity score models
  ################################################################################
  pi.list = lapply(pi.fit.list, function(pi.fit) pi.fit$model)
  alpha.list = lapply(pi.fit.list, function(pi.fit) pi.fit$theta.hat)
  propensity.collected = collect.propensity(pi.list, pi.fit.list, alpha.list)
  pi.list = propensity.collected$pi.list
  pi.fit.list = propensity.collected$pi.fit.list
  alpha.list = propensity.collected$alpha.list
  pi.m = propensity.collected$pi.m
  J = length(pi.list)
  ################################################################################

  dim.alpha = unlist(lapply(alpha.list, length))

  ################################################################################
  # Compress the propensity score models
  ################################################################################
  w.pi.fit = find.w.pi(auxilliary, pi.m, r, init = rep(1/J, J), perturb = perturb)
  w.pi = w.pi.fit$theta.hat
  w.pi = as.matrix((w.pi^2)/as.numeric(t(w.pi)%*%w.pi))
  compressed.pi = rep(1, N)
  compressed.pi[r == 1] = pi.m%*%w.pi
  ################################################################################

  ################################################################################
  # calculate pi.x.alpha.m, beta.pi.prime ... (for estimating the ASE)
  ################################################################################
  dim.alpha = unlist(lapply(alpha.list, length))

  imputed.pi.m = matrix(1, N, J)
  imputed.pi.m[r == 1,] =  get.pi.m(pi.fit.list, pi.list, alpha.list, nonrespondent = FALSE)

  imputed.pi = rep(1, N)
  imputed.pi[r == 1] = imputed.pi.m[r == 1,]%*%w.pi

  pi.x.alpha.m = matrix(0, N, sum(dim.alpha))
  for(j in 1:J){
    pi.x.alpha.m[r == 1, (sum(dim.alpha[0:(j-1)])+1):sum(dim.alpha[1:j])] = jacobian(function(alpha.v) pi.list[[j]](y[r == 1],  as.matrix(dat[r == 1, pi.fit.list[[j]]$model.x.names]), alpha.v, n), alpha.list[[j]])
  }
  beta.pi.prime = -(w.pi.fit$GammaW%*%w.pi.fit$h%*%(-as.vector(r*((imputed.pi)^(-2)))*t(t(pi.x.alpha.m)*rep(w.pi, dim.alpha))))/N

  norm.adj = function(w){
    (diag(2*w)*sum(w^2)-2*(w^2)%*%t(w))/(sum(w^2)^2)
  }

  # norm.adj = function(w){
  #   (diag(2*w)*sum(w^2)-2*(w)%*%t(w^2))/(sum(w^2)^2)
  # }

  norm.w.pi.adj = 0; if(length(w.pi) > 1) norm.w.pi.adj = norm.adj(w.pi.fit$theta.hat)

  alpha.iid =  (c(lapply(pi.fit.list, function(pi.fit) -pi.fit$GammaW)) %>% bdiag() %*%
                  do.call(rbind, lapply(pi.fit.list, function(pi.fit) pi.fit$g.m)) %>% as.matrix())[1:sum(dim.alpha),]

  beta.pi.iid = -w.pi.fit$GammaW%*%w.pi.fit$g.m
  ################################################################################

  ################################################################################
  # IPW estimator for the population mean mu
  ################################################################################
  mu.IPW = mean(perturb*r/compressed.pi*y)
  mu.IPW.true = NULL

  if(!is.null(true.pi)){
    mu.IPW.true = mean(r/true.pi*y)
  }
  ################################################################################

  ################################################################################
  # Estimate the asymptotic variance of IPW
  ################################################################################
  ipw.alpha = apply(-r*y/(imputed.pi^2)*(t(t(pi.x.alpha.m)*rep(w.pi, dim.alpha))), 2, mean)
  ipw.beta = apply(-r*y/(imputed.pi^2)*(imputed.pi.m%*%norm.w.pi.adj), 2, mean)
  E.H = ipw.alpha+ipw.beta%*%beta.pi.prime

  mu.IPW.iid = as.vector(r/imputed.pi*y+E.H%*%alpha.iid+ipw.beta%*%beta.pi.iid)
  se.IPW = sqrt(var(mu.IPW.iid)/N)

  mu.IPW.true.iid = as.vector(r/true.pi*y)
  se.IPW.true = sqrt(var(mu.IPW.true.iid)/N)
  ################################################################################

  return(list(mu.IPW = mu.IPW,
              mu.IPW.true = mu.IPW.true,
              se.IPW = se.IPW,
              se.IPW.true = se.IPW.true,
              pi.fit.list = pi.fit.list, pi.m = pi.m, dim.alpha = dim.alpha,
              nu = w.pi.fit$theta.hat,
              imbalance = sum(apply(w.pi.fit$g.m, 1, mean)^2),
              w.pi = w.pi))
}

mild.est = function(pi.m, true.index, exp.tilt, exp.tilt.x.names, auxilliary, family = "gaussian"){
  J = ncol(pi.m)
  pi.m[, true.index] = exp.tilt(y[r == 1], dat[r == 1, exp.tilt.x.names])*pi.m[, true.index]
  imputed.pi.m = matrix(1, N, J)
  imputed.pi.m[r == 1, ] = pi.m
  ################################################################################
  # Compress the propensity score models
  ################################################################################
  w.pi.fit = find.w.pi(auxilliary, pi.m, r, ortho = FALSE, init = rep(1/J, J))
  w.pi = w.pi.fit$theta.hat
  w.pi = as.matrix((w.pi^2)/as.numeric(t(w.pi)%*%w.pi))
  compressed.pi = rep(1, N)
  compressed.pi[r == 1] = pi.m%*%w.pi
  ################################################################################

  ################################################################################
  # calculate norm.w.pi.adj, beta.pi.iid ... (for estimating the ASE)
  ################################################################################
  norm.adj = function(w){
    (diag(2*w)*sum(w^2)-2*(w^2)%*%t(w))/(sum(w^2)^2)
  }
  norm.w.pi.adj = 0; if(length(w.pi) > 1) norm.w.pi.adj = norm.adj(w.pi.fit$theta.hat)

  beta.pi.iid = -w.pi.fit$GammaW%*%w.pi.fit$g.m
  ################################################################################

  ################################################################################
  # IPW estimator for the population mean mu
  ################################################################################
  mu.IPW = mean(r/compressed.pi*y)
  ################################################################################

  ################################################################################
  # Estimate the asymptotic variance of IPW
  ################################################################################
  ipw.beta = apply(-r*y/(compressed.pi^2)*(imputed.pi.m%*%norm.w.pi.adj), 2, mean)

  mu.IPW.iid = as.vector(r/compressed.pi*y+ipw.beta%*%beta.pi.iid)
  se.IPW = sqrt(var(mu.IPW.iid)/N)

  mu.IPW.true.iid = as.vector(r/compressed.pi*y)
  se.IPW.true = sqrt(var(mu.IPW.true.iid)/N)
  ################################################################################

  return(list(mu.IPW = mu.IPW,
              se.IPW = se.IPW,
              pi.m = pi.m,
              nu = w.pi.fit$theta.hat,
              w.pi = w.pi))
}

