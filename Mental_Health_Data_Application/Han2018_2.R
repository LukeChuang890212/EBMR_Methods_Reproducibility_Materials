Han2018 = function(outcome.list, pi.fit.list, family = "gaussian"){
  N = length(r)
  n = sum(r == 1)
  K = length(outcome.list)
  J = length(pi.fit.list)

  m1.m = NULL
  outcome.list = outcome.list[1:(wrong.num+1)]
  K = length(outcome.list)
  m1.m = matrix(NA, N, K)
  for(k in 1:K){
    m1.m[, k] = predict(outcome.list[[k]], newdata = dat, type = "response")
  }

  pi.list = lapply(pi.fit.list, function(pi.fit) pi.fit$model)
  alpha.list = lapply(pi.fit.list, function(pi.fit) pi.fit$theta.hat)

  mu.hat = rep(NA, J)
  for(j in 1:J){
    g = matrix(NA, n, K)
    for(k in 1:K){
      L = 2000
      m0 = rep(NA, N-n)
      sigma = ifelse(family == "gaussian", sd(y[r == 1]-m1.m[r == 1, k]), NA)
      nonrespondent.id = which(r == 0)
      for(i in 1:(N-n)){
        dat.value = dat[nonrespondent.id[i],]
        # z1.value = as.numeric(z1[nonrespondent.id[i]])
        # z2.value = ifelse(is.null(z2[nonrespondent.id[i]]), NA, z2[nonrespondent.id[i]])
        # u1.value = ifelse(is.null(u1[nonrespondent.id[i]]), NA, u1[nonrespondent.id[i]])
        # u2.value = ifelse(is.null(u2[nonrespondent.id[i]]), NA, u2[nonrespondent.id[i]])

        m1.x = predict(outcome.list[[k]],
                       dat.value,
                       type = "response")
        y.s = NULL
        if(family == "gaussian"){
          y.s = rnorm(L, mean = m1.x, sd = sigma)
        }else if(family == "binomial"){
          y.s = rbinom(L, size = 1, prob = m1.x)
        }

        pi.x = pi.list[[j]](y.s, rep(1, L)%*%as.matrix(dat.value[pi.fit.list[[j]]$model.x.names]), alpha.list[[j]], L)
        a = mean((1-pi.x)/pi.x)
        m0[i] = mean(((1-pi.x)/pi.x)/a*(y.s))
      }
      g[, k] = m1.m[r == 1, k]-mean(na.omit(c(m1.m[r == 1, k], m0)))
    }

    Fn = function(rho){
      # value = -mean(log(1+g%*%rho))
      # return(ifelse(is.infinite(value), 10^8, value))
      return(-mean(log(1+g%*%rho)))
    }

    fn = function(rho){
      return(-apply(g, 2, function(gj) mean(gj/(1+g%*%rho))))
    }

    opt = constrOptim(rep(0, K), f = Fn, grad = fn, ui = g, ci = rep(-1, n))
    rho = opt$par
    w = c(1/(n*(1+g%*%rho)))
    mu.hat[j] = sum(w*y[r == 1])
  }

  return(list(mu.hat = mu.hat, w = w))
}

