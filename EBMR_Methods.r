library(R6)

MyClass <- R6Class("MyClass",
                   public = list(
                     # Public fields (variables)
                     h = NULL,
                     candidate_ps = NULL,
                     r = NULL,
                     init = NULL,
                     discrete_dim = 0,
                     continuous_dim = 0,
                     J = NULL,
                     nu_sol_path = NULL,
                     
                     # Constructor to initialize fields
                     initialize = function(field1 = NULL, field2 = NULL) {
                       self$field1 <- field1
                       self$J <- field2
                     },
                     
                     # Public methods (functions)
                     method1 = function() {
                       print("Method 1 called!")
                     },
                     
                     method2 = function(x) {
                       return(self$field1 + x)
                     }
                   )
)

estimate_nu = function(h, candidate_ps, r, init){
  x1 = h[[1]]; x2 = h[[2]];
  if(!is.null(x1)) x1 = as.matrix(x1)
  if(!is.null(x2)) x2 = as.matrix(x2)

  N = length(r)
  discrete_dim = 0
  if(!is.null(x1)) for(j in 1:ncol(x1)) discrete_dim = discrete_dim + length(unique(x1[,j]))
  discrete_dim = ifelse(discrete_dim == 0, 0, discrete_dim - ncol(x1) + 1) # to avoid collinearity
  continuous_dim = ifelse(!is.null(x2), ncol(x2), 0)
  h_dim = discrete_dim + continuous_dim

  J = ncol(candidate_ps)
  if(is.null(init)) init = rep(0, J)

  d = NULL
  if(!is.null(x1)){
    x1 =  as.data.frame(x1)
    for(j in 1:ncol(x1)) x1[, j] = as.factor(x1[, j])
    d = model.matrix(lm(rep(1, N)~., data =  x1))
  }
  
  g = function(nu){
    g.matrix = matrix(NA, h_dim, N)
    rw = rep(0, N)
    rw[r == 1] = 1/(candidate_ps%*%nu)
    if(discrete_dim > 0){
      for(l in 1:discrete_dim){g.matrix[l,] = t(d)[l,]*(rw-1)}
    }
    if(continuous_dim > 0){
      for(l in (discrete_dim+1):(discrete_dim+continuous_dim)) g.matrix[l,] = t(x2)[l-discrete_dim,]*(rw-1)
    }
    return(g.matrix)
  }

  G = function(g.matrix){
    return(matrix(apply(g.matrix, 1, mean), h_dim, 1))
  }

  W = function(g.matrix){
    return(solve(g.matrix%*%t(g.matrix)/N))
  }

  Gamma = function(nu){
    gamma = array(NA, dim = c(h_dim, J, N))
    for(i in 1:h_dim){
      gamma[i,,] = t(jacobian(function(nu.v) g(nu.v)[i,], nu))
    }
    return(apply(gamma, c(1, 2), mean))
  }

  obj = function(nu){
    g.matrix = g(nu)
    G.hat = G(g.matrix)
    value = t(G.hat)%*%G.hat
    return(ifelse(is.infinite(value) || is.na(value), 10^8, value))
  }

  nu_sol_path = matrix(init, J)
  conv_err = 10^8
  t = 1

  while (conv_err > 10^(-8) & t < 1000){
    opt = optim(nu_sol_path[, t], obj, method = "L-BFGS-B",
                lower = rep(-Inf, ncol(candidate_ps)), upper = rep(Inf, ncol(candidate_ps)))
    nu_sol_path = cbind(nu_sol_path, opt$par)
    g.matrix = g(nu_sol_path[,t+1]); W.hat = W(g.matrix);
    obj = function(nu){
      g.matrix = g(nu)
      G.hat = G(g.matrix)
      value = t(G.hat)%*%W.hat%*%G.hat
      return(ifelse(is.infinite(value) || is.na(value), 10^8, value))
    }
    conv_err = max(abs(nu_sol_path[,t+1]-nu_sol_path[,t]))
    t = t + 1
  }

  nu.hat = nu_sol_path[, t]
  Gamma.hat = Gamma(nu.hat)
  g.matrix = g(nu.hat); W.hat = W(g);
  S = var(t(g.matrix))
  # cov.hat = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)/N
  cov.hat = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)%*%t(Gamma.hat)%*%W.hat%*%S%*%W.hat%*%Gamma.hat%*%solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)/N
  K_nu = eta%*%solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)%*%t(Gamma.hat)%*%W.hat
  se = sqrt(diag(cov.hat))

  return(list(sol.path = nu_sol_path, nu.hat = nu.hat, cov.hat = cov.hat,
              se = se, lower = nu.hat-qnorm(0.975)*se, upper = nu.hat+qnorm(0.975)*se,
              g.matrix = g.matrix, K_nu = K_nu, h = t(cbind(d, x2))))
}

collect.propensity = function(pi.list, pi.fit.list, alpha.list){
  n = sum(r)
  candidate_ps = NULL
  #------------------------------------------------------------------------------#
  # Collect the propensity score models
  #------------------------------------------------------------------------------#
  J = length(pi.list)
  candidate_ps = matrix(NA, n, J)
  for(j in 1:J){
    candidate_ps[, j] = pi.list[[j]](y[r == 1], as.matrix(dat[r == 1, pi.fit.list[[j]]$model.x.names]), alpha.list[[j]], n)
  }
  #------------------------------------------------------------------------------#
  return(list(pi.list = pi.list, pi.fit.list = pi.fit.list, alpha.list = alpha.list, candidate_ps = candidate_ps))
}

get.candidate_ps = function(pi.fit.list, pi.list, alpha.list, nonrespondent = TRUE){
  J = length(pi.fit.list)
  id = which(r == ifelse(nonrespondent, 0, 1))
  n = length(id)
  candidate_ps = matrix(NA, n, J)
  for(j in 1:J){
    candidate_ps[, j] = pi.list[[j]](y[id],  as.matrix(dat[id, pi.fit.list[[j]]$model.x.names]), alpha.list[[j]], n)
  }
  return(candidate_ps)
}

EBMR_IPW = function(pi.fit.list, h, family = "gaussian", ortho = FALSE, true.pi = NULL, perturb = rep(1, N)){
  ################################################################################
  # Collect the propensity score models
  ################################################################################
  pi.list = lapply(pi.fit.list, function(pi.fit) pi.fit$model)
  alpha.list = lapply(pi.fit.list, function(pi.fit) pi.fit$nu.hat)
  propensity.collected = collect.propensity(pi.list, pi.fit.list, alpha.list)
  pi.list = propensity.collected$pi.list
  pi.fit.list = propensity.collected$pi.fit.list
  alpha.list = propensity.collected$alpha.list
  candidate_ps = propensity.collected$candidate_ps
  J = length(pi.list)
  ################################################################################

  dim.alpha = unlist(lapply(alpha.list, length))

  ################################################################################
  # Compress the propensity score models
  ################################################################################
  nu.fit = estimate_nu(h, candidate_ps, r, init = rep(1/J, J), perturb = perturb)
  nu = nu.fit$nu.hat
  nu = as.matrix((nu^2)/as.numeric(t(nu)%*%nu))
  compressed.pi = rep(1, N)
  compressed.pi[r == 1] = candidate_ps%*%nu
  ################################################################################

  ################################################################################
  # calculate pi.x.alpha.m, beta.pi.prime ... (for estimating the ASE)
  ################################################################################
  dim.alpha = unlist(lapply(alpha.list, length))

  imputed.candidate_ps = matrix(1, N, J)
  imputed.candidate_ps[r == 1,] =  get.candidate_ps(pi.fit.list, pi.list, alpha.list, nonrespondent = FALSE)

  imputed.pi = rep(1, N)
  imputed.pi[r == 1] = imputed.candidate_ps[r == 1,]%*%nu

  pi.x.alpha.m = matrix(0, N, sum(dim.alpha))
  for(j in 1:J){
    pi.x.alpha.m[r == 1, (sum(dim.alpha[0:(j-1)])+1):sum(dim.alpha[1:j])] = jacobian(function(alpha.v) pi.list[[j]](y[r == 1],  as.matrix(dat[r == 1, pi.fit.list[[j]]$model.x.names]), alpha.v, n), alpha.list[[j]])
  }
  beta.pi.prime = -(nu.fit$K_nu%*%nu.fit$h%*%(-as.vector(r*((imputed.pi)^(-2)))*t(t(pi.x.alpha.m)*rep(nu, dim.alpha))))/N

  norm.adj = function(w){
    (diag(2*w)*sum(w^2)-2*(w^2)%*%t(w))/(sum(w^2)^2)
  }

  # norm.adj = function(w){
  #   (diag(2*w)*sum(w^2)-2*(w)%*%t(w^2))/(sum(w^2)^2)
  # }

  norm.nu.adj = 0; if(length(nu) > 1) norm.nu.adj = norm.adj(nu.fit$nu.hat)

  alpha.iid =  (c(lapply(pi.fit.list, function(pi.fit) -pi.fit$K)) %>% bdiag() %*%
                  do.call(rbind, lapply(pi.fit.list, function(pi.fit) pi.fit$g)) %>% as.matrix())[1:sum(dim.alpha),]

  beta.pi.iid = -nu.fit$K_nu%*%nu.fit$g
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
  ipw.alpha = apply(-r*y/(imputed.pi^2)*(t(t(pi.x.alpha.m)*rep(nu, dim.alpha))), 2, mean)
  ipw.beta = apply(-r*y/(imputed.pi^2)*(imputed.candidate_ps%*%norm.nu.adj), 2, mean)
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
              pi.fit.list = pi.fit.list, candidate_ps = candidate_ps, dim.alpha = dim.alpha,
              nu = nu.fit$nu.hat,
              imbalance = sum(apply(nu.fit$g, 1, mean)^2),
              nu = nu))
}

EBMR_IPW_with_locally_misspecified_model = function(candidate_ps, perturb.index, exp.tilt, exp.tilt.x.names, h){
  J = ncol(candidate_ps)
  candidate_ps[, perturb.index] = exp.tilt(y[r == 1], dat[r == 1, exp.tilt.x.names])*candidate_ps[, perturb.index]
  imputed.candidate_ps = matrix(1, N, J)
  imputed.candidate_ps[r == 1, ] = candidate_ps
  ################################################################################
  # Compress the propensity score models
  ################################################################################
  nu.fit = estimate_nu(h, candidate_ps, r, ortho = FALSE, init = rep(1/J, J))
  nu = nu.fit$nu.hat
  nu = as.matrix((nu^2)/as.numeric(t(nu)%*%nu))
  compressed.pi = rep(1, N)
  compressed.pi[r == 1] = candidate_ps%*%nu
  ################################################################################
  
  ################################################################################
  # calculate norm.nu.adj, beta.pi.iid ... (for estimating the ASE)
  ################################################################################
  norm.adj = function(w){
    (diag(2*w)*sum(w^2)-2*(w^2)%*%t(w))/(sum(w^2)^2)
  }
  norm.nu.adj = 0; if(length(nu) > 1) norm.nu.adj = norm.adj(nu.fit$nu.hat)
  
  beta.pi.iid = -nu.fit$K_nu%*%nu.fit$g
  ################################################################################
  
  ################################################################################
  # IPW estimator for the population mean mu
  ################################################################################
  mu.IPW = mean(r/compressed.pi*y)
  ################################################################################
  
  ################################################################################
  # Estimate the asymptotic variance of IPW
  ################################################################################
  ipw.beta = apply(-r*y/(compressed.pi^2)*(imputed.candidate_ps%*%norm.nu.adj), 2, mean)
  
  mu.IPW.iid = as.vector(r/compressed.pi*y+ipw.beta%*%beta.pi.iid)
  se.IPW = sqrt(var(mu.IPW.iid)/N)
  
  mu.IPW.true.iid = as.vector(r/compressed.pi*y)
  se.IPW.true = sqrt(var(mu.IPW.true.iid)/N)
  ################################################################################
  
  return(list(mu.IPW = mu.IPW,
              se.IPW = se.IPW,
              candidate_ps = candidate_ps,
              nu = nu.fit$nu.hat,
              nu = nu))
}

