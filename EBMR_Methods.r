library(Matrix)
library(R6)

MyClass <- R6Class("MyClass",
                   public = list(
                     # Public fields (variables)
                     h = NULL,
                     ps.matrix = NULL,
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

estimate_nu = function(ps.matrix, h_x, init = NULL){
  result = separate_variable_types(h_x)
  h_x1 = result$x1
  h_x2 = result$x2
  
  J = ncol(ps.matrix)
  
  d = NULL
  if(ncol(h_x1) > 0){
    for(j in 1:ncol(h_x1)) h_x1[, j] = as.factor(h_x1[, j])
    d = model.matrix(lm(rep(1, n)~., data =  h_x1))
  }
  
  discrete_dim = ncol(d)
  continuous_dim = ncol(h_x2)
  h_dim = discrete_dim + continuous_dim
  
  r = as.matrix(r)
  y = as.matrix(y)
  x = as.matrix(x)
  h_x1 = as.matrix(h_x1)
  h_x2 = as.matrix(h_x2)
  
  g = function(nu){
    g.matrix = matrix(NA, n, h_dim)
    rw = r/(ps.matrix%*%nu)
    if(discrete_dim > 0){
      for(l in 1:discrete_dim){g.matrix[, l] = d[, l]*(rw-1)}
    }
    if(continuous_dim > 0){
      for(l in (discrete_dim+1):(discrete_dim+continuous_dim)) g.matrix[, l] = h_x2[, l-discrete_dim]*(rw-1)
    }
    return(g.matrix)
  }

  G = function(g.matrix){
    return(matrix(apply(g.matrix, 2, mean), h_dim, 1))
  }
  
  W = function(g.matrix){
    return(solve(t(g.matrix)%*%g.matrix/n))
  }
  
  Gamma = function(nu){
    Gamma.arr = array(NA, dim = c(h_dim, n, J))
    for(l in 1:h_dim){
      Gamma.arr[l,,] = jacobian(function(nu) g(nu)[, l], nu)
    }
    return(apply(Gamma.arr, c(1, 3), mean))
  }

  obj = function(nu){
    g.matrix = g(nu)
    G.hat = G(g.matrix)
    value = t(G.hat)%*%G.hat
    return(ifelse(is.infinite(value) || is.na(value), 10^8, value))
  }
  
  if(is.null(init)) init = rep(0, J)
  nu_sol_path = matrix(init, J)
  conv_err = 10^8
  t = 1

  while (conv_err > 10^(-8) & t < 1000){
    opt = optim(nu_sol_path[, t], obj, method = "L-BFGS-B",
                lower = rep(-Inf, ncol(ps.matrix)), upper = rep(Inf, ncol(ps.matrix)))
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
  g.matrix = g(nu.hat)
  W.hat = W(g.matrix)
  S = var(g.matrix)
  # cov.hat = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)/N
  cov.hat = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)%*%t(Gamma.hat)%*%W.hat%*%S%*%W.hat%*%Gamma.hat%*%solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)/n
  se = sqrt(diag(cov.hat))
  
  results = list(coefficients = nu.hat, 
                 sol.path = nu_sol_path, 
                 cov.hat = cov.hat,
                 se = se, 
                 lower = nu.hat-qnorm(0.975)*se, 
                 upper = nu.hat+qnorm(0.975)*se,
                 g.matrix = g.matrix, 
                 K = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)%*%t(Gamma.hat)%*%W.hat,
                 h_x = cbind(d, h_x2))
  
  return(results)
}

collect.propensity = function(ps_model.list, ps_fit.list, alpha.list){
  n = sum(r)
  ps.matrix = NULL
  #------------------------------------------------------------------------------#
  # Collect the propensity score models
  #------------------------------------------------------------------------------#
  J = length(ps_model.list)
  ps.matrix = matrix(NA, n, J)
  for(j in 1:J){
    ps.matrix[, j] = ps_model.list[[j]](y[r == 1], as.matrix(data[r == 1, ps_fit.list[[j]]$model.x.names]), alpha.list[[j]], n)
  }
  #------------------------------------------------------------------------------#
  return(list(ps_model.list = ps_model.list, ps_fit.list = ps_fit.list, alpha.list = alpha.list, ps.matrix = ps.matrix))
}

get.ps.matrix = function(ps_fit.list, ps_model.list, alpha.list, nonrespondent = TRUE){
  J = length(ps_fit.list)
  id = which(r == ifelse(nonrespondent, 0, 1))
  n = length(id)
  ps.matrix = matrix(NA, n, J)
  for(j in 1:J){
    ps.matrix[, j] = ps_model.list[[j]](y[id],  as.matrix(data[id, ps_fit.list[[j]]$model.x.names]), alpha.list[[j]], n)
  }
  return(ps.matrix)
}

EBMR_IPW = function(ps_fit.list, h_x, true_ps = NULL){
  ################################################################################
  # Collect the propensity score models
  ################################################################################
  J = length(ps_fit.list)
  ps_model.list = lapply(ps_fit.list, function(ps_fit) ps_fit$model)
  alpha.list = lapply(ps_fit.list, function(ps_fit) ps_fit$coefficients)
  alpha_dim = unlist(lapply(alpha.list, length))
  ps.matrix = do.call(cbind, lapply(ps_fit.list, function(ps_fit) ps_fit$fitted.values))
  ################################################################################


  ################################################################################
  # Compress the propensity score models
  ################################################################################
  ensemble_fit = estimate_nu(ps.matrix, h_x, init = rep(1/J, J))
  nu.hat = ensemble_fit$coefficients
  w.hat = nu.hat^2/sum(nu.hat^2)
  ensemble_ps = ps.matrix%*%w.hat
  ################################################################################

  ################################################################################
  # calculate dot_pi, beta.pi.prime ... (for estimating the ASE)
  ################################################################################
  # imputed.ps.matrix = matrix(1, N, J)
  # imputed.ps.matrix[r == 1,] =  get.ps.matrix(ps_fit.list, ps_model.list, alpha.list, nonrespondent = FALSE)
  # 
  # imputed.pi = rep(1, N)
  # imputed.pi[r == 1] = imputed.ps.matrix[r == 1,]%*%nu
    
  dot_pi = matrix(NA, n, sum(alpha_dim))
  for(j in 1:J){
    x = as.matrix(data[ps_fit.list[[j]]$model.x.names])
    dot_pi[, (sum(alpha_dim[0:(j-1)])+1):sum(alpha_dim[1:j])] = jacobian(function(alpha) ps_model.list[[j]](x, y, alpha), alpha.list[[j]])
  }
  E_dot_g = -(t(dot_pi)*rep(w.hat, alpha_dim))%*%(ensemble_fit$h_x*as.vector(r*((ensemble_ps)^(-2))))/n
  
  dot_W = function(w){
    (diag(2*w)*sum(w^2)-2*(w)%*%t(w^2))/(sum(w^2)^2)
  }
  dot_W_nu_hat = 0; if(length(nu.hat) > 1) dot_W_nu_hat = dot_W(nu.hat)
  
  H_alpha.w = apply(t((t(dot_pi)*rep(w.hat, alpha_dim)))*as.vector(r*y*((ensemble_ps)^(-2))), 2, mean)
  w.H_nu = apply(ps.matrix%*%t(dot_W_nu_hat)*as.vector(r*y*((ensemble_ps)^(-2))), 2, mean)
  
  K_alpha = c(lapply(ps_fit.list, function(ps_fit) ps_fit$K)) %>% bdiag()
  g_all = do.call(rbind, lapply(ps_fit.list, function(ps_fit) t(ps_fit$g.matrix)))
  K_nu = ensemble_fit$K
  g = t(ensemble_fit$g.matrix)
  
  # alpha.iid = -K_alpha%*%g_all 
  # 
  # nu.iid = -ensemble_fit$K%*%t(ensemble_fit$g.matrix)
  ################################################################################

  ################################################################################
  # IPW estimator for the population mean mu
  ################################################################################
  mu.IPW = mean(r/ensemble_ps*y)
  mu.IPW.true = NULL
  if(!is.null(true_ps)){
    mu.IPW.true = mean(r/true_ps*y)
  }
  ################################################################################

  ################################################################################
  # Estimate the asymptotic variance of IPW
  ################################################################################
  mu.IPW.iid = as.vector(t(r/ensemble_ps*y)
                         +(t(H_alpha.w)+t(w.H_nu)%*%K_nu%*%t(E_dot_g))%*%K_alpha%*%g_all
                         +t(w.H_nu)%*%K_nu%*%g)
  se.IPW = sqrt(var(mu.IPW.iid)/n)

  mu.IPW.true.iid = as.vector(r/true_ps*y)
  se.IPW.true = sqrt(var(mu.IPW.true.iid)/N)
  ################################################################################

  return(list(mu.IPW = mu.IPW,
              mu.IPW.true = mu.IPW.true,
              se.IPW = se.IPW,
              se.IPW.true = se.IPW.true,
              ps_fit.list = ps_fit.list, ps.matrix = ps.matrix, alpha_dim = alpha_dim,
              nu = ensemble_fit$nu.hat,
              imbalance = sum(apply(ensemble_fit$g, 1, mean)^2),
              nu = nu))
}

EBMR_IPW_with_locally_misspecified_model = function(ps.matrix, perturb.index, exp.tilt, exp.tilt.x.names, h){
  ################################################################################
  # Collect the propensity score models
  ################################################################################
  J = ncol(ps.matrix)
  ps.matrix[, perturb.index] = exp.tilt(y[r == 1], data[r == 1, exp.tilt.x.names])*ps.matrix[, perturb.index]
  imputed.ps.matrix = matrix(1, N, J)
  imputed.ps.matrix[r == 1, ] = ps.matrix
  ################################################################################
  
  ################################################################################
  # Compress the propensity score models
  ################################################################################
  ensemble_fit = estimate_nu(h, ps.matrix, r, ortho = FALSE, init = rep(1/J, J))
  nu = ensemble_fit$nu.hat
  nu = as.matrix((nu^2)/as.numeric(t(nu)%*%nu))
  compressed.pi = rep(1, N)
  compressed.pi[r == 1] = ps.matrix%*%nu
  ################################################################################
  
  ################################################################################
  # calculate dot_W_nu_hat, nu.iid ... (for estimating the ASE)
  ################################################################################
  dot_W = function(w){
    (diag(2*w)*sum(w^2)-2*(w^2)%*%t(w))/(sum(w^2)^2)
  }
  dot_W_nu_hat = 0; if(length(nu) > 1) dot_W_nu_hat = dot_W(ensemble_fit$nu.hat)
  
  nu.iid = -ensemble_fit$K_nu%*%ensemble_fit$g
  ################################################################################
  
  ################################################################################
  # IPW estimator for the population mean mu
  ################################################################################
  mu.IPW = mean(r/compressed.pi*y)
  ################################################################################
  
  ################################################################################
  # Estimate the asymptotic variance of IPW
  ################################################################################
  ipw.beta = apply(-r*y/(compressed.pi^2)*(imputed.ps.matrix%*%dot_W_nu_hat), 2, mean)
  
  mu.IPW.iid = as.vector(r/compressed.pi*y+ipw.beta%*%nu.iid)
  se.IPW = sqrt(var(mu.IPW.iid)/N)
  
  mu.IPW.true.iid = as.vector(r/compressed.pi*y)
  se.IPW.true = sqrt(var(mu.IPW.true.iid)/N)
  ################################################################################
  
  return(list(mu.IPW = mu.IPW,
              se.IPW = se.IPW,
              ps.matrix = ps.matrix,
              nu = ensemble_fit$nu.hat,
              nu = nu))
}

