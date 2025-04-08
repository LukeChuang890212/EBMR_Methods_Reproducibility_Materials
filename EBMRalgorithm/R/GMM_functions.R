# Phi_alpha = function(param){
#   rw = r/model(x, y, param)
#   g.matrix = as.vector(rw-1)*h_x
#   return(g.matrix)
# }
#
# Phi_nu = function(param){
#   rw = r/(ps.matrix%*%param)
#   g.matrix = as.vector(rw-1)*h_x
#   return(g.matrix)
# }

gmm = function(g, W, n, h_dim, param_dim, init, se.fit = T){
  G = function(param){
    g.matrix = g(param)
    return(matrix(apply(g.matrix, 2, mean), h_dim, 1))
  }

  t_Gamma_i = function(param){
    Gamma.arr = array(NA, dim = c(n, param_dim, h_dim))
    for(l in 1:h_dim){
      Gamma.arr[,,l] = jacobian(function(param) g(param)[, l], param)
    }
    return(Gamma.arr)
  }

  Gamma = function(param){
    return(t(apply(t_Gamma_i(param), c(2, 3), mean)))
  }

  Gamma_2 = function(param){
    return(jacobian(function(param) as.vector(t(Gamma(param))), param))
  }

  obj = function(param){
    g.matrix = g(param)
    G.hat = G(param)
    value = t(G.hat)%*%G.hat
    return(ifelse(is.infinite(value) || is.na(value), 10^8, value))
  }

  if(is.null(init)) init = rep(0, param_dim)
  sol_path = matrix(init, param_dim)
  conv_err = 10^8
  t = 1

  while (conv_err > 10^(-8) & t < 1000){
    opt = optim(sol_path[, t], obj, method = "L-BFGS-B",
                lower = rep(-Inf, param_dim), upper = rep(Inf, param_dim))
    sol_path = cbind(sol_path, opt$par)
    g.matrix = g(sol_path[,t+1]); W.hat = W(g.matrix);
    obj = function(param){
      g.matrix = g(param)
      G.hat = G(param)
      value = t(G.hat)%*%W.hat%*%G.hat
      return(ifelse(is.infinite(value) || is.na(value), 10^8, value))
    }
    conv_err = max(abs(sol_path[,t+1]-sol_path[, t]))
    t = t + 1
  }

  estimates = sol_path[, t]

  Gamma.hat = W.hat = g.matrix = eta_s = Q = h_x = psi = se = NA
  if(se.fit){
    Gamma.hat = Gamma(estimates)
    g.matrix = g(estimates)
    W.hat = W(g.matrix)

    eta_s = apply(g.matrix, 2, mean)
    M_n = kronecker(eta_s%*%W.hat, diag(param_dim))%*%Gamma_2(estimates)
    H_0n = -solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)
    H_1n = t(Gamma.hat)%*%W.hat%*%(t(g.matrix)-eta_s)
    H_2n_2 = apply(t_Gamma_i(estimates), 1, function(m) (m-t(Gamma.hat))%*%W.hat%*%eta_s)
    H_2n_3 = sapply(1:n, function(i) t(Gamma.hat)%*%(g.matrix[i, ]%*%t(g.matrix[i, ])-W.hat)%*%eta_s)
    Q = solve(diag(param_dim)-H_0n%*%M_n)%*%H_0n
    psi = Q%*%(H_1n+H_2n_2+H_2n_3)
    # psi_alpha = solve(diag(alpha_dim)-H_0n%*%M_n)%*%H_0n%*%(H_1n+H_2n_2+H_2n_3)
    cov.hat = var(t(psi))/n
    se = sqrt(diag(cov.hat))
  }

  result = list(
    estimates = estimates,
    se = se,
    Gamma.hat = Gamma.hat,
    W.hat = W.hat,
    g.matrix = g.matrix,
    eta_s = eta_s,
    Q = Q,
    h_x = h_x,
    psi = psi
  )

  return(result)
}
