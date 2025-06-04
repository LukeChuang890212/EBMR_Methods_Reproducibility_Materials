formula = ps_specifications$formula.list[[1]]
h_x_names = ps_specifications$h_x_names.list[[1]]
alpha_init = ps_specifications$alpha_init.list[[1]]
inv_link = ps_specifications$inv_link

data = subdat

# self$EBMR_IPW = ifelse(is.perturb, EBMR_IPW_perturb, EBMR_IPW)
r = data$r
y = data[y_names]
n = nrow(data)
J = length(ps_specifications$formula.list)
ps_fit.list =list()
for(j in 1:J){
  formula = ps_specifications$formula.list[[j]]
  h_x_names = ps_specifications$h_x_names.list[[j]]
  alpha_init = ps_specifications$alpha_init.list[[j]]
  inv_link = ps_specifications$inv_link
  if(is.null(wt)){
    ps_fit.list[[j]] = WangShaoKim2014(formula, h_x_names, inv_link, alpha_init)
  }else{
    ps_fit.list[[j]] = WangShaoKim2014(formula, h_x_names, inv_link, alpha_init, se.fit = F, wt)
  }
}

WangShaoKim2014 = function(formula, h_x_names, inv_link, init = NULL, se.fit = T, wt = NULL) {
  # Basic setup
  result = parse_formula(formula)
  r = as.matrix(data[result$r_names])
  y = as.matrix(data[result$y_names])
  x = as.matrix(data[result$x_names])
  n = nrow(data)
  model_x_names = colnames(x)

  result = separate_variable_types(data[h_x_names])
  h_x1 = result$x1
  h_x2 = result$x2

  is.mnar = ifelse(ncol(y) == 0, FALSE, TRUE)
  alpha_dim = 1 + as.numeric(is.mnar) + ncol(x)

  d = NULL
  if(ncol(h_x1) > 0){
    for(j in 1:ncol(h_x1)) h_x1[, j] = as.factor(h_x1[, j])
    d = model.matrix(lm(rep(1, n)~., data = h_x1))
  }

  if(ncol(h_x2) == 0){
    h_x2 = NULL
  }else{
    h_x2 = as.matrix(h_x2)
  }

  h_x = cbind(d, h_x2)
  h_dim = ncol(h_x)

  model = function(x, y, alpha){
    if(!is.mnar) y = NULL
    inv_link(cbind(rep(1, n), y, x)%*%alpha)
  }

  wt <- if (is.null(wt)) 1 else wt

  # single, elegant definition of phi_alpha
  Phi_alpha <- function(param) {
    rw = r/model(x, y, param)
    g.matrix = as.vector(rw-1)*h_x
    return(wt*g.matrix)
  }
  # if(!is.null(wt)) Phi_alpha = function(param) wt*Phi_alpha(param)
  gmm_fit = gmm(Phi_alpha, W, n, h_dim, alpha_dim, init, se.fit)
  gmm_fit$se

  results = list(coefficients = gmm_fit$estimates,
                 se = gmm_fit$se,
                 fitted.values = model(x, y, gmm_fit$estimates),
                 model = model,
                 model_x_names = model_x_names,
                 h_x = h_x,
                 gmm_fit = gmm_fit)

  return(results)
}

g = Phi_alpha
param_dim = alpha_dim

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
    H_2n_3 = sapply(1:n, function(i) t(Gamma.hat)%*%(-W.hat%*%(g.matrix[i, ]%*%t(g.matrix[i, ]))%*%W.hat+W.hat)%*%eta_s)
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
