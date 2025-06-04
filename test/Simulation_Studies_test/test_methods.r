parse_formula = function(formula) {
  # Convert the formula to a string
  formula <- as.character(formula)
  
  # Extract the left-hand side of the formula (everything before the ~)
  lhs <- formula[2]
  
  # Isolate r
  r_names <- lhs
  
  # Extract the right-hand side of the formula (everything after the ~)
  rhs <- formula[3]
  
  # Extract the variables inside o() using regular expressions
  y_names <- str_extract(rhs, "o\\(([^)]+)\\)")   # Variables in o()
  
  if(is.na(y_names)){
    y_names = character(0)
  }else{
    # Clean up: Remove 'o()' and split the variables (there's only one in this case)
    y_names <- gsub("o\\(|\\)", "", y_names)
    y_names <- unlist(strsplit(y_names, split = "\\+"))
  }
  
  
  # Extract the rest of the variables (excluding o())
  x_names <- gsub("o\\([^)]+\\)\\s*\\+?", "", rhs)
  
  # Clean up: Split the rest of the variables by '+' and trim whitespace
  x_names <- unlist(strsplit(x_names, split = "\\+"))
  x_names <- trimws(x_names)
  
  # Return a list with r_names, y_names, x_names
  return(list(r_names = r_names, y_names = y_names, x_names = x_names))
}

separate_variable_types = function(x) {
  # Initialize empty lists to store continuous and discrete columns
  x1 <- list()
  x2 <- list()
  
  col_names = colnames(x)
  x1_names = c()
  x2_names = c()
  
  # Iterate through each column in the matrix 'x'
  for (i in 1:ncol(x)) {
    column <- x[, i]
    
    # Check if the column is numeric
    if (is.numeric(column)) {
      # If the column has more than 10 unique values, it's treated as continuous
      if (length(unique(column)) > 10) {
        x2[[colnames(x)[i]]] <- column
        x2_names = c(x2_names, col_names[i])
      } else {
        x1[[colnames(x)[i]]] <- column
        x1_names = c(x1_names, col_names[i])
      }
    } else {
      # For non-numeric columns, treat them as discrete
      x1[[colnames(x)[i]]] <- column
      x1_names = c(x1_names, col_names[i])
    }
  }
  
  # Convert the lists to matrices (if you want matrices, otherwise keep them as lists)
  x1 <- as.data.frame(x1)
  x2 <- as.data.frame(x2)
  
  # Return the results
  return(list(x1 = x1, x2 = x2, x1_names = x1_names, x2_names = x2_names))
}

WangShaoKim2014 = function(formula, h_x_names, inv_link, W, data, init = NULL) {
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
    d = model.matrix(lm(rep(1, n)~., data =  h_x1))
  }

  discrete_dim = ncol(d)
  continuous_dim = ncol(h_x2)
  h_dim = discrete_dim + continuous_dim

  if(continuous_dim == 0){
    h_x2 = NULL
  }else{
    h_x2 = as.matrix(h_x2)
  }

  model = function(x, y, alpha){
    if(!is.mnar) y = NULL
    inv_link(cbind(rep(1, n), y, x)%*%alpha)
  }
  w = function(x, y, alpha) 1/model(x, y, alpha)

  g = function(alpha){
    g.matrix = matrix(NA, n, h_dim)
    rw = r*w(x, y, alpha)
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

  # W = function(g.matrix){
  #   return(solve(t(g.matrix)%*%g.matrix/n))
  # }

  t_Gamma_i = function(alpha){
    Gamma.arr = array(NA, dim = c(n, alpha_dim, h_dim))
    for(l in 1:h_dim){
      Gamma.arr[,,l] = jacobian(function(alpha) g(alpha)[, l], alpha)
    }
    return(Gamma.arr)
  }

  Gamma = function(alpha){
    return(t(apply(t_Gamma_i(alpha), c(2, 3), mean)))
  }

  Gamma_2 = function(alpha){
    return(pracma::jacobian(function(alpha) as.vector(t(Gamma(alpha))), alpha))
  }

  obj = function(alpha){
    g.matrix = g(alpha)
    G.hat = G(g.matrix)
    value = t(G.hat)%*%G.hat
    return(ifelse(is.infinite(value) || is.na(value), 10^8, value))
  }

  if(is.null(init)) init = rep(0, alpha_dim)
  alpha_sol_path = matrix(init, alpha_dim)
  conv_err = 10^8
  t = 1

  while (conv_err > 10^(-8) & t < 1000){
    opt = optim(alpha_sol_path[,t], obj, method = "L-BFGS-B")
    alpha_sol_path = cbind(alpha_sol_path, opt$par)
    g.matrix = g(alpha_sol_path[,t+1]); W.hat = W(g.matrix);
    obj = function(alpha){
      g.matrix = g(alpha); G.hat = G(g.matrix);
      value = t(G.hat)%*%W.hat%*%G.hat
      return(ifelse(is.infinite(value) || is.na(value), 10^8, value))
    }
    conv_err = max(abs(alpha_sol_path[,t+1]-alpha_sol_path[,t]))
    t = t + 1
  }

  alpha.hat = alpha_sol_path[, t]
  fitted_values = model(x, y, alpha.hat)
  Gamma.hat = Gamma(alpha.hat)
  g.matrix = g(alpha.hat)
  W.hat = W(g.matrix)

  mu_s = apply(g.matrix, 2, mean)
  M_n = kronecker(mu_s%*%W.hat, diag(alpha_dim))%*%Gamma_2(alpha.hat)
  H_0n = -solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)
  H_1n = t(Gamma.hat)%*%W.hat%*%(t(g.matrix)-mu_s)
  H_2n_2 = apply(t_Gamma_i(alpha.hat), 1, function(m) (m-t(Gamma.hat))%*%W.hat%*%mu_s)
  H_2n_3 = sapply(1:n, function(i) t(Gamma.hat)%*%(g.matrix[i, ]%*%t(g.matrix[i, ])-W.hat)%*%mu_s)
  psi_alpha = solve(diag(alpha_dim)-H_0n%*%M_n)%*%H_0n%*%(H_1n+H_2n_2+H_2n_3)
  # psi_alpha = solve(diag(alpha_dim)-H_0n%*%M_n)%*%H_0n%*%(H_1n+H_2n_2+H_2n_3)
  cov.hat = var(t(psi_alpha))/n
  S = var(g.matrix)
  # cov.hat = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)/N
  cov.hat = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)%*%t(Gamma.hat)%*%W.hat%*%S%*%W.hat%*%Gamma.hat%*%solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)/n
  se = sqrt(diag(cov.hat))

  results = list(coefficients = alpha.hat,
                 fitted.values = fitted_values,
                 sol.path = alpha_sol_path,
                 cov.hat = cov.hat,
                 se = se,
                 lower = alpha.hat-qnorm(0.975)*se,
                 upper = alpha.hat+qnorm(0.975)*se,
                 g.matrix = g.matrix,
                 K = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)%*%t(Gamma.hat)%*%W.hat,
                 psi_alpha = psi_alpha,
                 model = model,
                 model_x_names = model_x_names,
                 h_x = cbind(d, h_x2))

  return(results)
}

estimate_nu = function(alpha, alpha_dim, ps_model.list, h_x, W, data, init = NULL){
  r = as.matrix(data$r)
  y = as.matrix(data$y)
  n = nrow(data)
  J = length(ps_model.list )

  h_dim = ncol(h_x)

  ps.matrix = matrix(NA, n, J)
  for(j in 1:J){
    x = as.matrix(data[ps_fit.list[[j]]$model_x_names])
    ps.matrix[, j] = ps_model.list[[j]](x, y, alpha[(1+sum(alpha_dim[0:(j-1)])):sum(alpha_dim[1:j])])
  }

  g = function(nu){
    rw = r/(ps.matrix%*%nu)
    g.matrix = as.vector(rw-1)*h_x
    return(g.matrix)
  }

  G = function(g.matrix){
    return(matrix(apply(g.matrix, 2, mean), h_dim, 1))
  }

  # W = function(g.matrix){
  #   return(solve(t(g.matrix)%*%g.matrix/n))
  # }

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
                lower = rep(-Inf, J), upper = rep(Inf, J))
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

  return(nu.hat)
}

ensemble = function(ps_fit.list, h_x, W, data, init = NULL, se.fit) {
  alpha.list = lapply(ps_fit.list, function(ps_fit) ps_fit$coefficients)
  alpha.hat = unlist(alpha.list)
  alpha_dim = unlist(lapply(alpha.list, length))
  ps_model.list = lapply(ps_fit.list, function(ps_fit) ps_fit$model)
  ps.matrix = do.call(cbind, lapply(ps_fit.list, function(ps_fit) ps_fit$fitted.values))

  # Basic setup
  r = as.matrix(data$r)
  n = nrow(data)
  # J = ncol(ps.matrix)

  result = separate_variable_types(h_x)
  h_x1 = result$x1
  h_x2 = result$x2

  d = NULL
  if(ncol(h_x1) > 0){
    for(j in 1:ncol(h_x1)) h_x1[, j] = as.factor(h_x1[, j])
    d = model.matrix(lm(rep(1, n)~., data =  h_x1))
  }

  if(ncol(h_x2) == 0){
    h_x2 = NULL
  }else{
    h_x2 = as.matrix(h_x2)
  }

  h_x = cbind(d, h_x2)
  h_dim = ncol(h_x)

  g = function(nu){
    rw = r/(ps.matrix%*%nu)
    g.matrix = as.vector(rw-1)*h_x
    return(g.matrix)
  }

  G = function(g.matrix){
    return(matrix(apply(g.matrix, 2, mean), h_dim, 1))
  }

  # W = function(g.matrix){
  #   return(solve(t(g.matrix)%*%g.matrix/n))
  # }

  t_Gamma_i = function(nu){
    Gamma.arr = array(NA, dim = c(n, J, h_dim))
    for(l in 1:h_dim){
      Gamma.arr[,,l] = jacobian(function(nu) g(nu)[, l], nu)
    }
    return(Gamma.arr)
  }

  Gamma = function(nu){
    return(t(apply(t_Gamma_i(nu), c(2, 3), mean)))
  }

  Gamma_2 = function(nu){
    return(jacobian(function(nu) as.vector(t(Gamma(nu))), nu))
  }

  nu.hat = estimate_nu(alpha.hat, alpha_dim, ps_model.list, h_x, W, data, init)

  # dot_nu.hat = NA
  # if(se.fit) dot_nu.hat = pracma::jacobian(function(alpha) estimate_nu(alpha, alpha_dim, ps_model.list, h_x, data, init), alpha.hat)

  Gamma.hat = Gamma(nu.hat)
  g.matrix = g(nu.hat)
  W.hat = W(g.matrix)
  mu_s = apply(g.matrix, 2, mean)
  M_n = kronecker(mu_s%*%W.hat, diag(J))%*%Gamma_2(nu.hat)
  H_0n = -solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)
  H_1n = t(Gamma.hat)%*%W.hat%*%(t(g.matrix)-mu_s)
  H_2n_2 = apply(t_Gamma_i(nu.hat), 1, function(m) (m-t(Gamma.hat))%*%W.hat%*%mu_s)
  H_2n_3 = sapply(1:n, function(i) t(Gamma.hat)%*%(g.matrix[i, ]%*%t(g.matrix[i, ])-W.hat)%*%mu_s)
  psi_nu = solve(diag(J)-H_0n%*%M_n)%*%H_0n%*%(H_1n+H_2n_2+H_2n_3)
  cov.hat = var(t(psi_nu))/n
  # S = var(g.matrix)
  # # cov.hat = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)/N
  # cov.hat = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)%*%t(Gamma.hat)%*%W.hat%*%S%*%W.hat%*%Gamma.hat%*%solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)/n
  se = sqrt(diag(cov.hat))

  results = list(coefficients = nu.hat,
                 # sol.path = nu_sol_path,
                 cov.hat = cov.hat,
                 se = se,
                 lower = nu.hat-qnorm(0.975)*se,
                 upper = nu.hat+qnorm(0.975)*se,
                 g.matrix = g.matrix,
                 Gamma = Gamma.hat,
                 W = W.hat,
                 # dot_nu.hat = dot_nu.hat,
                 K = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)%*%t(Gamma.hat)%*%W.hat,
                 Q = solve(diag(J)-H_0n%*%M_n)%*%H_0n,
                 mu_s = mu_s,
                 psi_nu = psi_nu,
                 h_x = cbind(d, h_x2))

  return(results)
}

EBMR_IPW = function(h_x_names, W, data, se.fit = TRUE, true_ps = NULL) {
  # Basic setup
  r = as.matrix(data$r)
  y = as.matrix(data$y)
  n = nrow(data)

  ################################################################################
  # Collect the propensity score models
  ################################################################################
  J = length(ps_fit.list)
  ps_model.list = lapply(ps_fit.list, function(ps_fit) ps_fit$model)
  alpha.list = lapply(ps_fit.list, function(ps_fit) ps_fit$coefficients)
  # alpha = unlist(lapply(ps_fit.list, function(ps_fit) ps_fit$coefficients))
  # alpha_var.list = lapply(ps_fit.list, function(ps_fit) (ps_fit$se)^2)
  alpha_dim = unlist(lapply(alpha.list, length))
  ps.matrix = do.call(cbind, lapply(ps_fit.list, function(ps_fit) ps_fit$fitted.values))
  ################################################################################

  ################################################################################
  # Ensemble step
  ################################################################################
  ensemble_fit = ensemble(ps_fit.list, data[h_x_names], W, data, init = rep(1/J, J), se.fit)
  # ensemble_fit = estimate_nu(ps.matrix, data[h_x_names], data, init = rep(1/J, J))
  nu.hat = ensemble_fit$coefficients
  # W = function(v) v^2/sum(v^2)
  # V_inv = if (J == 1) {
  #   1/sum(alpha_var.list[[1]])
  # } else {
  #   diag(sapply(alpha_var.list, function(alpha_var) 1/sum(alpha_var)))
  # }
  # w.hat = W(W(V_inv%*%nu.hat))
  w.hat = nu.hat^2/sum(nu.hat^2)
  ensemble_ps = ps.matrix%*%w.hat
  ################################################################################

  ################################################################################
  # Compute necessary quantities to estimate the influence function:
  # \psi(\bm{\alpha}_*, \bm{\nu}_*)
  ################################################################################
  dot_pi = matrix(NA, n, sum(alpha_dim))
  for(j in 1:J){
    x = as.matrix(data[ps_fit.list[[j]]$model_x_names])
    dot_pi[, (sum(alpha_dim[0:(j-1)])+1):sum(alpha_dim[1:j])] = jacobian(function(alpha) ps_model.list[[j]](x, y, alpha), alpha.list[[j]])
  }
  # E_dot_g = -(t(dot_pi)*rep(nu.hat, alpha_dim))%*%(ensemble_fit$h_x*as.vector(r*((ensemble_ps)^(-2))))/n


  dot_W = function(nu){
    nu = as.vector(nu)
    (diag(2*nu)*sum(nu^2)-2*(nu)%*%t(nu^2))/(sum(nu^2)^2)
  }
  # dot_W_nu_hat = 0; if(length(nu.hat) > 1) dot_W_nu_hat = dot_W(W(V_inv%*%nu.hat))%*%dot_W(V_inv%*%nu.hat)%*%V_inv
  dot_W_nu_hat = 0; if(length(nu.hat) > 1) dot_W_nu_hat = dot_W(nu.hat)

  H_alpha.w = apply(t((t(dot_pi)*rep(w.hat, alpha_dim)))*as.vector(r*y*((ensemble_ps)^(-2))), 2, mean)
  w.H_nu = apply(ps.matrix%*%t(dot_W_nu_hat)*as.vector(r*y*((ensemble_ps)^(-2))), 2, mean)

  # K_alpha = c(lapply(ps_fit.list, function(ps_fit) ps_fit$K)) %>% bdiag()
  # g_all = do.call(rbind, lapply(ps_fit.list, function(ps_fit) t(ps_fit$g.matrix)))
  psi_alpha = do.call(rbind, lapply(ps_fit.list, function(ps_fit) ps_fit$psi_alpha))
  # K_nu = ensemble_fit$K
  # g = t(ensemble_fit$g.matrix)
  psi_nu = ensemble_fit$psi_nu

  Gamma_nu = ensemble_fit$Gamma
  W_nu = ensemble_fit$W
  f_n = ensemble_fit$g.matrix
  h_nu = ensemble_fit$h_x
  mu_s = ensemble_fit$mu_s
  K_1 = t(Gamma_nu)%*%W_nu%*%t((t(dot_pi)*rep(nu.hat, alpha_dim))%*%(h_nu*as.vector(-r*((ps.matrix%*%nu.hat)^(-2))))/n)
  K_2_1 = (t(dot_pi)*rep(nu.hat, alpha_dim))%*%(as.vector(2*((ps.matrix%*%nu.hat)^(-3)))*as.vector(r*h_nu%*%W_nu%*%mu_s)*ps.matrix)/n
  K_2_2 = apply(as.vector((ps.matrix%*%nu.hat)^(-2))*as.vector(r*h_nu%*%W_nu%*%mu_s)*dot_pi, 2, mean)
  K_2_2 = bdiag(lapply(1:length(alpha_dim), function(j) matrix(K_2_2[(1+sum(alpha_dim[0:(j-1)])):sum(alpha_dim[0:(j)])], ncol = 1)))
  K_2 = t(K_2_1 - as.matrix(K_2_2))
  K_3 = t(Gamma_nu)%*%t((t(dot_pi)*rep(nu.hat, alpha_dim))%*%(h_nu*as.vector(-r*((ps.matrix%*%nu.hat)^(-2)))*as.vector(f_n%*%mu_s))/n)
  ################################################################################

  ################################################################################
  # IPW estimator for the population mean mu_0 with propensity score being estimated
  # by the methods of Wang, Shao and Kim (2014).
  ################################################################################
  mu_ipw = mean(r/ensemble_ps*y)
  # mu_ipw.iid = as.vector(t(r/ensemble_ps*y)
  #                        +(t(H_alpha.w)-t(w.H_nu)%*%K_nu%*%t(E_dot_g))%*%K_alpha%*%g_all
  #                        +t(w.H_nu)%*%K_nu%*%g)
  # mu_ipw.iid = as.vector(t(r/ensemble_ps*y)
  #                        +(t(H_alpha.w)-2*t(w.H_nu)%*%K_nu%*%t(E_dot_g))%*%K_alpha%*%g_all
  #                        +t(w.H_nu)%*%K_nu%*%g)
  se_ipw = NA
  if(se.fit){
    mu_ipw.iid = as.vector(t(r/ensemble_ps*y)
                           -(t(H_alpha.w)+t(w.H_nu)%*%ensemble_fit$Q%*%(K_1+K_2+K_3))%*%psi_alpha
                           -t(w.H_nu)%*%psi_nu)
    se_ipw = sqrt(var(mu_ipw.iid)/n)
  }
  ################################################################################

  ################################################################################
  # IPW estimator for the population mean mu_0 with known propensity score.
  ################################################################################
  mu_ipw.true = NA
  se_ipw.true = NA
  if(!is.null(true_ps)){
    mu_ipw.true = mean(r/true_ps*y)
    mu_ipw.true.iid = as.vector(r/true_ps*y)
    se_ipw.true = sqrt(var(mu_ipw.true.iid)/n)
  }
  ################################################################################

  result = list(mu_ipw = mu_ipw,
                mu_ipw.true = mu_ipw.true,
                se_ipw = se_ipw,
                se_ipw.true = se_ipw.true,
                ps.matrix = ps.matrix,
                nu.hat = nu.hat,
                w.hat = w.hat,
                imbalance = sum(apply(ensemble_fit$g, 1, mean)^2)
  )

  return(result)
}

WangShaoKim2014_perturb = function(formula, h_x_names, inv_link, W, data, wt, init = NULL, se.fit = TRUE) {
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
    d = model.matrix(lm(rep(1, n)~., data =  h_x1))
  }

  discrete_dim = ncol(d)
  continuous_dim = ncol(h_x2)
  h_dim = discrete_dim + continuous_dim

  if(continuous_dim == 0){
    h_x2 = NULL
  }else{
    h_x2 = as.matrix(h_x2)
  }

  model = function(x, y, alpha){
    if(!is.mnar) y = NULL
    inv_link(cbind(rep(1, n), y, x)%*%alpha)
  }
  w = function(x, y, alpha) 1/model(x, y, alpha)

  g = function(alpha){
    g.matrix = matrix(NA, n, h_dim)
    rw = r*w(x, y, alpha)
    if(discrete_dim > 0){
      for(l in 1:discrete_dim){g.matrix[, l] = d[, l]*(rw-1)}
    }
    if(continuous_dim > 0){
      for(l in (discrete_dim+1):(discrete_dim+continuous_dim)) g.matrix[, l] = h_x2[, l-discrete_dim]*(rw-1)
    }
    return(wt*g.matrix)
  }

  G = function(g.matrix){
    return(matrix(apply(g.matrix, 2, mean), h_dim, 1))
  }

  # W = function(g.matrix){
  #   return(solve(t(g.matrix)%*%g.matrix/n))
  # }

  t_Gamma_i = function(alpha){
    Gamma.arr = array(NA, dim = c(n, alpha_dim, h_dim))
    for(l in 1:h_dim){
      Gamma.arr[,,l] = jacobian(function(alpha) g(alpha)[, l], alpha)
    }
    return(Gamma.arr)
  }

  Gamma = function(alpha){
    return(t(apply(t_Gamma_i(alpha), c(2, 3), mean)))
  }

  Gamma_2 = function(alpha){
    return(jacobian(function(alpha) as.vector(Gamma(alpha)), alpha))
  }

  obj = function(alpha){
    g.matrix = g(alpha)
    G.hat = G(g.matrix)
    value = t(G.hat)%*%G.hat
    return(ifelse(is.infinite(value) || is.na(value), 10^8, value))
  }

  if(is.null(init)) init = rep(0, alpha_dim)
  alpha_sol_path = matrix(init, alpha_dim)
  conv_err = 10^8
  t = 1

  while (conv_err > 10^(-8) & t < 1000){
    opt = optim(alpha_sol_path[,t], obj, method = "L-BFGS-B")
    alpha_sol_path = cbind(alpha_sol_path, opt$par)
    g.matrix = g(alpha_sol_path[,t+1]); W.hat = W(g.matrix);
    obj = function(alpha){
      g.matrix = g(alpha); G.hat = G(g.matrix);
      value = t(G.hat)%*%W.hat%*%G.hat
      return(ifelse(is.infinite(value) || is.na(value), 10^8, value))
    }
    conv_err = max(abs(alpha_sol_path[,t+1]-alpha_sol_path[,t]))
    t = t + 1
  }

  alpha.hat = alpha_sol_path[, t]
  fitted_values = model(x, y, alpha.hat)
  
  Gamma.hat = Gamma(alpha.hat)
  g.matrix = g(alpha.hat)
  W.hat = W(g.matrix)
  
  Gamma.hat = W.hat = K = mu_s = psi_alpha = cov.hat = se = NA
  if(se.fit){
    mu_s = apply(g.matrix, 2, mean)
    M_n = kronecker(mu_s%*%W.hat, diag(alpha_dim))%*%Gamma_2(alpha.hat)
    H_0n = -solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)
    H_1n = t(Gamma.hat)%*%W.hat%*%(t(g.matrix)-mu_s)
    H_2n_2 = apply(t_Gamma_i(alpha.hat), 1, function(m) (m-t(Gamma.hat))%*%W.hat%*%mu_s)
    H_2n_3 = sapply(1:n, function(i) t(Gamma.hat)%*%(g.matrix[i, ]%*%t(g.matrix[i, ])-W.hat)%*%mu_s)
    psi_alpha = solve(diag(alpha_dim)-H_0n%*%M_n)%*%H_0n%*%(H_1n+H_2n_2+H_2n_3)
    # psi_alpha = solve(diag(alpha_dim)-H_0n%*%M_n)%*%H_0n%*%(H_1n+H_2n_2+H_2n_3)
    cov.hat = var(t(psi_alpha))/n
    S = var(g.matrix)
    # cov.hat = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)/N
    cov.hat = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)%*%t(Gamma.hat)%*%W.hat%*%S%*%W.hat%*%Gamma.hat%*%solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)/n
    se = sqrt(diag(cov.hat))
  }
  
  results = list(coefficients = alpha.hat,
                 fitted.values = fitted_values,
                 sol.path = alpha_sol_path,
                 cov.hat = cov.hat,
                 se = se,
                 lower = alpha.hat-qnorm(0.975)*se,
                 upper = alpha.hat+qnorm(0.975)*se,
                 g.matrix = g.matrix,
                 K = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)%*%t(Gamma.hat)%*%W.hat,
                 psi_alpha = psi_alpha,
                 model = model,
                 model_x_names = model_x_names,
                 h_x = cbind(d, h_x2))

  return(results)
}

estimate_nu_perturb = function(alpha, alpha_dim, ps_model.list, h_x, W, data, wt, init = NULL, se.fit){
  r = as.matrix(data$r)
  y = as.matrix(data$y)
  n = nrow(data)
  J = length(ps_model.list )

  h_dim = ncol(h_x)

  ps.matrix = matrix(NA, n, J)
  for(j in 1:J){
    x = as.matrix(data[ps_fit.list[[j]]$model_x_names])
    ps.matrix[, j] = ps_model.list[[j]](x, y, alpha[(1+sum(alpha_dim[0:(j-1)])):sum(alpha_dim[1:j])])
  }

  g = function(nu){
    rw = r/(ps.matrix%*%nu)
    g.matrix = as.vector(rw-1)*h_x
    return(wt*g.matrix)
  }

  G = function(g.matrix){
    return(matrix(apply(g.matrix, 2, mean), h_dim, 1))
  }

  # W = function(g.matrix){
  #   return(solve(t(g.matrix)%*%g.matrix/n))
  # }

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
                lower = rep(-Inf, J), upper = rep(Inf, J))
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

  return(nu.hat)
}

ensemble_perturb = function(ps_fit.list, h_x, W, data, wt, init = NULL, se.fit) {
  alpha.list = lapply(ps_fit.list, function(ps_fit) ps_fit$coefficients)
  alpha.hat = unlist(alpha.list)
  alpha_dim = unlist(lapply(alpha.list, length))
  ps_model.list = lapply(ps_fit.list, function(ps_fit) ps_fit$model)
  ps.matrix = do.call(cbind, lapply(ps_fit.list, function(ps_fit) ps_fit$fitted.values))

  # Basic setup
  r = as.matrix(data$r)
  n = nrow(data)
  # J = ncol(ps.matrix)

  result = separate_variable_types(h_x)
  h_x1 = result$x1
  h_x2 = result$x2

  d = NULL
  if(ncol(h_x1) > 0){
    for(j in 1:ncol(h_x1)) h_x1[, j] = as.factor(h_x1[, j])
    d = model.matrix(lm(rep(1, n)~., data =  h_x1))
  }

  if(ncol(h_x2) == 0){
    h_x2 = NULL
  }else{
    h_x2 = as.matrix(h_x2)
  }

  h_x = cbind(d, h_x2)
  h_dim = ncol(h_x)

  g = function(nu){
    rw = r/(ps.matrix%*%nu)
    g.matrix = as.vector(rw-1)*h_x
    return(wt*g.matrix)
  }

  G = function(g.matrix){
    return(matrix(apply(g.matrix, 2, mean), h_dim, 1))
  }

  # W = function(g.matrix){
  #   return(solve(t(g.matrix)%*%g.matrix/n))
  # }

  t_Gamma_i = function(nu){
    Gamma.arr = array(NA, dim = c(n, J, h_dim))
    for(l in 1:h_dim){
      Gamma.arr[,,l] = jacobian(function(nu) g(nu)[, l], nu)
    }
    return(Gamma.arr)
  }

  Gamma = function(nu){
    return(t(apply(t_Gamma_i(nu), c(2, 3), mean)))
  }

  Gamma_2 = function(nu){
    return(jacobian(function(nu) as.vector(Gamma(nu)), nu))
  }

  nu.hat = estimate_nu_perturb(alpha.hat, alpha_dim, ps_model.list, h_x, W, data, wt, init, se.fit)
  
  Gamma.hat = Gamma(nu.hat)
  g.matrix = g(nu.hat)
  W.hat = W(g.matrix)
  
  Gamma.hat = W.hat = K = Q = mu_s = psi_nu = cov.hat = se = NA
  if(se.fit){
    mu_s = apply(g.matrix, 2, mean)
    M_n = kronecker(mu_s%*%W.hat, diag(J))%*%Gamma_2(nu.hat)
    H_0n = -solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)
    H_1n = t(Gamma.hat)%*%W.hat%*%(t(g.matrix)-mu_s)
    H_2n_2 = apply(t_Gamma_i(nu.hat), 1, function(m) (m-t(Gamma.hat))%*%W.hat%*%mu_s)
    H_2n_3 = sapply(1:n, function(i) t(Gamma.hat)%*%(g.matrix[i, ]%*%t(g.matrix[i, ])-W.hat)%*%mu_s)
    psi_nu = solve(diag(J)-H_0n%*%M_n)%*%H_0n%*%(H_1n+H_2n_2+H_2n_3)
    cov.hat = var(t(psi_nu))/n
    # S = var(g.matrix)
    # # cov.hat = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)/N
    # cov.hat = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)%*%t(Gamma.hat)%*%W.hat%*%S%*%W.hat%*%Gamma.hat%*%solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)/n
    se = sqrt(diag(cov.hat))
  }

  results = list(coefficients = nu.hat,
                 # sol.path = nu_sol_path,
                 cov.hat = cov.hat,
                 se = se,
                 lower = nu.hat-qnorm(0.975)*se,
                 upper = nu.hat+qnorm(0.975)*se,
                 g.matrix = g.matrix,
                 Gamma = Gamma.hat,
                 W = W.hat,
                 # dot_nu.hat = dot_nu.hat,
                 # K = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)%*%t(Gamma.hat)%*%W.hat,
                 # Q = solve(diag(J)-H_0n%*%M_n)%*%H_0n,
                 mu_s = mu_s,
                 psi_nu = psi_nu,
                 h_x = cbind(d, h_x2))

  return(results)
}

EBMR_IPW_perturb = function(h_x_names, W, data, wt, se.fit = TRUE, true_ps = NULL) {
  # Basic setup
  r = as.matrix(data$r)
  y = as.matrix(data$y)
  n = nrow(data)

  ################################################################################
  # Collect the propensity score models
  ################################################################################
  J = length(ps_fit.list)
  ps_model.list = lapply(ps_fit.list, function(ps_fit) ps_fit$model)
  alpha.list = lapply(ps_fit.list, function(ps_fit) ps_fit$coefficients)
  # alpha = unlist(lapply(ps_fit.list, function(ps_fit) ps_fit$coefficients))
  # alpha_var.list = lapply(ps_fit.list, function(ps_fit) (ps_fit$se)^2)
  alpha_dim = unlist(lapply(alpha.list, length))
  ps.matrix = do.call(cbind, lapply(ps_fit.list, function(ps_fit) ps_fit$fitted.values))
  ################################################################################

  ################################################################################
  # Ensemble step
  ################################################################################
  ensemble_fit = ensemble_perturb(ps_fit.list, data[h_x_names], W, data, wt, init = rep(1/J, J), se.fit)
  # ensemble_fit = ensemble(ps_fit.list, data[h_x_names], W, data, init = rep(1/J, J), se.fit)
  # ensemble_fit = estimate_nu(ps.matrix, data[h_x_names], data, init = rep(1/J, J))
  nu.hat = ensemble_fit$coefficients
  # W = function(v) v^2/sum(v^2)
  # V_inv = if (J == 1) {
  #   1/sum(alpha_var.list[[1]])
  # } else {
  #   diag(sapply(alpha_var.list, function(alpha_var) 1/sum(alpha_var)))
  # }
  # w.hat = W(W(V_inv%*%nu.hat))
  w.hat = nu.hat^2/sum(nu.hat^2)
  ensemble_ps = ps.matrix%*%w.hat
  ################################################################################
  
  ################################################################################
  # IPW estimator for the population mean mu_0 with propensity score being estimated
  # by the methods of Wang, Shao and Kim (2014).
  ################################################################################
  mu_ipw = mean(wt*r/ensemble_ps*y)
  # mu_ipw.iid = as.vector(t(r/ensemble_ps*y)
  #                        +(t(H_alpha.w)-t(w.H_nu)%*%K_nu%*%t(E_dot_g))%*%K_alpha%*%g_all
  #                        +t(w.H_nu)%*%K_nu%*%g)
  # mu_ipw.iid = as.vector(t(r/ensemble_ps*y)
  #                        +(t(H_alpha.w)-2*t(w.H_nu)%*%K_nu%*%t(E_dot_g))%*%K_alpha%*%g_all
  #                        +t(w.H_nu)%*%K_nu%*%g)
  se_ipw = NA
  if(se.fit){
    ################################################################################
    # Compute necessary quantities to estimate the influence function:
    # \psi(\bm{\alpha}_*, \bm{\nu}_*)
    ################################################################################
    dot_pi = matrix(NA, n, sum(alpha_dim))
    for(j in 1:J){
      x = as.matrix(data[ps_fit.list[[j]]$model_x_names])
      dot_pi[, (sum(alpha_dim[0:(j-1)])+1):sum(alpha_dim[1:j])] = jacobian(function(alpha) ps_model.list[[j]](x, y, alpha), alpha.list[[j]])
    }
    # E_dot_g = -(t(dot_pi)*rep(nu.hat, alpha_dim))%*%(ensemble_fit$h_x*as.vector(r*((ensemble_ps)^(-2))))/n
    
    
    dot_W = function(nu){
      nu = as.vector(nu)
      (diag(2*nu)*sum(nu^2)-2*(nu)%*%t(nu^2))/(sum(nu^2)^2)
    }
    # dot_W_nu_hat = 0; if(length(nu.hat) > 1) dot_W_nu_hat = dot_W(W(V_inv%*%nu.hat))%*%dot_W(V_inv%*%nu.hat)%*%V_inv
    dot_W_nu_hat = 0; if(length(nu.hat) > 1) dot_W_nu_hat = dot_W(nu.hat)
    
    H_alpha.w = apply(t((t(dot_pi)*rep(w.hat, alpha_dim)))*as.vector(r*y*((ensemble_ps)^(-2))), 2, mean)
    w.H_nu = apply(ps.matrix%*%t(dot_W_nu_hat)*as.vector(r*y*((ensemble_ps)^(-2))), 2, mean)
    
    # K_alpha = c(lapply(ps_fit.list, function(ps_fit) ps_fit$K)) %>% bdiag()
    # g_all = do.call(rbind, lapply(ps_fit.list, function(ps_fit) t(ps_fit$g.matrix)))
    psi_alpha = do.call(rbind, lapply(ps_fit.list, function(ps_fit) ps_fit$psi_alpha))
    # K_nu = ensemble_fit$K
    # g = t(ensemble_fit$g.matrix)
    psi_nu = ensemble_fit$psi_nu
    
    Gamma_nu = ensemble_fit$Gamma
    W_nu = ensemble_fit$W
    f_n = ensemble_fit$g.matrix
    h_nu = ensemble_fit$h_x
    mu_s = ensemble_fit$mu_s
    K_1 = t(Gamma_nu)%*%W_nu%*%t((t(dot_pi)*rep(nu.hat, alpha_dim))%*%(h_nu*as.vector(-r*((ps.matrix%*%nu.hat)^(-2))))/n)
    K_2_1 = (t(dot_pi)*rep(nu.hat, alpha_dim))%*%(as.vector(2*((ps.matrix%*%nu.hat)^(-3)))*as.vector(r*h_nu%*%W_nu%*%mu_s)*ps.matrix)/n
    K_2_2 = apply(as.vector((ps.matrix%*%nu.hat)^(-2))*as.vector(r*h_nu%*%W_nu%*%mu_s)*dot_pi, 2, mean)
    K_2_2 = bdiag(lapply(1:length(alpha_dim), function(j) matrix(K_2_2[(1+sum(alpha_dim[0:(j-1)])):sum(alpha_dim[0:(j)])], ncol = 1)))
    K_2 = t(K_2_1 - as.matrix(K_2_2))
    K_3 = t(Gamma_nu)%*%t((t(dot_pi)*rep(nu.hat, alpha_dim))%*%(h_nu*as.vector(-r*((ps.matrix%*%nu.hat)^(-2)))*as.vector(f_n%*%mu_s))/n)
    ################################################################################
    
    mu_ipw.iid = as.vector(t(r/ensemble_ps*y)
                           -(t(H_alpha.w)+t(w.H_nu)%*%ensemble_fit$Q%*%(K_1+K_2+K_3))%*%psi_alpha
                           -t(w.H_nu)%*%psi_nu)
    se_ipw = sqrt(var(wt*mu_ipw.iid)/n)
  }
  ################################################################################
  
  ################################################################################
  # IPW estimator for the population mean mu_0 with known propensity score.
  ################################################################################
  mu_ipw.true = NA
  se_ipw.true = NA
  if(!is.null(true_ps)){
    mu_ipw.true = mean(r/true_ps*y)
    mu_ipw.true.iid = as.vector(r/true_ps*y)
    se_ipw.true = sqrt(var(mu_ipw.true.iid)/n)
  }
  ################################################################################

  result = list(mu_ipw = mu_ipw,
                mu_ipw.true = mu_ipw.true,
                se_ipw = se_ipw,
                se_ipw.true = se_ipw.true,
                ps.matrix = ps.matrix,
                nu.hat = nu.hat,
                w.hat = w.hat
                # imbalance = sum(apply(ensemble_fit$g, 1, mean)^2)
  )

  return(result)
}
