WangShaoKim2014 = function(formula, h_x_names, inv_link, W, wt = NULL, se.fit = T, init = NULL) {
  # Basic setup
  result = private$parse_formula(formula)
  r = as.matrix(self$data[result$r_names])
  y = as.matrix(self$data[result$y_names])
  x = as.matrix(self$data[result$x_names])
  n = nrow(self$data)
  model_x_names = colnames(x)

  result = private$separate_variable_types(self$data[h_x_names])
  h_x1 = result$x1
  h_x2 = result$x2

  is.mnar = ifelse(ncol(y) == 0, FALSE, TRUE)
  alpha_dim = 1 + as.numeric(is.mnar) + ncol(x)

  d = NULL
  if(ncol(h_x1) > 0){
    for(j in 1:ncol(h_x1)) h_x1[, j] = as.factor(h_x1[, j])
    d = model.matrix(lm(rep(1, n)~., self$data =  h_x1))
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

  Phi_alpha = private$Phi_alpha
  if(!is.null(wt)) Phi_alpha = function(param) wt*private$Phi_alpha(param)
  gmm_fit = private$gmm(private$Phi_alpha, private$W, alpha_dim, init, se.fit)

  results = list(coefficients = alpha.hat,
                 fitted.values = fitted_values,
                 model = model,
                 model_x_names = model_x_names,
                 h_x = h_x,
                 gmm_fit = gmm_fit)

  return(results)
}

ensemble = function(ps_fit.list, h_x_names, W, init = NULL, se.fit = T, wt = NULL) {
  # Basic setup
  r = as.matrix(private$r)
  n = private$n
  J = private$J

  result = private$separate_variable_types(self$data[h_x_names])
  h_x1 = result$x1
  h_x2 = result$x2

  d = NULL
  if(ncol(h_x1) > 0){
    for(j in 1:ncol(h_x1)) h_x1[, j] = as.factor(h_x1[, j])
    d = model.matrix(lm(rep(1, n)~., self$data =  h_x1))
  }

  if(ncol(h_x2) == 0){
    h_x2 = NULL
  }else{
    h_x2 = as.matrix(h_x2)
  }

  h_x = cbind(d, h_x2)
  h_dim = ncol(h_x)

  Phi_nu = private$Phi_nu
  if(!is.null(wt)) Phi_nu = function(param) wt*private$Phi_nu(param)
  gmm_fit = private$gmm(Phi_nu, private$W, J, init, se.fit)

  results = list(coefficients = gmm_fit$estimates,
                 h_x = h_x,
                 gmm_fit = gmm_fit)

  return(results)
}

EBMR_IPW = function(h_x_names, W, se.fit = TRUE, true_ps = NULL, wt = NULL) {
  # Basic setup
  r = as.matrix(private$r)
  y = as.matrix(private$y)
  n = private$n
  J = private$J
  ps_fit.list = self$ps_fit.list

  #-----------------------------------------------------------------------------#
  # Collect the propensity score models
  #-----------------------------------------------------------------------------#
  J = length(ps_fit.list)
  alpha.list = lapply(ps_fit.list, function(ps_fit) ps_fit$coefficients)
  alpha.hat = unlist(alpha.list)
  alpha_dim = unlist(lapply(alpha.list, length))
  ps_model.list = lapply(ps_fit.list, function(ps_fit) ps_fit$model)
  ps.matrix = do.call(cbind, lapply(ps_fit.list, function(ps_fit) ps_fit$fitted.values))
  #-----------------------------------------------------------------------------#

  #-----------------------------------------------------------------------------#
  # Ensemble step
  #-----------------------------------------------------------------------------#
  ensemble_fit = ensemble(ps_fit.list, h_x_names, private$W, init = rep(1/J, J), se.fit)
  nu.hat = ensemble_fit$coefficients
  w.hat = nu.hat^2/sum(nu.hat^2)
  ensemble_ps = ps.matrix%*%w.hat
  #-----------------------------------------------------------------------------#

  #-----------------------------------------------------------------------------#
  # IPW estimator for the population mean mu_0 with propensity score being estimated
  # by the methods of Wang, Shao and Kim (2014).
  #-----------------------------------------------------------------------------#
  mu_ipw = ifelse(is.null(wt), mean(r/ensemble_ps*y), mean(wt*r/ensemble_ps*y))
  se_ipw = NA
  if(se.fit){
    #--------------------------------------------------------------------------#
    # Compute necessary quantities to estimate the influence function:
    # \psi(\bm{\alpha}_*, \bm{\nu}_*)
    #--------------------------------------------------------------------------#
    dot_pi = matrix(NA, n, sum(alpha_dim))
    for(j in 1:J){
      x = as.matrix(self$data[ps_fit.list[[j]]$model_x_names])
      dot_pi[, (sum(alpha_dim[0:(j-1)])+1):sum(alpha_dim[1:j])] = jacobian(function(alpha) ps_model.list[[j]](x, y, alpha), alpha.list[[j]])
    }

    dot_W = function(nu){
      nu = as.vector(nu)
      (diag(2*nu)*sum(nu^2)-2*(nu)%*%t(nu^2))/(sum(nu^2)^2)
    }
    dot_W_nu_hat = 0; if(length(nu.hat) > 1) dot_W_nu_hat = dot_W(nu.hat)

    H_alpha.w = apply(t((t(dot_pi)*rep(w.hat, alpha_dim)))*as.vector(r*y*((ensemble_ps)^(-2))), 2, mean)
    w.H_nu = apply(ps.matrix%*%t(dot_W_nu_hat)*as.vector(r*y*((ensemble_ps)^(-2))), 2, mean)

    psi_alpha = do.call(rbind, lapply(ps_fit.list, function(ps_fit) ps_fit$gmm_fit$psi))
    psi_nu = ensemble_fit$gmm_fit$psi_nu

    Gamma_nu = ensemble_fit$gmm_fit$Gamma.hat
    W_nu = ensemble_fit$gmm_fit$W.hat
    f_n = ensemble_fit$gmm_fit$g.matrix
    h_nu = ensemble_fit$gmm_fit$h_x
    eta_s = ensemble_fit$gmm_fit$eta_s
    K_1 = t(Gamma_nu)%*%W_nu%*%t((t(dot_pi)*rep(nu.hat, alpha_dim))%*%(h_nu*as.vector(-r*((ps.matrix%*%nu.hat)^(-2))))/n)
    K_2_1 = (t(dot_pi)*rep(nu.hat, alpha_dim))%*%(as.vector(2*((ps.matrix%*%nu.hat)^(-3)))*as.vector(r*h_nu%*%W_nu%*%eta_s)*ps.matrix)/n
    K_2_2 = apply(as.vector((ps.matrix%*%nu.hat)^(-2))*as.vector(r*h_nu%*%W_nu%*%eta_s)*dot_pi, 2, mean)
    K_2_2 = bdiag(lapply(1:length(alpha_dim), function(j) matrix(K_2_2[(1+sum(alpha_dim[0:(j-1)])):sum(alpha_dim[0:(j)])], ncol = 1)))
    K_2 = t(K_2_1 - as.matrix(K_2_2))
    K_3 = t(Gamma_nu)%*%t((t(dot_pi)*rep(nu.hat, alpha_dim))%*%(h_nu*as.vector(-r*((ps.matrix%*%nu.hat)^(-2)))*as.vector(f_n%*%eta_s))/n)
    #--------------------------------------------------------------------------#

    mu_ipw.iid = as.vector(t(r/ensemble_ps*y)
                           -(t(H_alpha.w)+t(w.H_nu)%*%ensemble_fit$gmm_fit$Q%*%(K_1+K_2+K_3))%*%psi_alpha
                           -t(w.H_nu)%*%psi_nu)
    se_ipw = sqrt(var(mu_ipw.iid)/n)
  }
  #-----------------------------------------------------------------------------#

  #-----------------------------------------------------------------------------#
  # IPW estimator for the population mean mu_0 with known propensity score.
  #-----------------------------------------------------------------------------#
  mu_ipw.true = NA
  se_ipw.true = NA
  if(!is.null(true_ps)){
    mu_ipw.true = mean(r/true_ps*y)
    mu_ipw.true.iid = as.vector(r/true_ps*y)
    se_ipw.true = sqrt(var(mu_ipw.true.iid)/n)
  }
  #-----------------------------------------------------------------------------#

  result = list(mu_ipw = mu_ipw,
                mu_ipw.true = mu_ipw.true,
                se_ipw = se_ipw,
                se_ipw.true = se_ipw.true,
                ps.matrix = ps.matrix,
                nu.hat = nu.hat,
                w.hat = w.hat
  )

  return(result)
}

EBMR_IPW_with_locally_misspecified_model = function(ps.matrix, perturb_ps, exp_tilt, exp_tilt_x_names, h_x_names){
  # Basic setup
  r = as.matrix(private$r)
  y = as.matrix(private$y)
  n = private$n

  ################################################################################
  # Perturb the designated propensity score model with the exponential tilt model
  ################################################################################
  J = ncol(ps.matrix)
  ps.matrix[, perturb_ps] = exp_tilt(y, self$data[, exp_tilt_x_names])*ps.matrix[, perturb_ps]
  ################################################################################

  #-----------------------------------------------------------------------------#
  # Ensemble step
  #-----------------------------------------------------------------------------#
  ensemble_fit = ensemble(ps_fit.list, h_x_names, private$W, init = rep(1/J, J), se.fit)
  nu.hat = ensemble_fit$coefficients
  w.hat = nu.hat^2/sum(nu.hat^2)
  ensemble_ps = ps.matrix%*%w.hat
  #-----------------------------------------------------------------------------#

  #-----------------------------------------------------------------------------#
  # IPW estimator for the population mean mu_0 with propensity score being estimated
  # by the methods of Wang, Shao and Kim (2014).
  #-----------------------------------------------------------------------------#
  mu_ipw = ifelse(is.null(wt), mean(r/ensemble_ps*y), mean(wt*r/ensemble_ps*y))
  se_ipw = NA
  if(se.fit){
    #--------------------------------------------------------------------------#
    # Compute necessary quantities to estimate the influence function:
    # \psi(\bm{\alpha}_*, \bm{\nu}_*)
    #--------------------------------------------------------------------------#
    dot_pi = matrix(NA, n, sum(alpha_dim))
    for(j in 1:J){
      x = as.matrix(self$data[ps_fit.list[[j]]$model_x_names])
      dot_pi[, (sum(alpha_dim[0:(j-1)])+1):sum(alpha_dim[1:j])] = jacobian(function(alpha) ps_model.list[[j]](x, y, alpha), alpha.list[[j]])
    }

    dot_W = function(nu){
      nu = as.vector(nu)
      (diag(2*nu)*sum(nu^2)-2*(nu)%*%t(nu^2))/(sum(nu^2)^2)
    }
    dot_W_nu_hat = 0; if(length(nu.hat) > 1) dot_W_nu_hat = dot_W(nu.hat)

    H_alpha.w = apply(t((t(dot_pi)*rep(w.hat, alpha_dim)))*as.vector(r*y*((ensemble_ps)^(-2))), 2, mean)
    w.H_nu = apply(ps.matrix%*%t(dot_W_nu_hat)*as.vector(r*y*((ensemble_ps)^(-2))), 2, mean)

    psi_alpha = do.call(rbind, lapply(ps_fit.list, function(ps_fit) ps_fit$gmm_fit$psi))
    psi_nu = ensemble_fit$gmm_fit$psi_nu

    Gamma_nu = ensemble_fit$gmm_fit$Gamma.hat
    W_nu = ensemble_fit$gmm_fit$W.hat
    f_n = ensemble_fit$gmm_fit$g.matrix
    h_nu = ensemble_fit$gmm_fit$h_x
    eta_s = ensemble_fit$gmm_fit$eta_s
    K_1 = t(Gamma_nu)%*%W_nu%*%t((t(dot_pi)*rep(nu.hat, alpha_dim))%*%(h_nu*as.vector(-r*((ps.matrix%*%nu.hat)^(-2))))/n)
    K_2_1 = (t(dot_pi)*rep(nu.hat, alpha_dim))%*%(as.vector(2*((ps.matrix%*%nu.hat)^(-3)))*as.vector(r*h_nu%*%W_nu%*%eta_s)*ps.matrix)/n
    K_2_2 = apply(as.vector((ps.matrix%*%nu.hat)^(-2))*as.vector(r*h_nu%*%W_nu%*%eta_s)*dot_pi, 2, mean)
    K_2_2 = bdiag(lapply(1:length(alpha_dim), function(j) matrix(K_2_2[(1+sum(alpha_dim[0:(j-1)])):sum(alpha_dim[0:(j)])], ncol = 1)))
    K_2 = t(K_2_1 - as.matrix(K_2_2))
    K_3 = t(Gamma_nu)%*%t((t(dot_pi)*rep(nu.hat, alpha_dim))%*%(h_nu*as.vector(-r*((ps.matrix%*%nu.hat)^(-2)))*as.vector(f_n%*%eta_s))/n)
    #--------------------------------------------------------------------------#

    mu_ipw.iid = as.vector(t(r/ensemble_ps*y)
                           -(t(H_alpha.w)+t(w.H_nu)%*%ensemble_fit$gmm_fit$Q%*%(K_1+K_2+K_3))%*%psi_alpha
                           -t(w.H_nu)%*%psi_nu)
    se_ipw = sqrt(var(mu_ipw.iid)/n)
  }
  #-----------------------------------------------------------------------------#

  #-----------------------------------------------------------------------------#
  # IPW estimator for the population mean mu_0 with known propensity score.
  #-----------------------------------------------------------------------------#
  mu_ipw.true = NA
  se_ipw.true = NA
  if(!is.null(true_ps)){
    mu_ipw.true = mean(r/true_ps*y)
    mu_ipw.true.iid = as.vector(r/true_ps*y)
    se_ipw.true = sqrt(var(mu_ipw.true.iid)/n)
  }
  #-----------------------------------------------------------------------------#

  result = list(mu_ipw = mu_ipw,
                mu_ipw.true = mu_ipw.true,
                se_ipw = se_ipw,
                se_ipw.true = se_ipw.true,
                ps.matrix = ps.matrix,
                nu.hat = nu.hat,
                w.hat = w.hat
  )

  return(result)
}
