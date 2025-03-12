#' Estimate the Coefficients for the Propensity Score Model
#'
#' This method estimates the coefficients \eqn{\nu} using an Generalized Method of Moments (GMM) procedure.
#'
#' @import Matrix
#' @import numDeriv
#'
#' @param ps.matrix A matrix of propensity scores for each observation and model.
#' @param h_x A matrix of covariates, which includes both continuous and discrete variables.
#' @param init (optional) A vector of initial values for the optimization. Default is \eqn{\bf{0}}.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{\code{coefficients}}{The estimated coefficients \eqn{\nu}.}
#'   \item{\code{sol.path}}{Solution path for \eqn{\nu} across iterations.}
#'   \item{\code{cov.hat}}{Estimated covariance matrix for \eqn{\nu}.}
#'   \item{\code{se}}{Standard errors of the estimated coefficients.}
#'   \item{\code{lower}}{Lower bound of the 95% confidence interval for \eqn{\nu}.}
#'   \item{\code{upper}}{Upper bound of the 95% confidence interval for \eqn{\nu}.}
#'   \item{\code{g.matrix}}{Matrix of moment conditions used in the estimation.}
#'   \item{\code{K}}{Matrix related to the influence function.}
#'   \item{\code{h_x}}{Matrix of covariates to be balanced.}
#' }
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' ebmr <- EBMRAlgorithm$new(data = data)
#' nu_estimates <- ebmr$estimate_nu(ps.matrix = ps_data, h_x = covariates)
#' print(nu_estimates)
#' }

estimate_nu = function(ps.matrix, h_x, init = NULL) {
  # Basic setup
  r = as.matrix(private$r)
  n = private$n
  J = ncol(ps.matrix)

  result = private$separate_variable_types(h_x)
  h_x1 = result$x1
  h_x2 = result$x2

  d = NULL
  if(ncol(h_x1) > 0){
    for(j in 1:ncol(h_x1)) h_x1[, j] = as.factor(h_x1[, j])
    d = model.matrix(lm(rep(1, n)~., data =  h_x1))
  }

  discrete_dim = ncol(d)
  continuous_dim = ncol(h_x2)
  h_dim = discrete_dim + continuous_dim

  # r = as.matrix(r)
  if(continuous_dim == 0){
    h_x2 = NULL
  }else{
    h_x2 = as.matrix(h_x2)
  }
  # h_x1 = as.matrix(h_x1)
  # h_x2 = as.matrix(h_x2)

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

estimate_nu_perturb = function(ps.matrix, h_x, init = NULL) {
  # Basic setup
  r = as.matrix(private$r)
  n = private$n
  # wt = rexp(n)
  J = ncol(ps.matrix)

  result = private$separate_variable_types(h_x)
  h_x1 = result$x1
  h_x2 = result$x2

  d = NULL
  if(ncol(h_x1) > 0){
    for(j in 1:ncol(h_x1)) h_x1[, j] = as.factor(h_x1[, j])
    d = model.matrix(lm(rep(1, n)~., data =  h_x1))
  }

  discrete_dim = ncol(d)
  continuous_dim = ncol(h_x2)
  h_dim = discrete_dim + continuous_dim

  # r = as.matrix(r)
  if(continuous_dim == 0){
    h_x2 = NULL
  }else{
    h_x2 = as.matrix(h_x2)
  }
  # h_x1 = as.matrix(h_x1)
  # h_x2 = as.matrix(h_x2)

  g = function(nu){
    g.matrix = matrix(NA, n, h_dim)
    rw = r/(ps.matrix%*%nu)
    if(discrete_dim > 0){
      for(l in 1:discrete_dim){g.matrix[, l] = d[, l]*(rw-1)}
    }
    if(continuous_dim > 0){
      for(l in (discrete_dim+1):(discrete_dim+continuous_dim)) g.matrix[, l] = h_x2[, l-discrete_dim]*(rw-1)
    }
    return(rexp(n)*g.matrix)
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

#' Compute the proposed Inverse Probability Weighting (IPW) Estimator
#'
#' This method computes the Inverse Probability Weighting (IPW) estimator for
#' the population mean \eqn{\mu_0}, using the propensity scores estimated by
#' the ensemble method. It also computes standard errors and other related
#' quantities for the estimator, including the estimator when the true propensity
#' score is provided.
#'
#' @param h_x_names A character vector of variable names to be balanced.
#' @param true_ps (optional) A vector of true propensity scores. If provided, the IPW estimator
#'        will also be computed using the true propensity scores.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{\code{mu_ipw}}{The IPW estimator for the population mean \eqn{\mu_0} using estimated propensity scores.}
#'   \item{\code{mu_ipw.true}}{The IPW estimator using the true propensity scores (if provided).}
#'   \item{\code{se_ipw}}{The standard error of the IPW estimator using estimated propensity scores.}
#'   \item{\code{se_ipw.true}}{The standard error of the IPW estimator using true propensity scores (if provided).}
#'   \item{\code{ps.matrix}}{The matrix of estimated propensity scores for each observation and model.}
#'   \item{\code{nu.hat}}{The estimated coefficients \eqn{\nu} from the ensemble method.}
#'   \item{\code{w.hat}}{The estimated weights for the propensity score models.}
#'   \item{\code{imbalance}}{A measure of imbalance provided by the ensemble propensity scores.}
#' }
#'
#' @examples
#' \dontrun{
#' ebmr <- EBMRAlgorithm$new(data = data)
#' ipw_estimates <- ebmr$EBMR_IPW(h_x_names = covariates, true_ps = true_ps_data)
#' print(ipw_estimates)
#' }

EBMR_IPW = function(h_x_names, true_ps = NULL) {
  # Basic setup
  r = as.matrix(private$r)
  y = as.matrix(private$y)
  n = private$n

  ################################################################################
  # Collect the propensity score models
  ################################################################################
  J = length(self$ps_fit.list)
  ps_model.list = lapply(self$ps_fit.list, function(ps_fit) ps_fit$model)
  alpha.list = lapply(self$ps_fit.list, function(ps_fit) ps_fit$coefficients)
  alpha_dim = unlist(lapply(alpha.list, length))
  ps.matrix = do.call(cbind, lapply(self$ps_fit.list, function(ps_fit) ps_fit$fitted.values))
  ################################################################################

  ################################################################################
  # Ensemble step
  ################################################################################
  ensemble_fit = private$estimate_nu(ps.matrix, self$data[h_x_names], init = rep(1/J, J))
  nu.hat = ensemble_fit$coefficients
  w.hat = nu.hat^2/sum(nu.hat^2)
  ensemble_ps = ps.matrix%*%w.hat
  ################################################################################

  ################################################################################
  # Compute necessary quantities to estimate the influence function:
  # \psi(\bm{\alpha}_*, \bm{\nu}_*)
  ################################################################################
  dot_pi = matrix(NA, n, sum(alpha_dim))
  for(j in 1:J){
    x = as.matrix(self$data[self$ps_fit.list[[j]]$model_x_names])
    dot_pi[, (sum(alpha_dim[0:(j-1)])+1):sum(alpha_dim[1:j])] = jacobian(function(alpha) ps_model.list[[j]](x, y, alpha), alpha.list[[j]])
  }
  E_dot_g = -(t(dot_pi)*rep(w.hat, alpha_dim))%*%(ensemble_fit$h_x*as.vector(r*((ensemble_ps)^(-2))))/n

  dot_W = function(nu){
    (diag(2*nu)*sum(nu^2)-2*(nu)%*%t(nu^2))/(sum(nu^2)^2)
  }
  dot_W_nu_hat = 0; if(length(nu.hat) > 1) dot_W_nu_hat = dot_W(nu.hat)

  H_alpha.w = apply(t((t(dot_pi)*rep(w.hat, alpha_dim)))*as.vector(r*y*((ensemble_ps)^(-2))), 2, mean)
  w.H_nu = apply(ps.matrix%*%t(dot_W_nu_hat)*as.vector(r*y*((ensemble_ps)^(-2))), 2, mean)

  K_alpha = c(lapply(self$ps_fit.list, function(ps_fit) ps_fit$K)) %>% bdiag()
  g_all = do.call(rbind, lapply(self$ps_fit.list, function(ps_fit) t(ps_fit$g.matrix)))
  K_nu = ensemble_fit$K
  g = t(ensemble_fit$g.matrix)
  ################################################################################

  ################################################################################
  # IPW estimator for the population mean mu_0 with propensity score being estimated
  # by the methods of Wang, Shao and Kim (2014).
  ################################################################################
  mu_ipw = mean(r/ensemble_ps*y)
  # mu_ipw.iid = as.vector(t(r/ensemble_ps*y)
  #                        +(t(H_alpha.w)-t(w.H_nu)%*%K_nu%*%t(E_dot_g))%*%K_alpha%*%g_all
  #                        +t(w.H_nu)%*%K_nu%*%g)
  mu_ipw.iid = as.vector(t(r/ensemble_ps*y)
                         +(t(H_alpha.w)-2*t(w.H_nu)%*%K_nu%*%t(E_dot_g))%*%K_alpha%*%g_all
                         +t(w.H_nu)%*%K_nu%*%g)
  se_ipw = sqrt(var(mu_ipw.iid)/n)
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

EBMR_IPW_perturb = function(h_x_names, true_ps = NULL) {
  # Basic setup
  r = as.matrix(private$r)
  y = as.matrix(private$y)
  n = private$n
  # wt = private$wt

  ################################################################################
  # Collect the propensity score models
  ################################################################################
  J = length(self$ps_fit.list)
  ps_model.list = lapply(self$ps_fit.list, function(ps_fit) ps_fit$model)
  alpha.list = lapply(self$ps_fit.list, function(ps_fit) ps_fit$coefficients)
  alpha_dim = unlist(lapply(alpha.list, length))
  ps.matrix = do.call(cbind, lapply(self$ps_fit.list, function(ps_fit) ps_fit$fitted.values))
  ################################################################################

  ################################################################################
  # Ensemble step
  ################################################################################
  ensemble_fit = private$estimate_nu_perturb(ps.matrix, self$data[h_x_names], init = rep(1/J, J))
  nu.hat = ensemble_fit$coefficients
  w.hat = nu.hat^2/sum(nu.hat^2)
  ensemble_ps = ps.matrix%*%w.hat
  ################################################################################

  ################################################################################
  # Compute necessary quantities to estimate the influence function:
  # \psi(\bm{\alpha}_*, \bm{\nu}_*)
  ################################################################################
  dot_pi = matrix(NA, n, sum(alpha_dim))
  for(j in 1:J){
    x = as.matrix(self$data[self$ps_fit.list[[j]]$model_x_names])
    dot_pi[, (sum(alpha_dim[0:(j-1)])+1):sum(alpha_dim[1:j])] = jacobian(function(alpha) ps_model.list[[j]](x, y, alpha), alpha.list[[j]])
  }
  E_dot_g = -(t(dot_pi)*rep(w.hat, alpha_dim))%*%(ensemble_fit$h_x*as.vector(r*((ensemble_ps)^(-2))))/n

  dot_W = function(nu){
    (diag(2*nu)*sum(nu^2)-2*(nu)%*%t(nu^2))/(sum(nu^2)^2)
  }
  dot_W_nu_hat = 0; if(length(nu.hat) > 1) dot_W_nu_hat = dot_W(nu.hat)

  H_alpha.w = apply(t((t(dot_pi)*rep(w.hat, alpha_dim)))*as.vector(r*y*((ensemble_ps)^(-2))), 2, mean)
  w.H_nu = apply(ps.matrix%*%t(dot_W_nu_hat)*as.vector(r*y*((ensemble_ps)^(-2))), 2, mean)

  K_alpha = c(lapply(self$ps_fit.list, function(ps_fit) ps_fit$K)) %>% bdiag()
  g_all = do.call(rbind, lapply(self$ps_fit.list, function(ps_fit) t(ps_fit$g.matrix)))
  K_nu = ensemble_fit$K
  g = t(ensemble_fit$g.matrix)
  ################################################################################

  ################################################################################
  # IPW estimator for the population mean mu_0 with propensity score being estimated
  # by the methods of Wang, Shao and Kim (2014).
  ################################################################################
  mu_ipw = mean(rexp(n)*r/ensemble_ps*y)
  # mu_ipw.iid = as.vector(t(r/ensemble_ps*y)
  #                        +(t(H_alpha.w)-t(w.H_nu)%*%K_nu%*%t(E_dot_g))%*%K_alpha%*%g_all
  #                        +t(w.H_nu)%*%K_nu%*%g)
  mu_ipw.iid = as.vector(t(r/ensemble_ps*y)
                         +(t(H_alpha.w)-2*t(w.H_nu)%*%K_nu%*%t(E_dot_g))%*%K_alpha%*%g_all
                         +t(w.H_nu)%*%K_nu%*%g)
  # mu_ipw = mean(mu_ipw.iid)
  se_ipw = sqrt(var(mu_ipw.iid)/n)
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

#' Compute the IPW Estimator with Locally Misspecified Model
#'
#' This function computes the Inverse Probability Weighting (IPW) estimator for
#' the population mean \eqn{\mu_0}, incorporating a sensitivity analysis by
#' perturbing a specified propensity score model using an exponential tilt model.
#' The method allows users to assess the impact of local model misspecification
#' on the estimator.
#'
#' @param ps.matrix A matrix of estimated propensity scores for each observation and model.
#' @param perturb_ps An integer specifying the index of the propensity score model to be perturbed.
#'        This model is typically assumed to be the well-specified model among the candidate models.
#' @param exp_tilt A function defining the exponential tilt model.
#' @param exp_tilt_x_names A character vector specifying the covariates used in the exponential tilt model.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{\code{mu_ipw}}{The IPW estimator for the population mean \eqn{\mu_0} using the perturbed propensity scores.}
#'   \item{\code{se_ipw}}{The standard error of the IPW estimator.}
#'   \item{\code{ps.matrix}}{The updated matrix of estimated propensity scores after perturbation.}
#'   \item{\code{nu.hat}}{The estimated coefficients \eqn{\nu} from the ensemble method.}
#'   \item{\code{w.hat}}{The estimated weights for the propensity score models.}
#'   \item{\code{imbalance}}{A measure of imbalance provided by the ensemble propensity scores.}
#' }
#'
#' @details
#' The function perturbs the designated propensity score model using an exponential tilt model,
#' which is a sensitivity analysis technique for assessing local model misspecification. The
#' impact of perturbation on the IPW estimator is then examined. The ensemble method is used
#' to estimate the optimal weights for the propensity score models, ensuring balance across covariates.
#'
#' @examples
#' \dontrun{
#' ebmr <- EBMRAlgorithm$new(data = data)
#' ipw_sensitivity <- ebmr$EBMR_IPW_with_locally_misspecified_model(
#'   ps.matrix = estimated_ps_matrix,
#'   perturb_ps = 2,
#'   exp_tilt = function(y, data) exp(y * data$covariate),
#'   exp_tilt_x_names = c("covariate")
#' )
#' print(ipw_sensitivity)
#' }

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

  ################################################################################
  # Ensemble step
  ################################################################################
  ensemble_fit = private$estimate_nu(ps.matrix, self$data[h_x_names], init = rep(1/J, J))
  nu.hat = ensemble_fit$coefficients
  w.hat = nu.hat^2/sum(nu.hat^2)
  ensemble_ps = ps.matrix%*%w.hat
  ################################################################################

  ################################################################################
  # Compute necessary quantities to estimate the influence function:
  # \psi(\bm{\alpha}_*, \bm{\nu}_*)
  ###############################################################################
  dot_W = function(nu){
    (diag(2*nu)*sum(nu^2)-2*(nu)%*%t(nu^2))/(sum(nu^2)^2)
  }
  dot_W_nu_hat = 0; if(length(nu.hat) > 1) dot_W_nu_hat = dot_W(nu.hat)

  w.H_nu = apply(ps.matrix%*%t(dot_W_nu_hat)*as.vector(r*y*((ensemble_ps)^(-2))), 2, mean)

  K_nu = ensemble_fit$K
  g = t(ensemble_fit$g.matrix)
  ################################################################################

  ################################################################################
  # IPW estimator for the population mean mu_0 with propensity score being estimated
  # by the methods of Wang, Shao and Kim (2014).
  ################################################################################
  mu_ipw = mean(r/ensemble_ps*y)
  mu_ipw.iid = as.vector(t(r/ensemble_ps*y)+t(w.H_nu)%*%K_nu%*%g)
  se_ipw = sqrt(var(mu_ipw.iid)/n)
  ################################################################################

  result = list(mu_ipw = mu_ipw,
                se_ipw = se_ipw,
                ps.matrix = ps.matrix,
                nu.hat = nu.hat,
                w.hat = w.hat,
                imbalance = sum(apply(ensemble_fit$g, 1, mean)^2)
  )

  return(result)
}

# estimate_nu = function(ps.matrix, h_x, init = NULL){
#   result = separate_variable_types(h_x)
#   h_x1 = result$x1
#   h_x2 = result$x2
#
#   J = ncol(ps.matrix)
#
#   d = NULL
#   if(ncol(h_x1) > 0){
#     for(j in 1:ncol(h_x1)) h_x1[, j] = as.factor(h_x1[, j])
#     d = model.matrix(lm(rep(1, n)~., data =  h_x1))
#   }
#
#   discrete_dim = ncol(d)
#   continuous_dim = ncol(h_x2)
#   h_dim = discrete_dim + continuous_dim
#
#   r = as.matrix(r)
#   h_x1 = as.matrix(h_x1)
#   h_x2 = as.matrix(h_x2)
#
#   g = function(nu){
#     g.matrix = matrix(NA, n, h_dim)
#     rw = r/(ps.matrix%*%nu)
#     if(discrete_dim > 0){
#       for(l in 1:discrete_dim){g.matrix[, l] = d[, l]*(rw-1)}
#     }
#     if(continuous_dim > 0){
#       for(l in (discrete_dim+1):(discrete_dim+continuous_dim)) g.matrix[, l] = h_x2[, l-discrete_dim]*(rw-1)
#     }
#     return(g.matrix)
#   }
#
#   G = function(g.matrix){
#     return(matrix(apply(g.matrix, 2, mean), h_dim, 1))
#   }
#
#   W = function(g.matrix){
#     return(solve(t(g.matrix)%*%g.matrix/n))
#   }
#
#   Gamma = function(nu){
#     Gamma.arr = array(NA, dim = c(h_dim, n, J))
#     for(l in 1:h_dim){
#       Gamma.arr[l,,] = jacobian(function(nu) g(nu)[, l], nu)
#     }
#     return(apply(Gamma.arr, c(1, 3), mean))
#   }
#
#   obj = function(nu){
#     g.matrix = g(nu)
#     G.hat = G(g.matrix)
#     value = t(G.hat)%*%G.hat
#     return(ifelse(is.infinite(value) || is.na(value), 10^8, value))
#   }
#
#   if(is.null(init)) init = rep(0, J)
#   nu_sol_path = matrix(init, J)
#   conv_err = 10^8
#   t = 1
#
#   while (conv_err > 10^(-8) & t < 1000){
#     opt = optim(nu_sol_path[, t], obj, method = "L-BFGS-B",
#                 lower = rep(-Inf, ncol(ps.matrix)), upper = rep(Inf, ncol(ps.matrix)))
#     nu_sol_path = cbind(nu_sol_path, opt$par)
#     g.matrix = g(nu_sol_path[,t+1]); W.hat = W(g.matrix);
#     obj = function(nu){
#       g.matrix = g(nu)
#       G.hat = G(g.matrix)
#       value = t(G.hat)%*%W.hat%*%G.hat
#       return(ifelse(is.infinite(value) || is.na(value), 10^8, value))
#     }
#     conv_err = max(abs(nu_sol_path[,t+1]-nu_sol_path[,t]))
#     t = t + 1
#   }
#
#   nu.hat = nu_sol_path[, t]
#   Gamma.hat = Gamma(nu.hat)
#   g.matrix = g(nu.hat)
#   W.hat = W(g.matrix)
#   S = var(g.matrix)
#   # cov.hat = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)/N
#   cov.hat = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)%*%t(Gamma.hat)%*%W.hat%*%S%*%W.hat%*%Gamma.hat%*%solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)/n
#   se = sqrt(diag(cov.hat))
#
#   results = list(coefficients = nu.hat,
#                  sol.path = nu_sol_path,
#                  cov.hat = cov.hat,
#                  se = se,
#                  lower = nu.hat-qnorm(0.975)*se,
#                  upper = nu.hat+qnorm(0.975)*se,
#                  g.matrix = g.matrix,
#                  K = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)%*%t(Gamma.hat)%*%W.hat,
#                  h_x = cbind(d, h_x2))
#
#   return(results)
# }

# EBMR_IPW = function(self$ps_fit.list, h_x, true_ps = NULL){
#   ################################################################################
#   # Collect the propensity score models
#   ################################################################################
#   J = length(self$ps_fit.list)
#   ps_model.list = lapply(self$ps_fit.list, function(ps_fit) ps_fit$model)
#   alpha.list = lapply(self$ps_fit.list, function(ps_fit) ps_fit$coefficients)
#   alpha_dim = unlist(lapply(alpha.list, length))
#   ps.matrix = do.call(cbind, lapply(self$ps_fit.list, function(ps_fit) ps_fit$fitted.values))
#   ################################################################################
#
#   ################################################################################
#   # Ensemble step
#   ################################################################################
#   ensemble_fit = estimate_nu(ps.matrix, h_x, init = rep(1/J, J))
#   nu.hat = ensemble_fit$coefficients
#   w.hat = nu.hat^2/sum(nu.hat^2)
#   ensemble_ps = ps.matrix%*%w.hat
#   ################################################################################
#
#   ################################################################################
#   # Compute necessary quantities to estimate the influence function:
#   # \psi(\bm{\alpha}_*, \bm{\nu}_*)
#   ################################################################################
#   dot_pi = matrix(NA, n, sum(alpha_dim))
#   for(j in 1:J){
#     x = as.matrix(data[self$ps_fit.list[[j]]$model_x_names])
#     dot_pi[, (sum(alpha_dim[0:(j-1)])+1):sum(alpha_dim[1:j])] = jacobian(function(alpha) ps_model.list[[j]](x, y, alpha), alpha.list[[j]])
#   }
#   E_dot_g = -(t(dot_pi)*rep(w.hat, alpha_dim))%*%(ensemble_fit$h_x*as.vector(r*((ensemble_ps)^(-2))))/n
#
#   dot_W = function(w){
#     (diag(2*w)*sum(w^2)-2*(w)%*%t(w^2))/(sum(w^2)^2)
#   }
#   dot_W_nu_hat = 0; if(length(nu.hat) > 1) dot_W_nu_hat = dot_W(nu.hat)
#
#   H_alpha.w = apply(t((t(dot_pi)*rep(w.hat, alpha_dim)))*as.vector(r*y*((ensemble_ps)^(-2))), 2, mean)
#   w.H_nu = apply(ps.matrix%*%t(dot_W_nu_hat)*as.vector(r*y*((ensemble_ps)^(-2))), 2, mean)
#
#   K_alpha = c(lapply(self$ps_fit.list, function(ps_fit) ps_fit$K)) %>% bdiag()
#   g_all = do.call(rbind, lapply(self$ps_fit.list, function(ps_fit) t(ps_fit$g.matrix)))
#   K_nu = ensemble_fit$K
#   g = t(ensemble_fit$g.matrix)
#   ################################################################################
#
#   ################################################################################
#   # IPW estimator for the population mean mu_0 with propensity score being estimated
#   # by the methods of Wang, Shao and Kim (2014).
#   ################################################################################
#   mu_ipw = mean(r/ensemble_ps*y)
#   mu_ipw.iid = as.vector(t(r/ensemble_ps*y)
#                          +(t(H_alpha.w)+t(w.H_nu)%*%K_nu%*%t(E_dot_g))%*%K_alpha%*%g_all
#                          +t(w.H_nu)%*%K_nu%*%g)
#   se_ipw = sqrt(var(mu_ipw.iid)/n)
#   ################################################################################
#
#   ################################################################################
#   # IPW estimator for the population mean mu_0 with known propensity score.
#   ################################################################################
#   mu_ipw.true = NA
#   if(!is.null(true_ps)){
#     mu_ipw.true = mean(r/true_ps*y)
#     mu_ipw.true.iid = as.vector(r/true_ps*y)
#     se_ipw.true = sqrt(var(mu_ipw.true.iid)/N)
#   }
#   ################################################################################
#
#   results = list(mu_ipw = mu_ipw,
#                  mu_ipw.true = mu_ipw.true,
#                  se_ipw = se_ipw,
#                  se_ipw.true = se_ipw.true,
#                  ps.matrix = ps.matrix,
#                  nu.hat = ensemble_fit$nu.hat,
#                  w.hat = w.hat,
#                  imbalance = sum(apply(ensemble_fit$g, 1, mean)^2)
#   )
#
#   return(results)
# }

# EBMR_IPW_with_locally_misspecified_model = function(ps.matrix, perturb_ps, exp_tilt, exp_tilt_x_names, h){
#   ################################################################################
#   # Collect the propensity score models
#   ################################################################################
#   J = ncol(ps.matrix)
#   ps.matrix[, perturb_ps] = exp_tilt(y[r == 1], data[r == 1, exp_tilt_x_names])*ps.matrix[, perturb_ps]
#   imputed.ps.matrix = matrix(1, N, J)
#   imputed.ps.matrix[r == 1, ] = ps.matrix
#   ################################################################################
#
#   ################################################################################
#   # Compress the propensity score models
#   ################################################################################
#   ensemble_fit = estimate_nu(h, ps.matrix, r, ortho = FALSE, init = rep(1/J, J))
#   nu = ensemble_fit$nu.hat
#   nu = as.matrix((nu^2)/as.numeric(t(nu)%*%nu))
#   compressed.pi = rep(1, N)
#   compressed.pi[r == 1] = ps.matrix%*%nu
#   ################################################################################
#
#   ################################################################################
#   # calculate dot_W_nu_hat, nu.iid ... (for estimating the ASE)
#   ################################################################################
#   dot_W = function(w){
#     (diag(2*w)*sum(w^2)-2*(w^2)%*%t(w))/(sum(w^2)^2)
#   }
#   dot_W_nu_hat = 0; if(length(nu) > 1) dot_W_nu_hat = dot_W(ensemble_fit$nu.hat)
#
#   nu.iid = -ensemble_fit$K_nu%*%ensemble_fit$g
#   ################################################################################
#
#   ################################################################################
#   # IPW estimator for the population mean mu
#   ################################################################################
#   mu_ipw = mean(r/compressed.pi*y)
#   ################################################################################
#
#   ################################################################################
#   # Estimate the asymptotic variance of IPW
#   ################################################################################
#   ipw.beta = apply(-r*y/(compressed.pi^2)*(imputed.ps.matrix%*%dot_W_nu_hat), 2, mean)
#
#   mu_ipw.iid = as.vector(r/compressed.pi*y+ipw.beta%*%nu.iid)
#   se_ipw = sqrt(var(mu_ipw.iid)/N)
#
#   mu_ipw.true.iid = as.vector(r/compressed.pi*y)
#   se_ipw.true = sqrt(var(mu_ipw.true.iid)/N)
#   ################################################################################
#
#   return(list(mu_ipw = mu_ipw,
#               se_ipw = se_ipw,
#               ps.matrix = ps.matrix,
#               nu = ensemble_fit$nu.hat,
#               nu = nu))
# }

