#' Wang, Shao, and Kim (2014) Method for Propensity Score Estimation
#'
#' Implements the Wang, Shao, and Kim (2014) approach for estimating propensity score models in the presence of missing-not-at-random (MNAR) data.
#'
#' @import numDeriv
#'
#' @param formula A formula specifying the relationship between the response and predictors.
#' @param h_x_names A character vector of variable names to be balanced.
#' @param inv_link An inverse link function applied to the linear predictor.
#' @param init (Optional) A numeric vector specifying the initial values for the model parameters. Defaults to `NULL`.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{\code{coefficients}}{Estimated model coefficients.}
#'   \item{\code{fitted.values}}{Predicted propensity scores.}
#'   \item{\code{sol.path}}{Solution path of parameter estimates across iterations.}
#'   \item{\code{cov.hat}}{Estimated covariance matrix of the coefficients.}
#'   \item{\code{se}}{Standard errors of the estimated coefficients.}
#'   \item{\code{lower}}{Lower bound of the 95% confidence intervals.}
#'   \item{\code{upper}}{Upper bound of the 95% confidence intervals.}
#'   \item{\code{g.matrix}}{Matrix of moment conditions used in the estimation.}
#'   \item{\code{K}}{Matrix related to the influence function.}
#'   \item{\code{model}}{The propensity score model function.}
#'   \item{\code{is.mnar}}{Logical indicator of whether the outcome is missing-not-at-random.}
#'   \item{\code{model_x_names}}{Names of the predictor variables used in the model.}
#'   \item{\code{h_x}}{Matrix of covariates to be balanced.}
#' }
#'
#' @examples
#' \dontrun{
#' data <- data.frame(x = rnorm(100), y = rnorm(100), r = rbinom(100, 1, 0.5))
#' ebmr <- EBMRAlgorithm$new(data = data)
#' result <- ebmr$WangShaoKim2014(
#'   formula = r ~ x,
#'   h_x_names = c("x"),
#'   inv_link = function(eta) 1 / (1 + exp(-eta))
#' )
#' print(result$coefficients)
#' }
#'
#' @references Wang, Shao, & Kim (2014). "An instrumental variable approach for identification and estimation with nonignorable nonresponse."

WangShaoKim2014 = function(formula, h_x_names, inv_link, wt = NULL, se.fit = T, init = NULL) {
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

  Phi_alpha = function(param){
    rw = r/model(x, y, param)
    g.matrix = as.vector(rw-1)*h_x
    return(g.matrix)
  }
  if(!is.null(wt)) Phi_alpha = function(param) wt*Phi_alpha(param)
  gmm_fit = gmm(Phi_alpha, W, n, h_dim, alpha_dim, init, se.fit)

  results = list(coefficients = gmm_fit$estimates,
                 se = gmm_fit$se,
                 fitted.values = model(x, y, gmm_fit$estimates),
                 model = model,
                 model_x_names = model_x_names,
                 h_x = h_x,
                 gmm_fit = gmm_fit)

  return(results)
}

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

ensemble = function(ps.matrix, h_x_names, init = NULL, se.fit = T, wt = NULL) {
  # Basic setup
  r = as.matrix(r)
  n = length(r)
  J = ncol(ps.matrix)

  result = separate_variable_types(data[h_x_names])
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

  Phi_nu = function(param){
    rw = r/(ps.matrix%*%param)
    g.matrix = as.vector(rw-1)*h_x
    return(g.matrix)
  }
  if(!is.null(wt)) Phi_nu = function(param) wt*Phi_nu(param)
  gmm_fit = gmm(Phi_nu, W, n, h_dim, J, init, se.fit)

  results = list(coefficients = gmm_fit$estimates,
                 h_x = h_x,
                 gmm_fit = gmm_fit)

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

EBMR_IPW = function(h_x_names, se.fit = TRUE, true_ps = NULL, wt = NULL) {
  # Basic setup
  r = as.matrix(r)
  y = as.matrix(y)
  n = length(y)
  J = J

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
  ensemble_fit = ensemble(ps.matrix, h_x_names, init = rep(1/J, J), se.fit)
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
      x = as.matrix(data[ps_fit.list[[j]]$model_x_names])
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
    psi_nu = ensemble_fit$gmm_fit$psi

    Gamma_nu = ensemble_fit$gmm_fit$Gamma.hat
    W_nu = ensemble_fit$gmm_fit$W.hat
    Phi_nu = ensemble_fit$gmm_fit$g.matrix
    h_nu = ensemble_fit$h_x

    # f_n = ensemble_fit$gmm_fit$g.matrix
    # eta_s = ensemble_fit$gmm_fit$eta_s
    # h_nu = ensemble_fit$h_x
    # K_1 = t(Gamma_nu)%*%W_nu%*%t((t(dot_pi)*rep(nu.hat, alpha_dim))%*%(h_nu*as.vector(-r*((ps.matrix%*%nu.hat)^(-2)))))/n
    # K_2_1 = (t(dot_pi)*rep(nu.hat, alpha_dim))%*%(as.vector(2*((ps.matrix%*%nu.hat)^(-3)))*as.vector(r*h_nu%*%W_nu%*%eta_s)*ps.matrix)/n
    # K_2_2 = apply(as.vector((ps.matrix%*%nu.hat)^(-2))*as.vector(r*h_nu%*%W_nu%*%eta_s)*dot_pi, 2, mean)
    # K_2_2 = bdiag(lapply(1:length(alpha_dim), function(j) matrix(K_2_2[(1+sum(alpha_dim[0:(j-1)])):sum(alpha_dim[0:(j)])], ncol = 1)))
    # K_2 = t(K_2_1 - as.matrix(K_2_2))
    # # K_3 = t(Gamma_nu)%*%t((t(dot_pi)*rep(nu.hat, alpha_dim))%*%(h_nu*as.vector(-r*((ps.matrix%*%nu.hat)^(-2)))*as.vector(f_n%*%eta_s))/n)
    # K_3 = -W_nu%*%t((t(dot_pi)*rep(nu.hat, alpha_dim))%*%(h_nu*as.vector(-r*((ps.matrix%*%nu.hat)^(-2)))*as.vector(f_n%*%W_nu%*%eta_s)))/n
    # K_3 = K_3-W_nu%*%t(f_n)%*%(t(t(dot_pi)*rep(nu.hat, alpha_dim))*as.vector(h_nu%*%W_nu%*%eta_s*(-r*((ps.matrix%*%nu.hat)^(-2)))))/n
    # K_3 = t(Gamma_nu)%*%K_3

    Phi_nu.alpha = ((t(h_nu*as.vector(-r*((ps.matrix%*%nu.hat)^(-2))))%*%t(t(dot_pi)*rep(nu.hat, alpha_dim)))/n)
    dot_nu = -solve(t(Gamma_nu)%*%W_nu%*%Gamma_nu)%*%t(Gamma_nu)%*%W_nu%*%Phi_nu.alpha
    #--------------------------------------------------------------------------#

    # mu_ipw.iid = as.vector(t(r/ensemble_ps*y)
    #                        -(t(H_alpha.w)+t(w.H_nu)%*%ensemble_fit$gmm_fit$Q%*%(K_1+K_2+K_3))%*%psi_alpha
    #                        -t(w.H_nu)%*%psi_nu)
    # se_ipw = sqrt(var(mu_ipw.iid)/n)

    mu_ipw.iid = as.vector(t(r/ensemble_ps*y)
                           -(t(H_alpha.w)+t(w.H_nu)%*%dot_nu)%*%psi_alpha
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

EBMR_IPW_with_locally_misspecified_model = function(ps.matrix, perturb_ps, exp_tilt, exp_tilt_x_names, h_x_names, se.fit = FALSE){
  # Basic setup
  r = as.matrix(r)
  y = as.matrix(y)
  n = length(r)

  #-----------------------------------------------------------------------------#
  # Perturb the designated propensity score model with the exponential tilt model
  #-----------------------------------------------------------------------------#
  J = ncol(ps.matrix)
  ps.matrix[, perturb_ps] = exp_tilt(y, data[, exp_tilt_x_names])*ps.matrix[, perturb_ps]
  #-----------------------------------------------------------------------------#

  #-----------------------------------------------------------------------------#
  # Ensemble step
  #-----------------------------------------------------------------------------#
  ensemble_fit = ensemble(ps.matrix, h_x_names, init = rep(1/J, J), se.fit)
  nu.hat = ensemble_fit$coefficients
  w.hat = nu.hat^2/sum(nu.hat^2)
  ensemble_ps = ps.matrix%*%w.hat
  #-----------------------------------------------------------------------------#

  #-----------------------------------------------------------------------------#
  # IPW estimator for the population mean mu_0 with propensity score being estimated
  # by the methods of Wang, Shao and Kim (2014).
  #-----------------------------------------------------------------------------#
  mu_ipw = mean(r/ensemble_ps*y)
  se_ipw = NA
  if(se.fit){
    #--------------------------------------------------------------------------#
    # Compute necessary quantities to estimate the influence function:
    # \psi(\bm{\alpha}_*, \bm{\nu}_*)
    #--------------------------------------------------------------------------#
    dot_pi = matrix(NA, n, sum(alpha_dim))
    for(j in 1:J){
      x = as.matrix(data[ps_fit.list[[j]]$model_x_names])
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
    psi_nu = ensemble_fit$gmm_fit$psi

    Gamma_nu = ensemble_fit$gmm_fit$Gamma.hat
    W_nu = ensemble_fit$gmm_fit$W.hat
    f_n = ensemble_fit$gmm_fit$g.matrix
    eta_s = ensemble_fit$gmm_fit$eta_s
    h_nu = ensemble_fit$h_x
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

  result = list(mu_ipw = mu_ipw,
                se_ipw = se_ipw,
                ps.matrix = ps.matrix,
                nu.hat = nu.hat,
                w.hat = w.hat
  )

  return(result)
}
