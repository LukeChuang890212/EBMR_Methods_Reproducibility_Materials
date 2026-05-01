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

WangShaoKim2014 = function(formula, h_x_names, inv_link, init = NULL, se.fit = T, wt = NULL) {
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

  results = list(coefficients = NULL,
                 se = NULL,
                 fitted.values = NULL,
                 model = NULL,
                 model_x_names = model_x_names,
                 h_x = NULL,
                 gmm_fit = NULL,
                 is.mnar = is.mnar,
                 design_matrix = NULL,  # Store design matrix for analytical derivatives
                 inv_link = inv_link)   # Store inverse link function for analytical gradient

  if(is.mnar){
    alpha_dim = 1 + as.numeric(is.mnar) + ncol(x)

    d = NULL
    if(ncol(h_x1) > 0){
      for(j in 1:ncol(h_x1)) h_x1[, j] = as.factor(h_x1[, j])
      d = model.matrix(lm(rep(1, n)~., data = h_x1))
    }else{
      d = as.matrix(rep(1, n))
    }

    if(ncol(h_x2) == 0){
      h_x2 = NULL
    }else{
      h_x2 = as.matrix(h_x2)
    }

    h_x = cbind(d, h_x2)
    h_dim = ncol(h_x)

    wt <- if (is.null(wt)) 1 else wt

    # Precompute design matrix ONCE (avoid cbind in hot loop)
    design_matrix = cbind(rep(1, n), y, x)

    # Detect link function type once and create optimized inline functions
    # Test at multiple points to distinguish logistic vs probit
    test_val_1 = inv_link(1)
    test_val_0 = inv_link(0)

    # Determine link type:
    # logistic_complement: inv_link(0) = 0.5, inv_link(1) ~ 0.269
    # logistic: inv_link(0) = 0.5, inv_link(1) ~ 0.731
    # probit_complement: inv_link(0) = 0.5, inv_link(1) ~ 0.159 (pnorm(-1))
    # probit: inv_link(0) = 0.5, inv_link(1) ~ 0.841 (pnorm(1))

    # Use inv_link(1) to distinguish: probit values are more extreme
    # logistic(1) ~ 0.731, probit(1) ~ 0.841
    # logistic_complement(1) ~ 0.269, probit_complement(1) ~ 0.159
    is_complement = (test_val_1 < 0.5)
    is_probit = (test_val_1 > 0.8 || test_val_1 < 0.2)

    # Store link type for use in derivatives
    link_type = if (is_probit && is_complement) {
      "probit_complement"
    } else if (is_probit && !is_complement) {
      "probit"
    } else if (!is_probit && is_complement) {
      "logistic_complement"
    } else {
      "logistic"
    }

    # Inlined pi computation (avoids function call overhead)
    # logistic_complement: pi = 1/(1+exp(eta))
    # logistic: pi = exp(eta)/(1+exp(eta)) = 1/(1+exp(-eta))
    # probit: pi = pnorm(eta)
    # probit_complement: pi = pnorm(-eta) = 1 - pnorm(eta)
    compute_pi = if (link_type == "logistic_complement") {
      function(eta) 1 / (1 + exp(eta))
    } else if (link_type == "logistic") {
      function(eta) 1 / (1 + exp(-eta))
    } else if (link_type == "probit") {
      function(eta) pnorm(eta)
    } else {  # probit_complement
      function(eta) pnorm(-eta)
    }

    # Optimized model function using precomputed design_matrix and inlined pi
    model = function(x, y, alpha){
      compute_pi(design_matrix %*% alpha)
    }

    # Precompute r as vector once
    r_vec = as.vector(r)

    # Optimized Phi_alpha: compute eta once and reuse
    Phi_alpha <- function(param) {
      eta = design_matrix %*% param
      pi_vec = compute_pi(eta)
      rw = r_vec / as.vector(pi_vec)
      g.matrix = (rw - 1) * h_x
      return(wt * g.matrix)
    }

    # Analytical gradient of g(alpha)
    # g(alpha) = (r/pi - 1) * h_x, where pi = inv_link(design_matrix %*% alpha)
    # d(g_l)/d(alpha_j) = -r/pi^2 * d(pi)/d(alpha_j) * h_x[,l]
    # For logistic_complement: d(pi)/d(eta) = -pi*(1-pi)
    # For logistic: d(pi)/d(eta) = pi*(1-pi)
    # For probit: d(pi)/d(eta) = dnorm(eta) (phi(eta))
    # For probit_complement: d(pi)/d(eta) = -dnorm(eta)
    dg_alpha <- function(param) {
      eta = as.vector(design_matrix %*% param)
      pi_vec = as.vector(compute_pi(eta))

      # d(g_l)/d(alpha_j) = -r/pi^2 * d(pi)/d(alpha) * h_x[,l]
      # where d(pi)/d(alpha) = d(pi)/d(eta) * d(eta)/d(alpha) = d(pi)/d(eta) * X
      if (link_type == "logistic_complement") {
        # d(pi)/d(eta) = -pi*(1-pi), so d(g)/d(alpha) = -r/pi^2 * (-pi*(1-pi)) * X * h_x
        #              = r*(1-pi)/pi * X * h_x
        common_factor = r_vec * (1 - pi_vec) / pi_vec
      } else if (link_type == "logistic") {
        # d(pi)/d(eta) = pi*(1-pi), so d(g)/d(alpha) = -r/pi^2 * pi*(1-pi) * X * h_x
        #              = -r*(1-pi)/pi * X * h_x
        common_factor = -r_vec * (1 - pi_vec) / pi_vec
      } else if (link_type == "probit") {
        # d(pi)/d(eta) = dnorm(eta), so d(g)/d(alpha) = -r/pi^2 * dnorm(eta) * X * h_x
        common_factor = -r_vec * dnorm(eta) / (pi_vec^2)
      } else {  # probit_complement
        # d(pi)/d(eta) = -dnorm(eta), so d(g)/d(alpha) = -r/pi^2 * (-dnorm(eta)) * X * h_x
        #              = r * dnorm(eta) / pi^2 * X * h_x
        common_factor = r_vec * dnorm(eta) / (pi_vec^2)
      }

      # Build the 3D array: (n, alpha_dim, h_dim)
      Gamma_arr = array(0, dim = c(n, alpha_dim, h_dim))
      for (l in 1:h_dim) {
        Gamma_arr[, , l] = (common_factor * h_x[, l]) * design_matrix
      }
      return(wt * Gamma_arr)
    }

    # Analytical second derivative (Hessian) for Gamma_2 - OPTIMIZED
    # d(Gamma[l,j])/d(alpha_k) = mean(c * X[,k] * X[,j] * h_x[,l])
    # Use crossprod for faster computation
    # For probit, the Hessian is more complex due to dnorm derivative
    d2g_alpha <- function(param) {
      eta = as.vector(design_matrix %*% param)
      pi_vec = as.vector(compute_pi(eta))

      # c_vec is the same as common_factor in dg_alpha
      # For logistic: c = -r*(1-pi)/pi
      # For logistic_complement: c = r*(1-pi)/pi
      # For probit: c = -r * dnorm(eta) / pi^2
      # For probit_complement: c = r * dnorm(eta) / pi^2
      if (link_type == "logistic_complement") {
        c_vec = r_vec * (1 - pi_vec) / pi_vec
      } else if (link_type == "logistic") {
        c_vec = -r_vec * (1 - pi_vec) / pi_vec
      } else if (link_type == "probit") {
        c_vec = -r_vec * dnorm(eta) / (pi_vec^2)
      } else {  # probit_complement
        c_vec = r_vec * dnorm(eta) / (pi_vec^2)
      }

      # Precompute c * design_matrix
      c_X = c_vec * design_matrix  # n x alpha_dim

      # result[row_idx, k] = mean(c_X[,k] * X[,j] * h_x[,l]) = crossprod(c_X, X[,j] * h_x[,l]) / n
      result = matrix(0, h_dim * alpha_dim, alpha_dim)

      for (j in 1:alpha_dim) {
        X_j_h_x = design_matrix[, j] * h_x  # n x h_dim matrix
        # crossprod(c_X, X_j_h_x) gives alpha_dim x h_dim matrix
        block = crossprod(c_X, X_j_h_x) / n
        for (l in 1:h_dim) {
          result[(l-1)*alpha_dim + j, ] = block[, l]
        }
      }
      return(result)
    }

    gmm_fit = private$gmm(Phi_alpha, private$W, n, h_dim, alpha_dim, init, se.fit, dg_alpha, d2g_alpha)

    pi.hat = as.vector(model(x, y, gmm_fit$estimates))

    results$coefficients = gmm_fit$estimates
    results$mu.hat = mean(r/pi.hat*as.matrix(private$y))
    results$se = gmm_fit$se
    results$fitted.values = pi.hat
    results$model = model
    results$h_x = h_x
    results$gmm_fit = gmm_fit
    results$design_matrix = design_matrix  # Reuse precomputed matrix
    results$link_type = link_type  # Store link type for use in EBMR_IPW
  }else{
    glm_fit = glm(r~x, family = binomial())

    model = function(x, y, alpha){
      eta = cbind(rep(1, n), x)%*%alpha
      exp(eta)/(1+exp(eta))
    }

    pi.hat = as.vector(model(x, y, glm_fit$coefficients))
    pi_alpha = pi.hat*(1-pi.hat)*cbind(rep(1, n), x)
    h_x = pi_alpha/(1-pi.hat)
    score = as.vector(r/pi.hat-1)*h_x

    Gamma.hat = (t(-1/(pi.hat*(1-pi.hat))*pi_alpha-as.vector(r-pi.hat)/((pi.hat*(1-pi.hat))^2)*(1-2*pi.hat)*pi_alpha)%*%pi_alpha
                 +t(as.vector(r-pi.hat)/(pi.hat*(1-pi.hat))*(1-2*pi.hat)*pi_alpha)%*%cbind(rep(1, n), x))/n

    gmm_fit = list(
      psi = -Gamma.hat%*%t(score),
      Gamma.hat = Gamma.hat,
      W.hat = diag(length(glm_fit$coefficients)),
      g.matrix = score
    )

    results$coefficients = glm_fit$coefficients
    results$mu.hat = mean(r/pi.hat*as.matrix(private$y))
    results$se = summary(glm_fit)$coefficients[, "Std. Error"]
    results$fitted.values = pi.hat
    results$model = model
    results$h_x = h_x
    results$gmm_fit = gmm_fit
    results$design_matrix = cbind(rep(1, n), x)  # Store design matrix for analytical derivatives
    # For MAR models, use standard logistic link (same as glm with binomial family)
    results$inv_link = function(eta) exp(eta)/(1+exp(eta))  # Standard logistic
    results$link_type = "logistic"  # MAR models use standard logistic
  }

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

ensemble = function(ps.matrix, h_nu, init = NULL, se.fit = T, wt = NULL) {
  # Basic setup
  r = as.matrix(private$r)
  n = private$n
  J = private$J

  result = private$separate_variable_types(h_nu(self$data))
  h_x1 = result$x1
  h_x2 = result$x2

  d = NULL
  if(ncol(h_x1) > 0){
    for(j in 1:ncol(h_x1)) h_x1[, j] = as.factor(h_x1[, j])
    d = model.matrix(lm(rep(1, n)~., data =  h_x1))
  }else{
    d = as.matrix(rep(1, n))
  }

  if(ncol(h_x2) == 0){
    h_x2 = NULL
  }else{
    h_x2 = as.matrix(h_x2)
  }

  h_x = cbind(d, h_x2)
  h_dim = ncol(h_x)

  wt <- if (is.null(wt)) 1 else wt

  # Precompute r as vector once
  r_vec = as.vector(r)

  # Optimized Phi_nu: avoid redundant as.vector calls
  Phi_nu = function(param){
    ps_nu = as.vector(ps.matrix %*% param)
    rw = r_vec / ps_nu
    g.matrix = (rw - 1) * h_x
    return(wt * g.matrix)
  }

  # Analytical gradient of g(nu) for ensemble - OPTIMIZED
  # g(nu) = (r/(ps.matrix %*% nu) - 1) * h_x
  # d(g_l)/d(nu_j) = -r/(ps.matrix %*% nu)^2 * ps.matrix[,j] * h_x[,l]
  #
  # Precompute constant parts outside the function
  dg_nu <- function(param) {
    ps_nu = as.vector(ps.matrix %*% param)
    neg_r_ps2 = -r_vec / (ps_nu * ps_nu)  # n vector

    # Build 3D array: Gamma_arr[i, j, l] = neg_r_ps2[i] * ps.matrix[i,j] * h_x[i,l]
    base_mat = neg_r_ps2 * ps.matrix  # n x J matrix

    # Efficient construction using array recycling
    Gamma_arr = array(0, dim = c(n, J, h_dim))
    for (l in 1:h_dim) {
      Gamma_arr[, , l] = base_mat * h_x[, l]
    }

    return(wt * Gamma_arr)
  }

  # Analytical second derivative (Hessian) for Gamma_2 in ensemble - OPTIMIZED
  # d(Gamma[l,j])/d(nu_k) = mean(dc_factor * ps.matrix[,k] * ps.matrix[,j] * h_x[,l])
  # Output is (h_dim * J) x J matrix
  d2g_nu <- function(param) {
    ps_nu = as.vector(ps.matrix %*% param)
    dc_factor = 2 * r_vec / (ps_nu * ps_nu * ps_nu)  # n vector

    # Precompute dc_factor * ps.matrix
    dc_ps = dc_factor * ps.matrix  # n x J

    # Use crossprod for faster computation
    # result[row_idx, k] = mean(dc_ps[,k] * ps_j * h_x[,l]) = crossprod(dc_ps, ps_j * h_x[,l]) / n
    result = matrix(0, h_dim * J, J)

    for (j in 1:J) {
      ps_j_h_x = ps.matrix[, j] * h_x  # n x h_dim matrix
      # For all l at once: crossprod(dc_ps, ps_j_h_x) gives J x h_dim matrix
      # Each column l gives result[(l-1)*J + j, ] for all k
      block = crossprod(dc_ps, ps_j_h_x) / n  # J x h_dim
      for (l in 1:h_dim) {
        result[(l-1)*J + j, ] = block[, l]
      }
    }
    return(result)
  }

  gmm_fit = private$gmm(Phi_nu, private$W, n, h_dim, J, init, se.fit, dg_nu, d2g_nu)
  # init.mat = diag(J)*0.95
  # min_tmp = 10^8
  # gmm_fit = NULL
  # for(j in 1:J){
  #   gmm_fit_tmp = private$gmm(Phi_nu, private$W, n, h_dim, J, init.mat[j,], se.fit)
  #   if(gmm_fit_tmp$opt$objective < min_tmp){
  #     gmm_fit = gmm_fit_tmp
  #     min_tmp = gmm_fit_tmp$opt$objective
  #   }
  # }

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

EBMR_IPW = function(h_nu, nu_init = c(0.9999, rep(0, J-1)), se.fit = TRUE, true_ps = NULL, wt = NULL) {
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
  # mu.vector = unlist(lapply(ps_fit.list, function(ps_fit) ps_fit$mu.hat))
  alpha.hat = unlist(alpha.list)
  alpha_dim = unlist(lapply(alpha.list, length))
  ps_model.list = lapply(ps_fit.list, function(ps_fit) ps_fit$model)
  ps.matrix = do.call(cbind, lapply(ps_fit.list, function(ps_fit) ps_fit$fitted.values))
  #-----------------------------------------------------------------------------#

  #-----------------------------------------------------------------------------#
  # Ensemble step
  #-----------------------------------------------------------------------------#
  ensemble_fit = private$ensemble(ps.matrix, h_nu, nu_init, se.fit)
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
    # Analytical gradient for dot_pi
    # For logistic_complement: pi = 1/(1+exp(eta)), d(pi)/d(alpha) = -pi*(1-pi)*X
    # For logistic: pi = exp(eta)/(1+exp(eta)), d(pi)/d(alpha) = pi*(1-pi)*X
    # For probit: pi = pnorm(eta), d(pi)/d(alpha) = dnorm(eta)*X
    # For probit_complement: pi = pnorm(-eta), d(pi)/d(alpha) = -dnorm(eta)*X
    dot_pi = matrix(NA, n, sum(alpha_dim))
    for(j in 1:J){
      pi_j = ps_fit.list[[j]]$fitted.values
      design_mat = ps_fit.list[[j]]$design_matrix  # cbind(1, y, x)

      # Get link type from stored value or detect from inv_link
      link_type_j = ps_fit.list[[j]]$link_type
      if (is.null(link_type_j)) {
        # Fallback: detect link type (for backwards compatibility)
        inv_link_j = ps_fit.list[[j]]$inv_link
        test_val = inv_link_j(1)
        is_complement = (test_val < 0.5)
        is_probit = (test_val > 0.8 || test_val < 0.2)
        link_type_j = if (is_probit && is_complement) "probit_complement"
                      else if (is_probit) "probit"
                      else if (is_complement) "logistic_complement"
                      else "logistic"
      }

      # Compute analytical gradient based on link type
      col_idx = (sum(alpha_dim[0:(j-1)])+1):sum(alpha_dim[1:j])
      if (link_type_j == "logistic_complement") {
        # d(pi)/d(alpha) = -pi*(1-pi)*X
        dot_pi[, col_idx] = -design_mat * (pi_j * (1 - pi_j))
      } else if (link_type_j == "logistic") {
        # d(pi)/d(alpha) = pi*(1-pi)*X
        dot_pi[, col_idx] = design_mat * (pi_j * (1 - pi_j))
      } else if (link_type_j == "probit") {
        # d(pi)/d(alpha) = dnorm(eta)*X, where eta = qnorm(pi)
        # For probit, eta = qnorm(pi), so dnorm(eta) = dnorm(qnorm(pi))
        eta_j = qnorm(pi_j)
        dot_pi[, col_idx] = design_mat * dnorm(eta_j)
      } else {  # probit_complement
        # d(pi)/d(alpha) = -dnorm(eta)*X, where pi = pnorm(-eta)
        # eta = -qnorm(pi), so dnorm(-eta) = dnorm(qnorm(pi))
        eta_j = -qnorm(pi_j)
        dot_pi[, col_idx] = -design_mat * dnorm(eta_j)
      }
    }

    dot_W = function(nu){
      nu = as.vector(nu)
      (diag(2*nu)*sum(nu^2)-2*(nu)%*%t(nu^2))/(sum(nu^2)^2)
    }
    dot_W_nu_hat = 0; if(length(nu.hat) > 1) dot_W_nu_hat = dot_W(nu.hat)

    # Precompute common factor used multiple times
    ry_ps_inv2 = as.vector(r*y*((ensemble_ps)^(-2)))

    # t(t(dot_pi)*vec) multiplies each COLUMN of dot_pi by corresponding element of vec
    # (sweep(dot_pi, 2, vec, "*") is equivalent but t(t(X)*vec) is faster for this)
    H_alpha.w = colMeans(t(t(dot_pi)*rep(w.hat, alpha_dim))*ry_ps_inv2)
    w.H_nu = colMeans(ps.matrix%*%t(dot_W_nu_hat) * ry_ps_inv2)

    psi_alpha = do.call(rbind, lapply(ps_fit.list, function(ps_fit) ps_fit$gmm_fit$psi))
    psi_nu = ensemble_fit$gmm_fit$psi

    Gamma_nu = ensemble_fit$gmm_fit$Gamma.hat
    W_nu = ensemble_fit$gmm_fit$W.hat
    h_nu = ensemble_fit$h_x

    # Precompute ps.matrix %*% nu.hat once
    ps_nu = as.vector(ps.matrix%*%nu.hat)
    r_ps_nu_inv2 = as.vector(-r*(ps_nu^(-2)))

    # t(t(dot_pi)*rep(nu.hat, alpha_dim)) multiplies each COLUMN by the corresponding nu.hat element
    Phi_nu.alpha = crossprod(cbind(h_nu)*r_ps_nu_inv2, t(t(dot_pi)*rep(nu.hat, alpha_dim)))/n

    # Use crossprod for t(Gamma_nu)%*%W_nu
    GtW_nu = crossprod(Gamma_nu, W_nu)  # t(Gamma_nu) %*% W_nu
    dot_nu = -solve(GtW_nu %*% Gamma_nu) %*% GtW_nu %*% Phi_nu.alpha
    #--------------------------------------------------------------------------#

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
                nu.hat = ensemble_fit$coefficients,
                w.hat = w.hat,
                ensemble_fit = ensemble_fit
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

EBMR_IPW_with_locally_misspecified_model = function(ps.matrix, perturb_ps, exp_tilt, exp_tilt_x_names, h_nu, nu_init = rep(1/J, J), se.fit = FALSE){
  # Basic setup
  r = as.matrix(private$r)
  y = as.matrix(private$y)
  n = private$n

  #-----------------------------------------------------------------------------#
  # Perturb the designated propensity score model with the exponential tilt model
  #-----------------------------------------------------------------------------#
  J = ncol(ps.matrix)
  ps.matrix[, perturb_ps] = exp_tilt(y, self$data[, exp_tilt_x_names])*ps.matrix[, perturb_ps]
  # mu.vector = apply(as.vector(r*y)/ps.matrix, 2, mean)
  #-----------------------------------------------------------------------------#

  #-----------------------------------------------------------------------------#
  # Ensemble step
  #-----------------------------------------------------------------------------#
  ensemble_fit = private$ensemble(ps.matrix, h_nu, nu_init, se.fit)
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
      x = as.matrix(self$data[ps_fit.list[[j]]$model_x_names])
      dot_pi[, (sum(alpha_dim[0:(j-1)])+1):sum(alpha_dim[1:j])] = jacobian(function(alpha) ps_model.list[[j]](x, y, alpha), alpha.list[[j]])
    }

    dot_W = function(nu){
      nu = as.vector(nu)
      (diag(2*nu)*sum(nu^2)-2*(nu)%*%t(nu^2))/(sum(nu^2)^2)
    }
    dot_W_nu_hat = 0; if(length(nu.hat) > 1) dot_W_nu_hat = dot_W(nu.hat)

    # Precompute common factor used multiple times
    ry_ps_inv2 = as.vector(r*y*((ensemble_ps)^(-2)))

    # Optimized: avoid double transpose
    H_alpha.w = colMeans(dot_pi * rep(w.hat, alpha_dim) * ry_ps_inv2)
    w.H_nu = colMeans(ps.matrix%*%t(dot_W_nu_hat) * ry_ps_inv2)

    psi_alpha = do.call(rbind, lapply(ps_fit.list, function(ps_fit) ps_fit$gmm_fit$psi))
    psi_nu = ensemble_fit$gmm_fit$psi

    Gamma_nu = ensemble_fit$gmm_fit$Gamma.hat
    W_nu = ensemble_fit$gmm_fit$W.hat
    h_nu = ensemble_fit$h_x

    # Precompute ps.matrix %*% nu.hat once
    ps_nu = as.vector(ps.matrix%*%nu.hat)
    r_ps_nu_inv2 = as.vector(-r*(ps_nu^(-2)))

    # Use crossprod for better performance: crossprod(A, B) = t(A) %*% B
    # Optimized: avoid double transpose in second argument
    Phi_nu.alpha = crossprod(cbind(h_nu)*r_ps_nu_inv2, dot_pi * rep(nu.hat, alpha_dim))/n

    # Use crossprod for t(Gamma_nu)%*%W_nu
    GtW_nu = crossprod(Gamma_nu, W_nu)  # t(Gamma_nu) %*% W_nu
    dot_nu = -solve(GtW_nu %*% Gamma_nu) %*% GtW_nu %*% Phi_nu.alpha
    #--------------------------------------------------------------------------#

    mu_ipw.iid = as.vector(t(r/ensemble_ps*y)
                           -(t(H_alpha.w)+t(w.H_nu)%*%dot_nu)%*%psi_alpha
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
