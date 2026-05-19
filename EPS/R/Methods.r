#' Wang, Shao, and Kim (2014) Method for Propensity Score Estimation
#'
#' Implements the Wang, Shao, and Kim (2014) approach for estimating propensity
#' score models in the presence of missing-not-at-random (MNAR) data.
#'
#' @import numDeriv
#'
#' @param formula A formula specifying the propensity score model.
#' Examples:
#' \itemize{
#'   \item Simple MNAR: `r ~ y + health + father`
#'   \item MNAR with interactions: `r ~ y + y:health + y:father + health + father`
#'   \item MAR model: `r ~ health + father`
#' }
#' @param outcome Character string specifying the outcome variable name (e.g., "y").
#' Terms containing this variable are automatically identified as MNAR terms.
#' If NULL, the model is treated as MAR.
#' @param h_alpha A function or character vector specifying the auxiliary variables
#' to be balanced. If a function, it should take a data frame and return a data
#' frame or matrix with named columns (e.g., \code{function(data) data[c("x1","x2")]}).
#' If a character vector, it is treated as column names (backwards compatible with
#' the old \code{h_x_names} parameter).
#' @param inv_link An inverse link function applied to the linear predictor.
#' @param init (Optional) Initial values for the model parameters. Defaults to NULL.
#' @param se.fit Logical, whether to compute standard errors. Defaults to TRUE.
#' @param wt (Optional) Observation weights.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{\code{coefficients}}{Estimated model coefficients.}
#'   \item{\code{fitted.values}}{Predicted propensity scores.}
#'   \item{\code{se}}{Standard errors of the estimated coefficients.}
#'   \item{\code{model}}{The propensity score model function.}
#'   \item{\code{is.mnar}}{Logical indicator of whether the model is MNAR.}
#'   \item{\code{h_x}}{Matrix of covariates to be balanced.}
#'   \item{\code{design_matrix}}{The design matrix used in estimation.}
#'   \item{\code{link_type}}{The type of link function used.}
#' }
#'
#' @examples
#' \dontrun{
#' ebmr <- EBMRAlgorithm$new(data = data)
#'
#' # MNAR model with interactions
#' result <- ebmr$WangShaoKim2014(
#'   formula = r ~ y + y:health + y:father + health + father,
#'   outcome = "y",
#'   h_alpha = c("health", "father"),
#'   inv_link = function(eta) 1 / (1 + exp(eta))
#' )
#'
#' # MAR model (using function form)
#' result <- ebmr$WangShaoKim2014(
#'   formula = r ~ health + father,
#'   outcome = NULL,
#'   h_alpha = function(data) data[c("health", "father")],
#'   inv_link = function(eta) 1 / (1 + exp(eta))
#' )
#' }
#'
#' @references Wang, Shao, & Kim (2014). "An instrumental variable approach
#' for identification and estimation with nonignorable nonresponse."

WangShaoKim2014 = function(formula, outcome = NULL, h_alpha, inv_link,
                           init = NULL, se.fit = TRUE, wt = NULL,
                           optimizer = "L-BFGS-B", cond_threshold = 1e8) {
  # Backwards compatibility: if h_alpha is a character vector, convert to function
  if (is.character(h_alpha)) {
    h_alpha_names <- h_alpha
    h_alpha <- function(data) data[h_alpha_names]
  }

  # Parse formula using new method (auto-identifies outcome-related terms)
  parsed <- private$parse_formula(formula, outcome)
  r <- as.matrix(self$data[parsed$r_name])
  n <- nrow(self$data)

  # Build design matrix using model.matrix (handles interactions properly)
  dm <- private$build_design_matrix(self$data, parsed, include_intercept = TRUE)
  design_matrix <- dm$design
  alpha_dim <- ncol(design_matrix)

  # For backwards compatibility: extract x_names for model_x_names
  model_x_names <- parsed$x_terms

  result_sep <- private$separate_variable_types(h_alpha(self$data))
  h_x1 <- result_sep$x1
  h_x2 <- result_sep$x2

  is.mnar <- parsed$is_mnar

  results <- list(
    coefficients = NULL,
    se = NULL,
    fitted.values = NULL,
    model = NULL,
    model_x_names = model_x_names,
    h_x = NULL,
    gmm_fit = NULL,
    is.mnar = is.mnar,
    design_matrix = NULL,
    inv_link = inv_link
  )

  if (is.mnar) {
    # alpha_dim already computed from design_matrix

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

    # design_matrix already built by build_design_matrix() above

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

    # Pi computation using R's built-in plogis/pnorm for numerical stability.
    # These handle extreme eta values without overflow and maintain consistency
    # with derivative computations (no clamping discontinuity).
    compute_pi = if (link_type == "logistic_complement") {
      function(eta) plogis(-eta)   # = 1/(1+exp(eta))
    } else if (link_type == "logistic") {
      function(eta) plogis(eta)    # = 1/(1+exp(-eta))
    } else if (link_type == "probit") {
      function(eta) pnorm(eta)
    } else {  # probit_complement
      function(eta) pnorm(-eta)
    }

    # Optimized model function using precomputed design_matrix and inlined pi
    # Note: x and y args kept for backwards compatibility but not used
    model = function(x = NULL, y = NULL, alpha) {
      eta <- as.vector(design_matrix %*% alpha)
      eta <- pmin(pmax(eta, -20), 20)  # clamp eta to prevent overflow
      compute_pi(eta)
    }

    # Precompute r as vector once
    r_vec = as.vector(r)

    # Model-level cache: share eta/pi_vec between Phi_alpha and Gamma_alpha.
    # Both functions are called at the same param in every grad/conv_grad step.
    # Caching avoids duplicate design_matrix %*% param and compute_pi calls.
    mc <- new.env(parent = emptyenv())
    mc$param   <- NULL
    mc$pi_vec  <- NULL

    get_pi_vec <- function(param) {
      if (!is.null(mc$param) && identical(param, mc$param)) return(mc$pi_vec)
      mc$param  <- param
      eta <- as.vector(design_matrix %*% param)
      eta <- pmin(pmax(eta, -20), 20)  # clamp eta to prevent overflow
      mc$pi_vec <- as.vector(compute_pi(eta))
      mc$pi_vec
    }

    # Optimized Phi_alpha: uses shared cache for pi_vec
    Phi_alpha <- function(param) {
      pi_vec = get_pi_vec(param)
      rw = r_vec / pi_vec
      wt * (rw - 1) * h_x
    }

    # Analytical gradient of g(alpha)
    # g(alpha) = (r/pi - 1) * h_x, where pi = f(eta) (raw model)
    # d(g_l)/d(alpha_j) = -r/pi^2 * f'(eta) * X[,j] * h_x[,l]
    dg_alpha <- function(param) {
      pi_vec = get_pi_vec(param)   # shared cache: no duplicate compute_pi

      # common_factor = wt * (-r/pi^2 * f'(eta)), with pi = f(eta)
      # wt is folded into common_factor to avoid broadcasting issues with 3D arrays
      if (link_type == "logistic_complement") {
        common_factor = wt * r_vec * (1 - pi_vec) / pi_vec
      } else if (link_type == "logistic") {
        common_factor = wt * (-r_vec) * (1 - pi_vec) / pi_vec
      } else if (link_type == "probit") {
        eta = as.vector(design_matrix %*% param)
        common_factor = wt * (-r_vec) * dnorm(eta) / (pi_vec^2)
      } else {  # probit_complement
        eta = as.vector(design_matrix %*% param)
        common_factor = wt * r_vec * dnorm(eta) / (pi_vec^2)
      }

      # Build the 3D array: (n, alpha_dim, h_dim)
      Gamma_arr = array(0, dim = c(n, alpha_dim, h_dim))
      for (l in 1:h_dim) {
        Gamma_arr[, , l] = (common_factor * h_x[, l]) * design_matrix
      }
      return(Gamma_arr)
    }

    # Analytical second derivative for raw model: pi = f(eta)
    # c_i = -r_i/pi_i^2 * f'(eta_i), dc_i/d(eta_i) needed for Hessian
    # With pi = f: logistic cases simplify to dc/deta = r*(1-pi)/pi
    d2g_alpha <- function(param) {
      pi_vec = get_pi_vec(param)   # shared cache
      eta = as.vector(design_matrix %*% param)  # needed for probit dc_deta

      if (link_type == "logistic_complement" || link_type == "logistic") {
        # Both logistic variants: dc/deta = r*(1-pi)/pi
        dc_deta = r_vec * (1 - pi_vec) / pi_vec
      } else if (link_type == "probit") {
        phi = dnorm(eta)
        dc_deta = r_vec * (eta * phi / (pi_vec^2) + 2 * phi^2 / (pi_vec^3))
      } else {  # probit_complement
        phi = dnorm(eta)
        dc_deta = r_vec * (eta * phi / (pi_vec^2) - 2 * phi^2 / (pi_vec^3))
      }

      # result[(l-1)*alpha_dim + j, k] = mean(dc_deta * X[,k] * X[,j] * h_x[,l])
      dc_X = dc_deta * design_matrix  # n x alpha_dim
      result = matrix(0, h_dim * alpha_dim, alpha_dim)
      for (j in 1:alpha_dim) {
        Xj_h = design_matrix[, j] * h_x  # n x h_dim
        block = crossprod(dc_X, Xj_h) / n  # alpha_dim x h_dim
        for (l in 1:h_dim) {
          result[(l-1)*alpha_dim + j, ] = block[, l]
        }
      }
      return(result)
    }

    # Fast Gamma function: returns E[dg/dalpha] = Gamma_mat (h_dim x alpha_dim) directly
    # using crossprod, avoiding the full (n x alpha_dim x h_dim) 3D array allocation.
    # Mathematically: Gamma_mat[l,j] = mean_i(common_factor_i * h_x[i,l] * X[i,j])
    #                                 = (1/n) * crossprod(h_x * common_factor, X)[l,j]
    # This is the fast path used in grad() and conv_grad() during optimization;
    # dg_alpha is still used for SE computation (needs full 3D array for H_2n_2).
    Gamma_alpha <- function(param) {
      pi_vec = get_pi_vec(param)   # shared cache: reuses eta/pi_vec from Phi_alpha
      # wt folded into common_factor (same as dg_alpha)
      if (link_type == "logistic_complement") {
        common_factor = wt * r_vec * (1 - pi_vec) / pi_vec
      } else if (link_type == "logistic") {
        common_factor = wt * (-r_vec) * (1 - pi_vec) / pi_vec
      } else if (link_type == "probit") {
        eta = as.vector(design_matrix %*% param)
        common_factor = wt * (-r_vec) * dnorm(eta) / (pi_vec^2)
      } else {  # probit_complement
        eta = as.vector(design_matrix %*% param)
        common_factor = wt * r_vec * dnorm(eta) / (pi_vec^2)
      }
      # h_dim x alpha_dim matrix via BLAS crossprod
      crossprod(h_x * common_factor, design_matrix) / n
    }

    if (is.null(init)) {
      init <- rep(0, alpha_dim)
    }

    gmm_fit = private$gmm(Phi_alpha, private$W, n, h_dim, alpha_dim, init, se.fit, dg_alpha, d2g_alpha,
                          Gamma_direct = Gamma_alpha, optimizer = optimizer,
                          cond_threshold = cond_threshold)

    pi.hat = as.vector(model(alpha = gmm_fit$estimates))

    results$coefficients = gmm_fit$estimates
    results$mu.hat = mean(r / pi.hat * as.matrix(private$y))
    results$se = gmm_fit$se
    results$fitted.values = pi.hat
    results$model = model
    results$h_x = h_x
    results$gmm_fit = gmm_fit
    results$opt = gmm_fit$opt
    results$design_matrix = design_matrix
    results$link_type = link_type

  } else {
    # MAR model: use standard glm
    # Build formula from x_terms
    if (length(parsed$x_terms) > 0) {
      glm_formula <- as.formula(
        paste(parsed$r_name, "~", paste(parsed$x_terms, collapse = " + "))
      )
    } else {
      glm_formula <- as.formula(paste(parsed$r_name, "~ 1"))
    }
    glm_fit <- glm(glm_formula, data = self$data, family = binomial())

    # Model function for MAR
    model = function(x = NULL, y = NULL, alpha) {
      eta <- design_matrix %*% alpha
      exp(eta) / (1 + exp(eta))
    }

    pi.hat <- as.vector(model(alpha = glm_fit$coefficients))
    pi_alpha <- pi.hat * (1 - pi.hat) * design_matrix
    h_x <- pi_alpha / (1 - pi.hat)
    score <- as.vector(r / pi.hat - 1) * h_x

    Gamma.hat <- (t(-1 / (pi.hat * (1 - pi.hat)) * pi_alpha -
                    as.vector(r - pi.hat) / ((pi.hat * (1 - pi.hat))^2) *
                    (1 - 2 * pi.hat) * pi_alpha) %*% pi_alpha +
                  t(as.vector(r - pi.hat) / (pi.hat * (1 - pi.hat)) *
                    (1 - 2 * pi.hat) * pi_alpha) %*% design_matrix) / n

    gmm_fit <- list(
      psi = -Gamma.hat %*% t(score),
      Gamma.hat = Gamma.hat,
      W.hat = diag(length(glm_fit$coefficients)),
      g.matrix = score
    )

    results$coefficients = glm_fit$coefficients
    results$mu.hat = mean(r / pi.hat * as.matrix(private$y))
    results$se = summary(glm_fit)$coefficients[, "Std. Error"]
    results$fitted.values = pi.hat
    results$model = model
    results$h_x = h_x
    results$gmm_fit = gmm_fit
    results$design_matrix = design_matrix
    results$inv_link = function(eta) exp(eta) / (1 + exp(eta))
    results$link_type = "logistic"
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

ensemble = function(ps.matrix, h_nu, init = NULL, se.fit = T, wt = NULL,
                    nu_optimizer = "L-BFGS-B", nu_cond_threshold = 1e8) {
  # Basic setup
  r = as.matrix(private$r)
  n = private$n
  J = ncol(ps.matrix)

  # J=1: no ensemble needed
  if (J == 1) {
    h_x = as.matrix(h_nu(self$data))
    h_dim = ncol(h_x)
    return(list(
      coefficients = 1,
      h_x = h_x,
      gmm_fit = list(
        estimates = 1,
        psi = matrix(0, h_dim, n),
        Gamma.hat = matrix(0, h_dim, 1),
        W.hat = diag(h_dim),
        opt = list(objective = 0, iterations = 0, converged = TRUE, final_grad_norm = 0)
      )
    ))
  }

  h_x2 = as.matrix(h_nu(self$data))

  # Use intercept + covariates as moment conditions (same structure as h_alpha)
  # PS differences are not used because they depend on y which is missing when r=0
  d = as.matrix(rep(1, n))

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
    # wt folded into neg_r_ps2 to avoid broadcasting issues with 3D arrays
    neg_r_ps2 = wt * (-r_vec / (ps_nu * ps_nu))  # n vector

    # Build 3D array: Gamma_arr[i, j, l] = neg_r_ps2[i] * ps.matrix[i,j] * h_x[i,l]
    base_mat = neg_r_ps2 * ps.matrix  # n x J matrix

    # Efficient construction using array recycling
    Gamma_arr = array(0, dim = c(n, J, h_dim))
    for (l in 1:h_dim) {
      Gamma_arr[, , l] = base_mat * h_x[, l]
    }

    return(Gamma_arr)
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

  # Fast Gamma function for ensemble: returns E[dg_nu/dnu] (h_dim x J) directly
  # Gamma_mat[l,j] = mean_i(neg_r_ps2[i] * ps.matrix[i,j] * h_x[i,l])
  #                = (1/n) * crossprod(h_x * neg_r_ps2, ps.matrix)[l,j]
  Gamma_nu_direct <- function(param) {
    ps_nu = as.vector(ps.matrix %*% param)
    neg_r_ps2 = wt * (-r_vec / (ps_nu * ps_nu))
    crossprod(h_x * neg_r_ps2, ps.matrix) / n  # h_dim x J
  }

  # Starting value: use init if provided, otherwise uniform
  uniform_init <- if (!is.null(init)) init else rep(1/J, J)

  gmm_fit <- private$gmm(Phi_nu, private$W_nu, n, h_dim, J, uniform_init, se.fit, dg_nu, d2g_nu,
                          Gamma_direct = Gamma_nu_direct,
                          optimizer = nu_optimizer, cond_threshold = nu_cond_threshold)

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
#' @param h_nu A function that takes a data frame and returns a matrix of
#'        covariates for the ensemble step (e.g., \code{function(d) cbind(x1=d$x1, x2=d$x2)}).
#' @param true_ps (optional) A vector of true propensity scores. If provided, the IPW estimator
#'        will also be computed using the true propensity scores.
#' @param type Character string specifying the type of IPW estimator: \code{"HT"} for
#'        Horvitz-Thompson (default) or \code{"Hajek"} for the Hajek ratio estimator.
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
#' ipw_estimates <- ebmr$EBMR_IPW(h_nu = function(d) cbind(x1=d$x1, x2=d$x2), true_ps = true_ps_data)
#' print(ipw_estimates)
#' }

EBMR_IPW = function(h_nu, model_indices = NULL, nu_init = rep(1/J, J), se.fit = TRUE, true_ps = NULL, wt = NULL, type = c("HT", "Hajek"),
                    nu_optimizer = "L-BFGS-B", nu_cond_threshold = 1e8) {
  # Basic setup
  r = as.matrix(private$r)
  y = as.matrix(private$y)
  n = private$n
  type = match.arg(type)

  # Subset PS models if model_indices specified
  if (is.null(model_indices)) {
    ps_fit.list = self$ps_fit.list
  } else {
    ps_fit.list = self$ps_fit.list[model_indices]
  }

  #-----------------------------------------------------------------------------#

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
  ensemble_fit = private$ensemble(ps.matrix, h_nu, nu_init, se.fit, wt = wt,
                                  nu_optimizer = nu_optimizer,
                                  nu_cond_threshold = nu_cond_threshold)
  nu.hat = ensemble_fit$coefficients
  w.hat = nu.hat^2/sum(nu.hat^2)
  ensemble_ps = ps.matrix%*%w.hat
  #-----------------------------------------------------------------------------#

  #-----------------------------------------------------------------------------#
  # IPW estimator for the population mean mu_0 with propensity score being estimated
  # by the methods of Wang, Shao and Kim (2014).
  #-----------------------------------------------------------------------------#
  if (type == "HT") {
    mu_ipw = ifelse(is.null(wt), mean(r/ensemble_ps*y), mean(wt*r/ensemble_ps*y))
  } else {
    A_hat = ifelse(is.null(wt), mean(r/ensemble_ps*y), mean(wt*r/ensemble_ps*y))
    B_hat = ifelse(is.null(wt), mean(r/ensemble_ps), mean(wt*r/ensemble_ps))
    mu_ipw = A_hat / B_hat
  }
  se_ipw = NA
  if(se.fit){
    #--------------------------------------------------------------------------#
    # Compute necessary quantities to estimate the influence function:
    # \psi(\bm{\alpha}_*, \bm{\nu}_*)
    #--------------------------------------------------------------------------#
    # Analytical gradient for dot_pi
    # Raw model: pi = f(eta), d(pi)/d(alpha) = f'(eta)*X
    dot_pi = matrix(NA, n, sum(alpha_dim))
    for(j in 1:J){
      pi_j = ps_fit.list[[j]]$fitted.values
      design_mat = ps_fit.list[[j]]$design_matrix

      link_type_j = ps_fit.list[[j]]$link_type
      if (is.null(link_type_j)) {
        inv_link_j = ps_fit.list[[j]]$inv_link
        test_val = inv_link_j(1)
        is_complement = (test_val < 0.5)
        is_probit = (test_val > 0.8 || test_val < 0.2)
        link_type_j = if (is_probit && is_complement) "probit_complement"
                      else if (is_probit) "probit"
                      else if (is_complement) "logistic_complement"
                      else "logistic"
      }

      col_idx = (sum(alpha_dim[0:(j-1)])+1):sum(alpha_dim[1:j])
      if (link_type_j == "logistic_complement") {
        # f'(eta) = -pi*(1-pi)
        dot_pi[, col_idx] = -design_mat * (pi_j * (1 - pi_j))
      } else if (link_type_j == "logistic") {
        # f'(eta) = pi*(1-pi)
        dot_pi[, col_idx] = design_mat * (pi_j * (1 - pi_j))
      } else if (link_type_j == "probit") {
        eta_j = qnorm(pi_j)
        dot_pi[, col_idx] = design_mat * dnorm(eta_j)
      } else {  # probit_complement
        eta_j = -qnorm(pi_j)
        dot_pi[, col_idx] = -design_mat * dnorm(eta_j)
      }
    }

    dot_W = function(nu){
      nu = as.vector(nu)
      (diag(2*nu)*sum(nu^2)-2*(nu)%*%t(nu^2))/(sum(nu^2)^2)
    }
    # Precompute common factor used multiple times
    if (type == "HT") {
      ry_ps_inv2 = as.vector(r*y*((ensemble_ps)^(-2)))
    } else {
      ry_ps_inv2 = as.vector(r*(y - mu_ipw)*((ensemble_ps)^(-2)) / B_hat)
    }

    H_alpha.w = colMeans(t(t(dot_pi)*rep(w.hat, alpha_dim))*ry_ps_inv2)
    psi_alpha = do.call(rbind, lapply(ps_fit.list, function(ps_fit) ps_fit$gmm_fit$psi))

    # Base influence function (alpha estimation uncertainty only)
    if (type == "HT") {
      mu_ipw.iid = as.vector(t(r/ensemble_ps*y) - t(H_alpha.w)%*%psi_alpha)
    } else {
      base_term = as.vector(r/ensemble_ps*(y - mu_ipw) / B_hat)
      mu_ipw.iid = as.vector(t(base_term) - t(H_alpha.w)%*%psi_alpha)
    }

    # Add ensemble (nu) estimation uncertainty when J > 1
    if (J > 1) {
      dot_W_nu_hat = dot_W(nu.hat)
      w.H_nu = colMeans(ps.matrix%*%t(dot_W_nu_hat) * ry_ps_inv2)
      psi_nu = ensemble_fit$gmm_fit$psi
      Gamma_nu = ensemble_fit$gmm_fit$Gamma.hat
      W_nu = ensemble_fit$gmm_fit$W.hat
      h_nu = ensemble_fit$h_x

      ps_nu = as.vector(ps.matrix%*%nu.hat)
      r_ps_nu_inv2 = as.vector(-r*(ps_nu^(-2)))
      Phi_nu.alpha = crossprod(cbind(h_nu)*r_ps_nu_inv2, t(t(dot_pi)*rep(nu.hat, alpha_dim)))/n

      GtW_nu = crossprod(Gamma_nu, W_nu)
      dot_nu = -solve(GtW_nu %*% Gamma_nu) %*% GtW_nu %*% Phi_nu.alpha

      mu_ipw.iid = mu_ipw.iid - as.vector(
        (t(w.H_nu)%*%dot_nu)%*%psi_alpha + t(w.H_nu)%*%psi_nu)
    }

    se_ipw = sqrt(var(mu_ipw.iid)/n)
  }
  #-----------------------------------------------------------------------------#

  #-----------------------------------------------------------------------------#
  # IPW estimator for the population mean mu_0 with known propensity score.
  #-----------------------------------------------------------------------------#
  mu_ipw.true = NA
  se_ipw.true = NA
  if(!is.null(true_ps)){
    if (type == "HT") {
      mu_ipw.true = mean(r/true_ps*y)
      mu_ipw.true.iid = as.vector(r/true_ps*y)
    } else {
      A_true = mean(r/true_ps*y)
      B_true = mean(r/true_ps)
      mu_ipw.true = A_true / B_true
      mu_ipw.true.iid = as.vector(r/true_ps*(y - mu_ipw.true) / B_true)
    }
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

EBMR_IPW_with_locally_misspecified_model = function(ps.matrix, perturb_ps, exp_tilt, exp_tilt_x_names, h_nu, nu_init = rep(1/J, J), se.fit = FALSE, type = c("HT", "Hajek")){
  # Basic setup
  r = as.matrix(private$r)
  y = as.matrix(private$y)
  n = private$n
  type = match.arg(type)

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
  if (type == "HT") {
    mu_ipw = mean(r/ensemble_ps*y)
  } else {
    A_hat = mean(r/ensemble_ps*y)
    B_hat = mean(r/ensemble_ps)
    mu_ipw = A_hat / B_hat
  }
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
    if (type == "HT") {
      ry_ps_inv2 = as.vector(r*y*((ensemble_ps)^(-2)))
    } else {
      ry_ps_inv2 = as.vector(r*(y - mu_ipw)*((ensemble_ps)^(-2)) / B_hat)
    }

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

    if (type == "HT") {
      mu_ipw.iid = as.vector(t(r/ensemble_ps*y)
                             -(t(H_alpha.w)+t(w.H_nu)%*%dot_nu)%*%psi_alpha
                             -t(w.H_nu)%*%psi_nu)
    } else {
      base_term = as.vector(r/ensemble_ps*(y - mu_ipw) / B_hat)
      mu_ipw.iid = as.vector(t(base_term)
                             -(t(H_alpha.w)+t(w.H_nu)%*%dot_nu)%*%psi_alpha
                             -t(w.H_nu)%*%psi_nu)
    }
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

#' Compute the IPW Estimator for Regression Coefficients
#'
#' This method estimates regression coefficients \eqn{\theta} by solving the
#' inverse probability weighted estimating equation:
#' \deqn{\frac{1}{n}\sum_{i=1}^{n}\frac{R_i}{\pi_i}\Phi_\theta(Y_i, X_i; \theta) = 0}
#' where \eqn{\Phi_\theta} is the estimating function (score function) for the
#' regression model, and \eqn{\pi_i} is the ensemble propensity score estimator.
#'
#' @param h_nu A function that takes the data frame and returns a matrix of
#'        covariates for the ensemble step (same as in EBMR_IPW).
#' @param reg_formula A formula specifying the regression model (e.g., y ~ x1 + x2).
#' @param family A character string specifying the regression family.
#'        Options: "gaussian" (linear regression), "binomial" (logistic regression).
#'        Default is "gaussian".
#' @param Phi_theta (Optional) A custom estimating function. If NULL, uses the
#'        canonical estimating function for the specified family.
#'        Should be a function(theta, Y, X) returning an n x p matrix where p = length(theta).
#' @param dPhi_theta (Optional) A custom gradient function for Phi_theta.
#'        If NULL, computed numerically using numDeriv::jacobian.
#'        Should be a function(theta, Y, X) returning the n x p x p array of derivatives.
#' @param theta_init (Optional) Initial values for theta. If NULL, uses OLS/GLM estimates.
#' @param nu_init (Optional) Initial values for nu. Default is c(0.9999, rep(0, J-1)).
#' @param se.fit Logical, whether to compute standard errors. Default is TRUE.
#' @param true_ps (Optional) A vector of true propensity scores for comparison.
#' @param wt (Optional) Observation weights.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{\code{theta_hat}}{The estimated regression coefficients.}
#'   \item{\code{se_theta}}{Standard errors of the estimated coefficients.}
#'   \item{\code{cov_theta}}{Estimated covariance matrix for theta.}
#'   \item{\code{ps.matrix}}{Matrix of estimated propensity scores.}
#'   \item{\code{nu.hat}}{The estimated ensemble coefficients.}
#'   \item{\code{w.hat}}{The estimated weights for propensity score models.}
#'   \item{\code{convergence}}{Convergence status from the optimization.}
#' }
#'
#' @details
#' The variance of \eqn{\hat{\theta}} is estimated using the sandwich formula:
#' \deqn{Var(\hat{\theta}) = A^{-1} B (A^{-1})^T}
#' where:
#' \itemize{
#'   \item \eqn{A = \frac{1}{n}\sum_{i=1}^{n}\frac{R_i}{\pi_i}\frac{\partial \Phi_\theta}{\partial \theta}}
#'   \item \eqn{B = Var\left(\frac{1}{n}\sum_{i=1}^{n}\frac{R_i}{\pi_i}\Phi_\theta\right)}
#' }
#' The variance B accounts for the estimation uncertainty in the propensity scores
#' using the influence function approach from EBMR_IPW.
#'
#' @examples
#' \dontrun{
#' ebmr <- EBMRAlgorithm$new(data = data)
#' # Linear regression
#' result <- ebmr$EBMR_IPW_regression(
#'   h_nu = function(d) cbind(d$x1, d$x2),
#'   reg_formula = y ~ x1 + x2,
#'   family = "gaussian"
#' )
#'
#' # Logistic regression
#' result <- ebmr$EBMR_IPW_regression(
#'   h_nu = function(d) cbind(d$x1, d$x2),
#'   reg_formula = y ~ x1 + x2,
#'   family = "binomial"
#' )
#' }

EBMR_IPW_regression = function(h_nu, reg_formula, family = "gaussian",
                                Phi_theta = NULL, dPhi_theta = NULL,
                                theta_init = NULL, nu_init = NULL,
                                se.fit = TRUE, true_ps = NULL, wt = NULL) {
  # Basic setup
  r = as.matrix(private$r)
  y = as.matrix(private$y)
  n = private$n
  J = private$J
  ps_fit.list = self$ps_fit.list

  # Set default nu_init
  if (is.null(nu_init)) {
    nu_init = c(0.9999, rep(0, J-1))
  }

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
  # Ensemble step (same as EBMR_IPW)
  #-----------------------------------------------------------------------------#
  ensemble_fit = private$ensemble(ps.matrix, h_nu, nu_init, se.fit)
  nu.hat = ensemble_fit$coefficients
  w.hat = nu.hat^2/sum(nu.hat^2)
  ensemble_ps = as.vector(ps.matrix %*% w.hat)
  #-----------------------------------------------------------------------------#

  #-----------------------------------------------------------------------------#
  # Build regression design matrix
  #-----------------------------------------------------------------------------#
  # Parse the regression formula
  reg_terms = terms(reg_formula, data = self$data)
  Y = model.response(model.frame(reg_formula, data = self$data))
  X = model.matrix(reg_formula, data = self$data)
  p = ncol(X)  # Number of regression parameters
  #-----------------------------------------------------------------------------#

  #-----------------------------------------------------------------------------#
  # Define estimating function Phi_theta if not provided
  #-----------------------------------------------------------------------------#
  if (is.null(Phi_theta)) {
    if (family == "gaussian") {
      # Linear regression: Phi_theta = X^T(Y - X*theta)
      # Returns n x p matrix where each row is X[i,] * (Y[i] - X[i,] %*% theta)
      Phi_theta = function(theta, Y, X) {
        resid = as.vector(Y - X %*% theta)
        return(X * resid)  # n x p matrix
      }
    } else if (family == "binomial") {
      # Logistic regression: Phi_theta = X^T(Y - pi(X*theta))
      # where pi(eta) = exp(eta)/(1+exp(eta))
      Phi_theta = function(theta, Y, X) {
        eta = as.vector(X %*% theta)
        pi_val = 1 / (1 + exp(-eta))
        resid = as.vector(Y) - pi_val
        return(X * resid)  # n x p matrix
      }
    } else {
      stop("family must be 'gaussian' or 'binomial', or provide custom Phi_theta")
    }
  }
  #-----------------------------------------------------------------------------#

  #-----------------------------------------------------------------------------#
  # Define gradient of Phi_theta if not provided
  #-----------------------------------------------------------------------------#
  if (is.null(dPhi_theta)) {
    if (family == "gaussian") {
      # For linear: d(Phi_theta)/d(theta) = -X^T X (constant)
      # Returns p x p matrix (averaged over observations)
      dPhi_theta = function(theta, Y, X) {
        return(-X)  # d(residual)/d(theta) = -X, so d(X*resid)/d(theta) needs outer product
      }
      # Actually, we need the Jacobian: d(Phi[i,j])/d(theta[k]) = -X[i,j] * X[i,k]
      # This is an n x p x p array, but we need the mean
      dPhi_theta_mean = function(theta, Y, X) {
        # Mean of d(Phi)/d(theta) = -t(X) %*% X / n
        return(-crossprod(X, X) / nrow(X))
      }
    } else if (family == "binomial") {
      # For logistic: d(Phi_theta)/d(theta) = -X^T diag(pi*(1-pi)) X
      dPhi_theta_mean = function(theta, Y, X) {
        eta = as.vector(X %*% theta)
        pi_val = 1 / (1 + exp(-eta))
        W_diag = pi_val * (1 - pi_val)
        return(-crossprod(X, X * W_diag) / nrow(X))
      }
    } else {
      # Use numerical differentiation
      dPhi_theta_mean = function(theta, Y, X) {
        numDeriv::jacobian(function(th) colMeans(Phi_theta(th, Y, X)), theta)
      }
    }
  } else {
    # User provided dPhi_theta
    dPhi_theta_mean = function(theta, Y, X) {
      dPhi_arr = dPhi_theta(theta, Y, X)
      if (length(dim(dPhi_arr)) == 3) {
        # n x p x p array, take mean over first dimension
        return(apply(dPhi_arr, c(2, 3), mean))
      } else {
        # Already p x p
        return(dPhi_arr)
      }
    }
  }
  #-----------------------------------------------------------------------------#

  #-----------------------------------------------------------------------------#
  # Get initial theta if not provided
  #-----------------------------------------------------------------------------#
  if (is.null(theta_init)) {
    # Use complete case GLM for initial values
    cc_idx = as.vector(r) == 1
    if (family == "gaussian") {
      cc_fit = lm(reg_formula, data = self$data[cc_idx, ])
      theta_init = coef(cc_fit)
    } else if (family == "binomial") {
      cc_fit = glm(reg_formula, data = self$data[cc_idx, ], family = binomial())
      theta_init = coef(cc_fit)
    } else {
      theta_init = rep(0, p)
    }
  }
  #-----------------------------------------------------------------------------#

  #-----------------------------------------------------------------------------#
  # Solve the IPW estimating equation
  # Find theta such that: (1/n) sum_i R_i/pi_i * Phi_theta(Y_i, X_i; theta) = 0
  #-----------------------------------------------------------------------------#
  r_vec = as.vector(r)
  wt_vec = if (is.null(wt)) rep(1, n) else wt

  # IPW weighted mean of Phi_theta
  ipw_Phi_mean = function(theta) {
    Phi_mat = Phi_theta(theta, Y, X)  # n x p matrix
    ipw_weights = wt_vec * r_vec / ensemble_ps
    return(colSums(Phi_mat * ipw_weights) / n)
  }

  # Objective function: sum of squared estimating equations
  obj_fn = function(theta) {
    eq = ipw_Phi_mean(theta)
    return(sum(eq^2))
  }

  # Gradient of objective function (for faster optimization)
  grad_obj_fn = function(theta) {
    eq = ipw_Phi_mean(theta)
    # d(obj)/d(theta) = 2 * eq^T * d(eq)/d(theta)
    # d(eq)/d(theta) = (1/n) sum_i R_i/pi_i * d(Phi)/d(theta)
    A_hat = dPhi_theta_mean(theta, Y, X)
    # Adjust A_hat for IPW weighting
    Phi_mat = Phi_theta(theta, Y, X)
    ipw_weights = wt_vec * r_vec / ensemble_ps

    # Numerical approximation for weighted mean of gradient
    A_ipw = matrix(0, p, p)
    for (k in 1:p) {
      for (j in 1:p) {
        if (family == "gaussian") {
          A_ipw[j, k] = -sum(ipw_weights * X[, j] * X[, k]) / n
        } else if (family == "binomial") {
          eta = as.vector(X %*% theta)
          pi_val = 1 / (1 + exp(-eta))
          W_diag = pi_val * (1 - pi_val)
          A_ipw[j, k] = -sum(ipw_weights * X[, j] * X[, k] * W_diag) / n
        }
      }
    }
    return(2 * as.vector(t(A_ipw) %*% eq))
  }

  # Optimize
  opt_result = optim(theta_init, obj_fn, method = "BFGS",
                     control = list(maxit = 1000, reltol = 1e-10))

  theta_hat = opt_result$par
  names(theta_hat) = colnames(X)
  #-----------------------------------------------------------------------------#

  #-----------------------------------------------------------------------------#
  # Compute standard errors using sandwich formula
  #-----------------------------------------------------------------------------#
  se_theta = rep(NA, p)
  cov_theta = NULL

  if (se.fit) {
    # A = (1/n) sum_i R_i/pi_i * d(Phi)/d(theta)|_{theta=theta_hat}
    # For linear: A = -(1/n) sum_i R_i/pi_i * X_i X_i^T
    # For logistic: A = -(1/n) sum_i R_i/pi_i * pi_i(1-pi_i) * X_i X_i^T

    ipw_weights = wt_vec * r_vec / ensemble_ps

    if (family == "gaussian") {
      A_hat = -crossprod(X * sqrt(ipw_weights), X * sqrt(ipw_weights)) / n
    } else if (family == "binomial") {
      eta = as.vector(X %*% theta_hat)
      pi_val = 1 / (1 + exp(-eta))
      W_diag = pi_val * (1 - pi_val)
      A_hat = -crossprod(X * sqrt(ipw_weights * W_diag), X * sqrt(ipw_weights * W_diag)) / n
    } else {
      # Numerical approximation
      A_hat = dPhi_theta_mean(theta_hat, Y, X)
    }

    # A_inv = inverse of A (the "bread")
    A_inv = tryCatch(solve(A_hat), error = function(e) {
      warning("A matrix is singular, using pseudo-inverse")
      MASS::ginv(A_hat)
    })

    #-------------------------------------------------------------------------#
    # Compute B = Var((1/n) sum_i R_i/pi_i * Phi_theta)
    # This requires accounting for PS estimation uncertainty
    #-------------------------------------------------------------------------#

    # Phi at theta_hat
    Phi_hat = Phi_theta(theta_hat, Y, X)  # n x p matrix

    # Similar to EBMR_IPW variance calculation, but for each component of Phi
    # We need the influence function for the IPW mean of Phi

    # Analytical gradient for dot_pi (same as EBMR_IPW)
    dot_pi = matrix(NA, n, sum(alpha_dim))
    for(j in 1:J){
      pi_j = ps_fit.list[[j]]$fitted.values
      design_mat = ps_fit.list[[j]]$design_matrix

      link_type_j = ps_fit.list[[j]]$link_type
      if (is.null(link_type_j)) {
        inv_link_j = ps_fit.list[[j]]$inv_link
        test_val = inv_link_j(1)
        is_complement = (test_val < 0.5)
        is_probit = (test_val > 0.8 || test_val < 0.2)
        link_type_j = if (is_probit && is_complement) "probit_complement"
                      else if (is_probit) "probit"
                      else if (is_complement) "logistic_complement"
                      else "logistic"
      }

      col_idx = (sum(alpha_dim[0:(j-1)])+1):sum(alpha_dim[1:j])
      if (link_type_j == "logistic_complement") {
        dot_pi[, col_idx] = -design_mat * (pi_j * (1 - pi_j))
      } else if (link_type_j == "logistic") {
        dot_pi[, col_idx] = design_mat * (pi_j * (1 - pi_j))
      } else if (link_type_j == "probit") {
        eta_j = qnorm(pi_j)
        dot_pi[, col_idx] = design_mat * dnorm(eta_j)
      } else {
        eta_j = -qnorm(pi_j)
        dot_pi[, col_idx] = -design_mat * dnorm(eta_j)
      }
    }

    dot_W = function(nu){
      nu = as.vector(nu)
      (diag(2*nu)*sum(nu^2)-2*(nu)%*%t(nu^2))/(sum(nu^2)^2)
    }
    dot_W_nu_hat = 0
    if(length(nu.hat) > 1) dot_W_nu_hat = dot_W(nu.hat)

    # For each component l of Phi, compute the influence function
    # psi_Phi[i, l] = R_i/pi_i * Phi_l - mu_Phi_l - adjustment for PS estimation

    psi_alpha = do.call(rbind, lapply(ps_fit.list, function(ps_fit) ps_fit$gmm_fit$psi))
    psi_nu = ensemble_fit$gmm_fit$psi

    Gamma_nu = ensemble_fit$gmm_fit$Gamma.hat
    W_nu = ensemble_fit$gmm_fit$W.hat
    h_nu_mat = ensemble_fit$h_x

    # Precompute ps.matrix %*% nu.hat once
    ps_nu = as.vector(ps.matrix %*% nu.hat)
    r_ps_nu_inv2 = as.vector(-r_vec * (ps_nu^(-2)))

    # Phi_nu.alpha for ensemble (same structure as EBMR_IPW)
    Phi_nu.alpha = crossprod(cbind(h_nu_mat) * r_ps_nu_inv2, t(t(dot_pi) * rep(nu.hat, alpha_dim))) / n

    GtW_nu = crossprod(Gamma_nu, W_nu)
    dot_nu = -solve(GtW_nu %*% Gamma_nu) %*% GtW_nu %*% Phi_nu.alpha

    # Influence function for each Phi component
    # mu_Phi_ipw.iid[i, l] = R_i/pi_i * Phi[i,l]
    #                       - H_alpha_Phi[l,] %*% psi_alpha[,i]
    #                       - H_nu_Phi[l,] %*% (dot_nu %*% psi_alpha[,i] + psi_nu[,i])

    mu_Phi_ipw = colMeans(Phi_hat * ipw_weights)  # p-vector

    # H_alpha_Phi[l, .] = E[R*Phi[l]/pi^2 * d(pi)/d(alpha) * w]
    # H_nu_Phi[l, .] = E[R*Phi[l]/pi^2 * ps.matrix %*% d(w)/d(nu)]

    rPhi_ps_inv2 = r_vec * (ensemble_ps^(-2))  # n-vector

    # For each Phi component l, compute H_alpha and H_nu
    H_alpha_Phi = matrix(0, p, sum(alpha_dim))
    H_nu_Phi = matrix(0, p, J)

    for (l in 1:p) {
      rPhi_l_ps_inv2 = rPhi_ps_inv2 * Phi_hat[, l]
      H_alpha_Phi[l, ] = colMeans(t(t(dot_pi) * rep(w.hat, alpha_dim)) * rPhi_l_ps_inv2)
      H_nu_Phi[l, ] = colMeans(ps.matrix %*% t(dot_W_nu_hat) * rPhi_l_ps_inv2)
    }

    # Compute influence function matrix (n x p)
    # iid_Phi[i, l] = R_i/pi_i * Phi[i,l]
    #               - H_alpha_Phi[l,] %*% psi_alpha[,i]
    #               - H_nu_Phi[l,] %*% (dot_nu %*% psi_alpha[,i] + psi_nu[,i])

    # psi_alpha is (sum(alpha_dim) x n), psi_nu is (J x n)
    adj_alpha = t(H_alpha_Phi %*% psi_alpha)  # n x p
    adj_nu_via_alpha = t(H_nu_Phi %*% dot_nu %*% psi_alpha)  # n x p
    adj_nu = t(H_nu_Phi %*% psi_nu)  # n x p

    iid_Phi = (Phi_hat * ipw_weights) - adj_alpha - adj_nu_via_alpha - adj_nu

    # B = Var(mean(iid_Phi)) = Cov(iid_Phi) / n
    B_hat = cov(iid_Phi) / n

    # Sandwich: Var(theta_hat) = A_inv %*% B %*% t(A_inv)
    cov_theta = A_inv %*% B_hat %*% t(A_inv)
    se_theta = sqrt(diag(cov_theta))
    names(se_theta) = colnames(X)
  }
  #-----------------------------------------------------------------------------#

  #-----------------------------------------------------------------------------#
  # Results with true PS (if provided)
  #-----------------------------------------------------------------------------#
  theta_hat_true = NULL
  se_theta_true = NULL

  if (!is.null(true_ps)) {
    ipw_weights_true = as.vector(wt_vec * r_vec / true_ps)

    ipw_Phi_mean_true = function(theta) {
      Phi_mat = Phi_theta(theta, Y, X)
      return(colSums(Phi_mat * ipw_weights_true) / n)
    }

    obj_fn_true = function(theta) {
      eq = ipw_Phi_mean_true(theta)
      return(sum(eq^2))
    }

    opt_result_true = optim(theta_init, obj_fn_true, method = "BFGS",
                            control = list(maxit = 1000, reltol = 1e-10))
    theta_hat_true = opt_result_true$par
    names(theta_hat_true) = colnames(X)

    # Simple variance for true PS (no PS estimation uncertainty)
    if (se.fit) {
      if (family == "gaussian") {
        A_hat_true = -crossprod(X * sqrt(ipw_weights_true), X * sqrt(ipw_weights_true)) / n
      } else if (family == "binomial") {
        eta = as.vector(X %*% theta_hat_true)
        pi_val = 1 / (1 + exp(-eta))
        W_diag = pi_val * (1 - pi_val)
        A_hat_true = -crossprod(X * sqrt(ipw_weights_true * W_diag),
                                 X * sqrt(ipw_weights_true * W_diag)) / n
      }
      A_inv_true = solve(A_hat_true)

      Phi_hat_true = Phi_theta(theta_hat_true, Y, X)
      iid_Phi_true = Phi_hat_true * ipw_weights_true
      B_hat_true = cov(iid_Phi_true) / n

      cov_theta_true = A_inv_true %*% B_hat_true %*% t(A_inv_true)
      se_theta_true = sqrt(diag(cov_theta_true))
      names(se_theta_true) = colnames(X)
    }
  }
  #-----------------------------------------------------------------------------#

  result = list(
    theta_hat = theta_hat,
    se_theta = se_theta,
    cov_theta = cov_theta,
    theta_hat_true = theta_hat_true,
    se_theta_true = se_theta_true,
    ps.matrix = ps.matrix,
    nu.hat = ensemble_fit$coefficients,
    w.hat = w.hat,
    ensemble_fit = ensemble_fit,
    convergence = opt_result$convergence,
    family = family,
    reg_formula = reg_formula
  )

  return(result)
}
