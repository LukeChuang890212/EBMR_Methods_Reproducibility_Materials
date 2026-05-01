#' EBMRAlgorithm Class
#'
#' The `EBMRAlgorithmFast` class provides an implementation for estimating
#' the coefficients in a propensity score model and related functionalities
#' based on the ensemble method for inverse probability weighting (IPW).
#' This class allows users to perform propensity score fitting and estimation
#' for the specified data and propensity score specifications.
#'
#' @import R6
#'
#' @section Public Methods:
#' The following methods are available in the EBMRAlgorithmFast class:
#'
#' \code{initialize(y_names, ps_specifications, data)}:
#' Initializes the class and estimates the propensity scores based on the provided formulas,
#' covariates to be balanced, inverse link function, and the specified data.
#'
#' \code{WangShaoKim2014(formula, h_alpha, inv_link, init = NULL)}:
#' Implements the Wang, Shao, and Kim (2014) approach for estimating propensity score models in the presence of missing-not-at-random (MNAR) data.
#' \code{h_alpha} can be a function (e.g., \code{function(data) data[c("x1","x2")]}) or a character vector of column names (backwards compatible).
#'
#' \code{EBMR_IPW(h_nu, true_ps = NULL)}:
#' Computes the Inverse Probability Weighting (IPW) estimator for the population mean \eqn{\mu_0},
#' using the propensity scores estimated by the ensemble method. It also computes standard errors and other related
#' quantities for the estimator, including the estimator when the true propensity score is provided.
#'
#' @section Private Fields:
#' \describe{
#'   \item{\code{r}}{The response indicator.}
#'   \item{\code{y}}{The outcome variable(s) used in the model.}
#'   \item{\code{n}}{The number of observations.}
#'   \item{\code{ps_fit.list}}{List of fitted propensity score models.}
#' }
#'
#' @examples
#' \dontrun{
#' # Define propensity score specifications (for illustration)
#' ps_specifications = list(
#'   formula.list = list(formula1, formula2),
#'   h_alpha.list = list(function(d) d[c("x1","x2")], c("x1")),  # function or character vector
#'   inv_link = inv_link_function
#' )
#'
#' # Initialize the EBMRAlgorithmFast class
#' ebmr = EBMRAlgorithmFast$new(y_names = "outcome", ps_specifications = ps_specifications, data = dataset)
#' }
#'
#' @export
#'
#' @name EBMRAlgorithmFast3

library(R6)
library(stringr)
library(Matrix)
library(dplyr)
library(numDeriv)
library(DEoptim)
library(parallel)

# source("./R/WangShaoKim2014.r")
source("./R/GMM_functions.R")
source("./R/Methods.r")
source("./R/Preprocessor.r")
source("./R/Fool_proofing.r")

EBMRAlgorithmFast3 <- R6Class("EBMRAlgorithmFast3",
                  public = list(
                    # Public fields (variables)
                    data = NULL,
                    ps_fit.list = list(),
                    WangShaoKim2014 = WangShaoKim2014,
                    EBMR_IPW = EBMR_IPW,
                    EBMR_IPW_regression = EBMR_IPW_regression,
                    EBMR_IPW_with_locally_misspecified_model = EBMR_IPW_with_locally_misspecified_model,

                    # Constructor to initialize fields
                    initialize = function(y_names, ps_specifications, data, W, wt = NULL, method = "GN") {
                      self$data = private$check_data(y_names, data)
                      private$r = self$data$r
                      private$y = self$data[y_names]
                      private$n = nrow(self$data)
                      private$J = length(ps_specifications$formula.list)
                      private$W = W
                      private$wt = wt
                      private$method = method

                      # Get outcome variable name (new in Fast2)
                      outcome <- ps_specifications$outcome

                      # Support both new h_alpha.list and old h_x_names.list
                      h_alpha_list <- if (!is.null(ps_specifications$h_alpha.list)) {
                        ps_specifications$h_alpha.list
                      } else {
                        ps_specifications$h_x_names.list
                      }

                      for (j in 1:private$J) {
                        formula <- ps_specifications$formula.list[[j]]
                        h_alpha <- h_alpha_list[[j]]
                        alpha_init <- ps_specifications$alpha_init.list[[j]]
                        inv_link <- ps_specifications$inv_link

                        if (is.null(wt)) {
                          self$ps_fit.list[[j]] <- self$WangShaoKim2014(
                            formula = formula,
                            outcome = outcome,
                            h_alpha = h_alpha,
                            inv_link = inv_link,
                            init = alpha_init,
                            method = method
                          )
                        } else {
                          self$ps_fit.list[[j]] <- self$WangShaoKim2014(
                            formula = formula,
                            outcome = outcome,
                            h_alpha = h_alpha,
                            inv_link = inv_link,
                            init = alpha_init,
                            se.fit = FALSE,
                            wt = wt,
                            method = method
                          )
                        }
                      }
                    }
                  ),

                  private = list(
                    # private fields (variables)
                    r = NULL,
                    y = NULL,
                    n = NULL,
                    J = NULL,
                    W = NULL,
                    wt = NULL,
                    method = NULL,

                    # private methods
                    check_data = check_data,
                    parse_formula = parse_formula,
                    build_design_matrix = build_design_matrix,
                    separate_variable_types = separate_variable_types,
                    ensemble = ensemble,
                    gmm = gmm
                  )

)
