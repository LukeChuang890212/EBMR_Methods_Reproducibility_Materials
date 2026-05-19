#' EBMRAlgorithm Class
#'
#' The `EIPS` class provides an implementation for estimating
#' the coefficients in a propensity score model and related functionalities
#' based on the ensemble method for inverse probability weighting (IPW).
#'
#' @import R6
#'
#' @section Public Methods:
#' The following methods are available in the EIPS class:
#'
#' \code{initialize(y_names, ps_specifications, data, W, wt = NULL, W_nu = NULL)}:
#' Initializes the class and estimates the propensity scores.
#'
#' \code{WangShaoKim2014(formula, h_alpha, inv_link, init = NULL)}:
#' Implements the Wang, Shao, and Kim (2014) approach.
#'
#' \code{EBMR_IPW(h_nu)}:
#' Computes the IPW estimator.
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
#' ps_specifications = list(
#'   formula.list = list(formula1, formula2),
#'   h_alpha.list = list(function(d) d[c("x1","x2")], c("x1")),
#'   inv_link = inv_link_function
#' )
#'
#' ebmr = EIPS$new("y", ps_specifications, data, W)
#' }
#'
#' @export
#'
#' @name EIPS

library(R6)
library(stringr)
library(Matrix)
library(dplyr)
library(numDeriv)

source("./R/GMM_functions.R")
source("./R/Methods.r")
source("./R/Preprocessor.r")
source("./R/Fool_proofing.r")

EIPS <- R6Class("EIPS",
                  public = list(
                    # Public fields (variables)
                    data = NULL,
                    ps_fit.list = list(),
                    WangShaoKim2014 = WangShaoKim2014,
                    EBMR_IPW = EBMR_IPW,
                    EBMR_IPW_regression = EBMR_IPW_regression,
                    EBMR_IPW_with_locally_misspecified_model = EBMR_IPW_with_locally_misspecified_model,

                    # Constructor to initialize fields
                    initialize = function(y_names, ps_specifications, data, W, wt = NULL, W_nu = NULL) {
                      self$data = private$check_data(y_names, data)
                      private$r = self$data$r
                      private$y = self$data[y_names]
                      private$n = nrow(self$data)
                      private$J = length(ps_specifications$formula.list)
                      private$W = W
                      private$W_nu = if (!is.null(W_nu)) W_nu else W
                      private$wt = wt
                      # Get outcome variable name (new in Fast2)
                      outcome <- ps_specifications$outcome
                      # Get optimizer and cond_threshold (support scalar or per-model list)
                      opt_spec <- ps_specifications$optimizer
                      cond_spec <- ps_specifications$cond_threshold

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

                        # Resolve per-model optimizer and cond_threshold
                        optimizer_j <- if (is.list(opt_spec)) opt_spec[[j]]
                                       else if (!is.null(opt_spec)) opt_spec
                                       else "L-BFGS-B"
                        cond_j <- if (is.list(cond_spec)) cond_spec[[j]]
                                  else if (!is.null(cond_spec)) cond_spec
                                  else 1e8

                        if (is.null(wt)) {
                          self$ps_fit.list[[j]] <- self$WangShaoKim2014(
                            formula = formula,
                            outcome = outcome,
                            h_alpha = h_alpha,
                            inv_link = inv_link,
                            init = alpha_init,
                            optimizer = optimizer_j,
                            cond_threshold = cond_j
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
                            optimizer = optimizer_j,
                            cond_threshold = cond_j
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
                    W_nu = NULL,
                    wt = NULL,
                    # private methods
                    check_data = check_data,
                    parse_formula = parse_formula,
                    build_design_matrix = build_design_matrix,
                    separate_variable_types = separate_variable_types,
                    ensemble = ensemble,
                    gmm = gmm
                  )

)
