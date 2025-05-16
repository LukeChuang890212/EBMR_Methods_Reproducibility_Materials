#' EBMRAlgorithm Class
#'
#' The `EBMRAlgorithm` class provides an implementation for estimating
#' the coefficients in a propensity score model and related functionalities
#' based on the ensemble method for inverse probability weighting (IPW).
#' This class allows users to perform propensity score fitting and estimation
#' for the specified data and propensity score specifications.
#'
#' @import R6
#'
#' @section Public Methods:
#' The following methods are available in the EBMRAlgorithm class:
#'
#' \code{initialize(y_names, ps_specifications, data)}:
#' Initializes the class and estimates the propensity scores based on the provided formulas,
#' covariates to be balanced, inverse link function, and the specified data.
#'
#' \code{WangShaoKim2014(formula, h_x_names, inv_link, init = NULL)}:
#' Implements the Wang, Shao, and Kim (2014) approach for estimating propensity score models in the presence of missing-not-at-random (MNAR) data.
#'
#' \code{EBMR_IPW(h_x_names, true_ps = NULL)}:
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
#'   h_x_names.list = list(h_x1, h_x2),
#'   inv_link = inv_link_function
#' )
#'
#' # Initialize the EBMRAlgorithm class
#' ebmr = EBMRAlgorithm$new(y_names = "outcome", ps_specifications = ps_specifications, data = dataset)
#' }
#'
#' @export
#'
#' @name EBMRAlgorithm

library(R6)
library(stringr)
library(Matrix)
library(dplyr)
library(numDeriv)

# source("./R/WangShaoKim2014.r")
source("./R/GMM_functions.r")
source("./R/test_Methods.r")
source("./R/Preprocessor.r")
source("./R/Fool_proofing.r")

EBMRAlgorithm <- R6Class("EBMRAlgorithm",
                   public = list(
                     # Public fields (variables)
                     data = NULL,
                     ps_fit.list = list(),
                     WangShaoKim2014 = WangShaoKim2014,
                     EBMR_IPW = EBMR_IPW,
                     EBMR_IPW_with_locally_misspecified_model = EBMR_IPW_with_locally_misspecified_model,

                     # Constructor to initialize fields
                     initialize = function(y_names, ps_specifications, data, W, wt = NULL) {
                       self$data = private$check_data(y_names, data)
                       # self$EBMR_IPW = ifelse(is.perturb, EBMR_IPW_perturb, EBMR_IPW)
                       private$r = self$data$r
                       private$y = self$data[y_names]
                       private$n = nrow(self$data)
                       private$J = length(ps_specifications$formula.list)
                       private$W = W

                       for(j in 1:private$J){
                         formula = ps_specifications$formula.list[[j]]
                         h_x_names = ps_specifications$h_x_names.list[[j]]
                         alpha_init = ps_specifications$alpha_init.list[[j]]
                         inv_link = ps_specifications$inv_link
                         if(is.null(wt)){
                           self$ps_fit.list[[j]] = self$WangShaoKim2014(formula, h_x_names, inv_link, alpha_init)
                         }else{
                           self$ps_fit.list[[j]] = self$WangShaoKim2014(formula, h_x_names, inv_link, alpha_init, se.fit = F, wt)
                         }
                       }

                       # print(self$ps_fit.list[[3]]$coefficients)
                       # print(self$ps_fit.list[[3]]$se)
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

                     # private methods
                     check_data = check_data,
                     parse_formula = parse_formula,
                     separate_variable_types = separate_variable_types,
                     ensemble = ensemble,
                     # Phi_alpha = Phi_alpha,
                     # Phi_nu = Phi_nu,
                     # G = G,
                     # t_Gamma_i = t_Gamma_i,
                     # Gamma = Gamma,
                     # Gamma_2 = Gamma_2,
                     # obj = obj,
                     gmm = gmm
                   )

)

