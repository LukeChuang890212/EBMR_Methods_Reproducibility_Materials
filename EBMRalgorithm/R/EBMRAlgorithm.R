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

source("./R/WangShaoKim2014.r")
source("./R/Methods.r")
source("./R/Preprocessor.r")
source("./R/Fool_proofing.r")

EBMRAlgorithm <- R6Class("EBMRAlgorithm",
                   public = list(
                     # Public fields (variables)
                     data = NULL,
                     ps_fit.list = list(),

                     # Constructor to initialize fields
                     initialize = function(y_names, ps_specifications, data) {
                       self$data = private$check_data(y_names, data)
                       private$r = self$data$r
                       private$y = self$data[y_names]
                       private$n = nrow(self$data)
                       private$wt = rexp(private$n)
                       private$wt = private$wt/sum(private$wt)*private$n

                       J = length(ps_specifications$formula.list)
                       for(j in 1:J){
                         formula = ps_specifications$formula.list[[j]]
                         h_x_names = ps_specifications$h_x_names.list[[j]]
                         inv_link = ps_specifications$inv_link
                         self$ps_fit.list[[j]] = self$WangShaoKim2014(formula, h_x_names, inv_link)
                       }
                       # print(self$ps_fit.list[[3]]$coefficients)
                       # print(self$ps_fit.list[[3]]$se)
                     },

                     # public fields (variables)
                     WangShaoKim2014 = WangShaoKim2014,
                     EBMR_IPW = EBMR_IPW,
                     WangShaoKim2014_perturb = WangShaoKim2014_perturb,
                     EBMR_IPW_perturb = EBMR_IPW_perturb,
                     EBMR_IPW_with_locally_misspecified_model = EBMR_IPW_with_locally_misspecified_model
                   ),

                   private = list(
                     # private fields (variables)
                     r = NULL,
                     y = NULL,
                     n = NULL,
                     wt = NULL,

                     # private methods
                     check_data = check_data,
                     parse_formula = parse_formula,
                     separate_variable_types = separate_variable_types,
                     estimate_nu = estimate_nu,
                     estimate_nu_perturb = estimate_nu_perturb
                   )

)

###################################
#  ---------- Test ----------------
###################################

# Chuang2023.1.1 = function(n){
#   z1 = rbinom(n, size = 1, prob = 0.3)
#   z2 = rnorm(n, mean = 0, sd = 2)
#
#   u1 = rbinom(n, size = 1, prob = 0.7)
#   u2 = rnorm(n, mean = 0, sd = 2)
#
#   m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
#   response.prob = function(y, u1, u2)  1/(1+exp(-0.2+0.1*y+0.1*u1+0.1*u2))
#
#   y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)
#   mean(y)
#
#   r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
#   mean(r)
#
#   mean(y[r == 1]); mean(y[r == 0]);
#
#   dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
#   return(dat)
# }
#
# data = Chuang2023.1.1(1000)
# r = data$r
# y = data$y
# u1 = data$u1
# u2 = data$u2
# z1 = data$z1
# z2 = data$z2
#
# J = 3
# dat = data
# alpha.true = c(-0.2, 0.1, 0.1, 0.1)
# true.pi.model = function(y, u1, u2, r, n, alpha.true) 1/(1+exp(cbind(rep(1, 1000), y, u1, u2)%*%alpha.true))
# true.pi = true.pi.model(y, u1, u2, r, 1000, alpha.true)
#
# source("E:\\Other computers\\我的電腦\\MNAR-Simulation\\MNAR_2023\\ChuangChao2023_SM.r")
# n = sum(r)
# N = 1000
# h.list = list(function(u1, u2, z1, z2) list(cbind(as.factor(u1), as.factor(z1)), cbind(u2, z2)),
#               function(u1, u2, z1, z2) list(cbind(as.factor(u1), as.factor(z1)), cbind(z2)),
#               function(u1, u2, z1, z2) list(cbind(as.factor(z1)), cbind(u2, z2)))
#
# propensity.list = list(list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
#                             w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
#                             model.y = function(y) y,
#                             model.x1.names = c("u1"),
#                             model.x2.names =c("u2")),
#                        list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
#                             w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
#                             model.y = function(y) y,
#                             model.x1.names = c("u1"),
#                             model.x2.names = NULL),
#                        list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
#                             w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
#                             model.y = function(y) y,
#                             model.x1.names = NULL,
#                             model.x2.names = c("u2")))
#
# start = Sys.time()
# pi.fit.list = list()
# for(j in 1:J){
#   pi.fit.list[[j]] = Wang2014.1(auxilliary = h.list[[j]](u1, u2, z1, z2),
#                                 model.y = propensity.list[[j]]$model.y(y),
#                                 model.x1.names = propensity.list[[j]]$model.x1.names,
#                                 model.x2.names = propensity.list[[j]]$model.x2.names,
#                                 w = propensity.list[[j]]$w, w.prime = propensity.list[[j]]$w.prime)
# }
# pi.fit.list[[3]]$theta.hat
# pi.fit.list[[3]]$se
# auxilliary = list(cbind(as.factor(u1)), cbind(u2), y)
# est.res1 = ChuangChao2023(pi.fit.list, auxilliary, family, ortho = TRUE, true.pi)
# est1 =  unlist(est.res1[1:4])
# est1
# Sys.time() - start
#
# # source("C:\\Users\\stat-pc\\Desktop\\NTHU_Research\\EBMR_Methods_Reproducibility_Materials\\EBMRalgorithm\\R\\Methods.r")
# ps_specifications = list(
#   formula.list = list(
#     r ~ o(y) + u1 + u2,
#     r ~ o(y) + u1,
#     r ~ o(y) + u2
#   ),
#   h_x_names.list = list(
#     c("u1", "u2", "z1", "z2"),
#     c("u1", "z1", "z2"),
#     c("u2", "z1", "z2")
#   ),
#   inv_link = function(eta) 1/(1+exp(eta))
# )
#
# start = Sys.time()
# ebmr <- EBMRAlgorithm$new("y", ps_specifications, data)
# result = ebmr$EBMR_IPW(h_x_names = c("u1", "u2"), true_ps = true.pi)
# result = unlist(result[1:4])
# result
# Sys.time() - start
