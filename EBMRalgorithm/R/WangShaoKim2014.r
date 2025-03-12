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

WangShaoKim2014 = function(formula, h_x_names, inv_link, init = NULL) {
  # Basic setup
  result = private$parse_formula(formula)
  r = as.matrix(self$data[result$r_names])
  y = as.matrix(self$data[result$y_names])
  x = as.matrix(self$data[result$x_names])
  n = private$n
  model_x_names = colnames(x)

  # result = separate_variable_types(x)
  # model_x1 = result$x1
  # model_x2 = result$x2
  # model_x1_names = result$x1_names
  # model_x2_names = result$x2_names

  result = private$separate_variable_types(self$data[h_x_names])
  h_x1 = result$x1
  h_x2 = result$x2

  # h_x1 = h[[1]]; h_x2 = h[[2]];
  # if(!is.null(x1)) x1 = as.matrix(x1)
  # if(!is.null(x2)) x2 = as.matrix(x2)
  #
  # model_x1 = NULL; model_x2 = NULL;
  # if(!is.null(model_x1_names)) model_x1 = as.matrix(dat[model_x1_names])
  # if(!is.null(model_x2_names)) model_x2 = as.matrix(dat[model_x2_names])

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

  # r = as.matrix(r)
  # y = as.matrix(y)
  # x = as.matrix(x)
  if(continuous_dim == 0){
    h_x2 = NULL
  }else{
    h_x2 = as.matrix(h_x2)
  }
  # h_x1 = as.matrix(h_x1)
  # h_x2 = as.matrix(h_x2)

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

  W = function(g.matrix){
    return(solve(t(g.matrix)%*%g.matrix/n))
  }

  Gamma = function(alpha){
    Gamma.arr = array(NA, dim = c(h_dim, n, alpha_dim))
    for(l in 1:h_dim){
      Gamma.arr[l,,] = jacobian(function(alpha) g(alpha)[, l], alpha)
    }
    return(apply(Gamma.arr, c(1, 3), mean))
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
                 model = model,
                 model_x_names = model_x_names,
                 h_x = cbind(d, h_x2))

  return(results)
}

WangShaoKim2014_perturb = function(formula, h_x_names, inv_link, init = NULL) {
  # Basic setup
  result = private$parse_formula(formula)
  r = as.matrix(self$data[result$r_names])
  y = as.matrix(self$data[result$y_names])
  x = as.matrix(self$data[result$x_names])
  n = private$n
  # wt = rexp(n)
  model_x_names = colnames(x)

  # result = separate_variable_types(x)
  # model_x1 = result$x1
  # model_x2 = result$x2
  # model_x1_names = result$x1_names
  # model_x2_names = result$x2_names

  result = private$separate_variable_types(self$data[h_x_names])
  h_x1 = result$x1
  h_x2 = result$x2

  # h_x1 = h[[1]]; h_x2 = h[[2]];
  # if(!is.null(x1)) x1 = as.matrix(x1)
  # if(!is.null(x2)) x2 = as.matrix(x2)
  #
  # model_x1 = NULL; model_x2 = NULL;
  # if(!is.null(model_x1_names)) model_x1 = as.matrix(dat[model_x1_names])
  # if(!is.null(model_x2_names)) model_x2 = as.matrix(dat[model_x2_names])

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

  # r = as.matrix(r)
  # y = as.matrix(y)
  # x = as.matrix(x)
  if(continuous_dim == 0){
    h_x2 = NULL
  }else{
    h_x2 = as.matrix(h_x2)
  }
  # h_x1 = as.matrix(h_x1)
  # h_x2 = as.matrix(h_x2)

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
    return(rexp(n)*g.matrix)
  }

  G = function(g.matrix){
    return(matrix(apply(g.matrix, 2, mean), h_dim, 1))
  }

  W = function(g.matrix){
    return(solve(t(g.matrix)%*%g.matrix/n))
  }

  Gamma = function(alpha){
    Gamma.arr = array(NA, dim = c(h_dim, n, alpha_dim))
    for(l in 1:h_dim){
      Gamma.arr[l,,] = jacobian(function(alpha) g(alpha)[, l], alpha)
    }
    return(apply(Gamma.arr, c(1, 3), mean))
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
                 model = model,
                 model_x_names = model_x_names,
                 h_x = cbind(d, h_x2))

  return(results)
}


# parse_formula <- function(formula) {
#   # Convert the formula to a string
#   formula <- as.character(formula)
#
#   # Extract the left-hand side of the formula (everything before the ~)
#   lhs <- formula[2]
#
#   # Isolate r
#   r <- lhs
#
#   # Extract the right-hand side of the formula (everything after the ~)
#   rhs <- formula[3]
#
#   # Extract the variables inside o() using regular expressions
#   y <- str_extract(rhs, "o\\(([^)]+)\\)")   # Variables in o()
#
#   # Clean up: Remove 'o()' and split the variables (there's only one in this case)
#   y <- gsub("o\\(|\\)", "", y)
#   y <- unlist(strsplit(y, split = "\\+"))
#
#   # Extract the rest of the variables (excluding o())
#   x <- gsub("o\\([^)]+\\)\\s*\\+?", "", rhs)
#
#   # Clean up: Split the rest of the variables by '+' and trim whitespace
#   x <- unlist(strsplit(x, split = "\\+"))
#   x <- trimws(x)
#
#   # Return a list with r, y, x
#   return(list(r = r, y = y, x = x))
# }

# # Function to separate continuous and discrete variable
# separate_variable_types <- function(x) {
#
#   # Initialize empty lists to store continuous and discrete columns
#   x1 <- list()
#   x2 <- list()
#
#   col_names = colnames(x)
#   x1_names = c()
#   x2_names = c()
#
#   # Iterate through each column in the matrix 'x'
#   for (i in 1:ncol(x)) {
#     column <- x[, i]
#
#     # Check if the column is numeric
#     if (is.numeric(column)) {
#       # If the column has more than 10 unique values, it's treated as continuous
#       if (length(unique(column)) > 10) {
#         x2[[colnames(x)[i]]] <- column
#         x2_names = c(x2_names, col_names[i])
#       } else {
#         x1[[colnames(x)[i]]] <- column
#         x1_names = c(x1_names, col_names[i])
#       }
#     } else {
#       # For non-numeric columns, treat them as discrete
#       x1[[colnames(x)[i]]] <- column
#       x1_names = c(x1_names, col_names[i])
#     }
#   }
#
#   # Convert the lists to matrices (if you want matrices, otherwise keep them as lists)
#   x1 <- as.data.frame(x1)
#   x2 <- as.data.frame(x2)
#
#   # Return the results
#   return(list(x1 = x1, x2 = x2, x1_names = x1_names, x2_names = x2_names))
# }

# WangShaoKim2014 = function(formula, h_x, inv_link, data, init = NULL){
#   # Parse the formula
#   result = parse_formula(formula)
#   r = data[result$r]
#   y = data[result$y]
#   x = data[result$x]
#   model_x_names = colnames(x)
#
#   n = nrow(data)
#
#   # result = separate_variable_types(x)
#   # model_x1 = result$x1
#   # model_x2 = result$x2
#   # model_x1_names = result$x1_names
#   # model_x2_names = result$x2_names
#
#   result = separate_variable_types(h_x)
#   h_x1 = result$x1
#   h_x2 = result$x2
#
#   model = function(x, y, alpha){
#     inv_link(cbind(rep(1, n), y, x)%*%alpha)
#   }
#   w = function(x, y, alpha) 1/model(x, y, alpha)
#
#   # h_x1 = h[[1]]; h_x2 = h[[2]];
#   # if(!is.null(x1)) x1 = as.matrix(x1)
#   # if(!is.null(x2)) x2 = as.matrix(x2)
#   #
#   # model_x1 = NULL; model_x2 = NULL;
#   # if(!is.null(model_x1_names)) model_x1 = as.matrix(dat[model_x1_names])
#   # if(!is.null(model_x2_names)) model_x2 = as.matrix(dat[model_x2_names])
#
#   is.mnar = ifelse(ncol(y) == 0, 0, 1)
#   alpha_dim = 1 + is.mnar + ncol(x)
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
#   y = as.matrix(y)
#   x = as.matrix(x)
#   h_x1 = as.matrix(h_x1)
#   h_x2 = as.matrix(h_x2)
#
#   g = function(alpha){
#     g.matrix = matrix(NA, n, h_dim)
#     rw = r*w(x, y, alpha)
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
#   Gamma = function(alpha){
#     Gamma.arr = array(NA, dim = c(h_dim, n, alpha_dim))
#     for(l in 1:h_dim){
#       Gamma.arr[l,,] = jacobian(function(alpha) g(alpha)[, l], alpha)
#     }
#     return(apply(Gamma.arr, c(1, 3), mean))
#   }
#
#   obj = function(alpha){
#     g.matrix = g(alpha)
#     G.hat = G(g.matrix)
#     value = t(G.hat)%*%G.hat
#     return(ifelse(is.infinite(value) || is.na(value), 10^8, value))
#   }
#
#   if(is.null(init)) init = rep(0, alpha_dim)
#   alpha_sol_path = matrix(init, alpha_dim)
#   conv_err = 10^8
#   t = 1
#
#   while (conv_err > 10^(-8) & t < 1000){
#     opt = optim(alpha_sol_path[,t], obj, method = "L-BFGS-B")
#     alpha_sol_path = cbind(alpha_sol_path, opt$par)
#     g.matrix = g(alpha_sol_path[,t+1]); W.hat = W(g.matrix);
#     obj = function(alpha){
#       g.matrix = g(alpha); G.hat = G(g.matrix);
#       value = t(G.hat)%*%W.hat%*%G.hat
#       return(ifelse(is.infinite(value) || is.na(value), 10^8, value))
#     }
#     conv_err = max(abs(alpha_sol_path[,t+1]-alpha_sol_path[,t]))
#     t = t + 1
#   }
#
#   alpha.hat = alpha_sol_path[, t]
#   fitted_values = model(x, y, alpha.hat)
#   Gamma.hat = Gamma(alpha.hat)
#   g.matrix = g(alpha.hat)
#   W.hat = W(g.matrix)
#   S = var(g.matrix)
#   # cov.hat = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)/N
#   cov.hat = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)%*%t(Gamma.hat)%*%W.hat%*%S%*%W.hat%*%Gamma.hat%*%solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)/n
#   se = sqrt(diag(cov.hat))
#
#   results = list(coefficients = alpha.hat,
#                  fitted.values = fitted_values,
#                  sol.path = alpha_sol_path,
#                  cov.hat = cov.hat,
#                  se = se,
#                  lower = alpha.hat-qnorm(0.975)*se,
#                  upper = alpha.hat+qnorm(0.975)*se,
#                  g.matrix = g.matrix,
#                  K = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)%*%t(Gamma.hat)%*%W.hat,
#                  model = model,
#                  model_x_names = model_x_names,
#                  h_x = cbind(d, h_x2))
#
#   return(results)
# }

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
# pi.fit.list[[1]]$theta.hat
# pi.fit.list[[1]]$se
# auxilliary = list(cbind(as.factor(u1)), cbind(u2), y)
# est.res1 = ChuangChao2023(pi.fit.list, auxilliary, family, ortho = TRUE, true.pi)
# est1 =  unlist(est.res1[1:4])
# est1
# Sys.time() - start
#
# source("C:\\Users\\stat-pc\\Desktop\\NTHU_Research\\EBMR_Methods_Reproducibility_Materials\\EBMRalgorithm\\R\\Methods.r")
# n = 1000
# formula.list = list(r ~ o(y) + u1 + u2, r ~ o(y) + u1, r ~ o(y) + u2)
# h_x.list = list(
#   c("u1", "u2", "z1", "z2"), c("u1", "z1", "z2"), c("u2", "z1", "z2")
# )
# inv_link = function(eta) 1/(1+exp(eta))
#
# start = Sys.time()
# ps_fit.list = list()
# for(j in 1:J){
#   ps_fit.list[[j]] = WangShaoKim2014(formula.list[[j]], data[h_x.list[[j]]], inv_link, data)
# }
# ps_fit.list[[1]]$coefficients
# ps_fit.list[[1]]$se
# h_x = data[c("u1", "u2")]
# est.res1 = EBMR_IPW(ps_fit.list, h_x, true_ps = true.pi)
# est1 =  unlist(est.res1[1:4])
# est1
# Sys.time() - start
