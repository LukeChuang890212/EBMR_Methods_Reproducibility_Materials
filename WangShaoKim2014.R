library(stringr)             # Load the package
library(numDeriv)
library(R6)

MyLinearModel <- R6Class("MyLinearModel",
                         public = list(
                           formula = NULL,
                           data = NULL,
                           coefficients = NULL,
                           residuals = NULL,
                           fitted.values = NULL,
                           model = NULL,
                           
                           # Constructor: Initialize with formula and data
                           initialize = function(formula, data) {
                             self$formula <- formula
                             self$data <- data
                           },
                           
                           # Fit the model
                           fit = function() {
                             self$model <- lm(self$formula, data = self$data)
                             self$coefficients <- coef(self$model)
                             self$residuals <- residuals(self$model)
                             self$fitted.values <- fitted(self$model)
                           },
                           
                           # Print method to display results like lm()
                           print = function() {
                             cat("Call:\n")
                             print(self$formula)
                             cat("\nCoefficients:\n")
                             print(self$coefficients)
                           },
                           
                           # Summary method to show detailed output
                           summary = function() {
                             summary(self$model)
                           }
                         )
)

parse_formula <- function(formula) {
  # Convert the formula to a string
  formula <- as.character(formula)
  
  # Extract the left-hand side of the formula (everything before the ~)
  lhs <- formula[2]
  
  # Isolate r 
  r <- lhs
  
  # Extract the right-hand side of the formula (everything after the ~)
  rhs <- formula[3]
  
  # Extract the variables inside o() using regular expressions
  y <- str_extract(rhs, "o\\(([^)]+)\\)")   # Variables in o()
  
  # Clean up: Remove 'o()' and split the variables (there's only one in this case)
  y <- gsub("o\\(|\\)", "", y)
  y <- unlist(strsplit(y, split = "\\+"))
  
  # Extract the rest of the variables (excluding o())
  x <- gsub("o\\([^)]+\\)\\s*\\+?", "", rhs)
  
  # Clean up: Split the rest of the variables by '+' and trim whitespace
  x <- unlist(strsplit(x, split = "\\+"))
  x <- trimws(x)
  
  # Return a list with r, y, x
  return(list(r = r, y = y, x = x))
}

# Function to separate continuous and discrete variable 
separate_variable_types <- function(x) {
  
  # Initialize empty lists to store continuous and discrete columns
  x1 <- list()
  x2 <- list()
  
  col_names = colnames(x)
  x1_names = c()
  x2_names = c()
  
  # Iterate through each column in the matrix 'x'
  for (i in 1:ncol(x)) {
    column <- x[, i]
    
    # Check if the column is numeric
    if (is.numeric(column)) {
      # If the column has more than 10 unique values, it's treated as continuous
      if (length(unique(column)) > 10) {
        x2[[colnames(x)[i]]] <- column
        x2_names = c(x2_names, col_names[i])
      } else {
        x1[[colnames(x)[i]]] <- column
        x1_names = c(x1_names, col_names[i])
      }
    } else {
      # For non-numeric columns, treat them as discrete
      x1[[colnames(x)[i]]] <- column
      x1_names = c(x1_names, col_names[i])
    }
  }
  
  # Convert the lists to matrices (if you want matrices, otherwise keep them as lists)
  x1 <- as.data.frame(x1)
  x2 <- as.data.frame(x2)
  
  # Return the results
  return(list(x1 = x1, x2 = x2, x1_names = x1_names, x2_names = x2_names))
}

WangShaoKim2014 = function(formula, h_x, w, data, init = NULL){
  # Parse the formula
  result = parse_formula(formula)
  r = data[result$r]
  y = data[result$y]
  x = data[result$x]
  
  n = nrow(data)
  
  # result = separate_variable_types(x)
  # model_x1 = result$x1
  # model_x2 = result$x2
  # model_x1_names = result$x1_names
  # model_x2_names = result$x2_names
  
  result = separate_variable_types(h_x)
  h_x1 = result$x1
  h_x2 = result$x2
  
  propensity = function(y, x, alpha, L) 1/w(alpha, y, x, L)

  # h_x1 = h[[1]]; h_x2 = h[[2]];
  # if(!is.null(x1)) x1 = as.matrix(x1)
  # if(!is.null(x2)) x2 = as.matrix(x2)
  # 
  # model_x1 = NULL; model_x2 = NULL;
  # if(!is.null(model_x1_names)) model_x1 = as.matrix(dat[model_x1_names])
  # if(!is.null(model_x2_names)) model_x2 = as.matrix(dat[model_x2_names])

  is.mnar = ifelse(ncol(y) == 0, 0, 1)
  alpha_dim = 1 + is.mnar + ncol(x)

  d = NULL
  if(ncol(h_x1) > 0){
    for(j in 1:ncol(h_x1)) h_x1[, j] = as.factor(h_x1[, j])
    d = model.matrix(lm(rep(1, n)~., data =  h_x1))
  }
  
  discrete_dim = ncol(d)
  continuous_dim = ncol(h_x2)
  h_dim = discrete_dim + continuous_dim
  
  r = as.matrix(r)
  y = as.matrix(y)
  x = as.matrix(x)
  
  g = function(alpha){
    g.matrix = matrix(NA, n, h_dim)
    rw = r*w(alpha, y, x, n)
    if(discrete_dim > 0){
      for(l in 1:discrete_dim){g.matrix[, l] = d[, l]*(rw-1)}
    }
    if(continuous_dim > 0){
      for(l in (discrete_dim+1):(discrete_dim+continuous_dim)) g.matrix[, l] = as.matrix(h_x2)[, l-discrete_dim]*(rw-1)
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
  Gamma.hat = Gamma(alpha.hat)
  g.matrix = g(alpha.hat); W.hat = W(g.matrix);
  S = var(g.matrix)
  # cov.hat = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)/N
  cov.hat = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)%*%t(Gamma.hat)%*%W.hat%*%S%*%W.hat%*%Gamma.hat%*%solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)/n
  se = sqrt(diag(cov.hat))

  return(list(sol.path = alpha_sol_path, alpha.hat = alpha.hat, cov.hat = cov.hat,
              se = se, lower = alpha.hat-qnorm(0.975)*se, upper = alpha.hat+qnorm(0.975)*se,
              g.matrix = g.matrix, K = solve(t(Gamma.hat)%*%W.hat%*%Gamma.hat)%*%t(Gamma.hat)%*%W.hat,
              model = propensity, model.x.names = c(model_x1_names, model_x2_names), h_x = t(cbind(d, h_x2))))
}

Chuang2023.1.1 = function(n){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)
  
  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)
  
  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2)  1/(1+exp(-0.2+0.1*y+0.1*u1+0.1*u2))
  
  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)
  mean(y)
  
  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)
  
  mean(y[r == 1]); mean(y[r == 0]);
  
  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

data = Chuang2023.1.1(1000)
y = data$y
u1 = data$u1
u2 = data$u2
z1 = data$z1
z2 = data$z2

J = 3
dat = data

h.list = list(function(u1, u2, z1, z2) list(cbind(as.factor(u1), as.factor(z1)), cbind(u2, z2)),
              function(u1, u2, z1, z2) list(cbind(as.factor(u1), as.factor(z1)), cbind(z2)),
              function(u1, u2, z1, z2) list(cbind(as.factor(z1)), cbind(u2, z2)))

propensity.list = list(list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = c("u1"),
                            model.x2.names =c("u2")),
                       list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = c("u1"),
                            model.x2.names = NULL),
                       list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = NULL,
                            model.x2.names = c("u2")))
pi.fit.list = list()

start = Sys.time()
pi.fit.list[[3]] = Wang2014.1(auxilliary = h.list[[3]](u1, u2, z1, z2),
                              model.y = propensity.list[[3]]$model.y(y),
                              model.x1.names = propensity.list[[3]]$model.x1.names,
                              model.x2.names = propensity.list[[3]]$model.x2.names,
                              w = propensity.list[[3]]$w, w.prime = propensity.list[[3]]$w.prime)
pi.fit.list[[3]]$theta.hat
pi.fit.list[[3]]$se
Sys.time() - start

formula = r ~ o(y) +  u2
h_x = data[c("u2", "z1", "z2")]

start = Sys.time()
pi.fit.list[[3]] = WangShaoKim2014(formula, h_x, w = propensity.list[[3]]$w, data)
pi.fit.list[[3]]$alpha.hat
pi.fit.list[[3]]$se
Sys.time() - start
