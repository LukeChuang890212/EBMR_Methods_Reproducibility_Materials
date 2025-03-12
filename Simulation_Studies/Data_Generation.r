setting1.A1 = function(n){
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

setting1.B1 = function(n){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)
  
  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)
  
  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2){
    propensity = 1/(1+exp(0.1-0.1*y-0.1*u1+0.1*u2))*exp(n^(-1/2)*(y))
    # propensity[propensity > 1] = 0.95
    return(propensity)
  }
  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)
  mean(y)
  
  propensity = response.prob(y, u1, u2)
  # propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)
  mean(r)
  
  mean(y[r == 1]); mean(y[r == 0]);
  
  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting1.A2 = function(n){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)
  
  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)
  
  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2)  1/(1+exp(-1.1+0.1*y+0.1*u1+0.1*u2))
  
  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)
  mean(y)
  
  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)
  
  mean(y[r == 1]); mean(y[r == 0]);
  
  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting1.B2 = function(n){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)
  
  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)
  
  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2){
    propensity = 1/(1+exp(-1+0.1*y+0.1*u1+0.1*u2))*exp(n^(-1/2)*(y))
    # propensity[propensity > 1] = 0.95
    return(propensity)
  }
  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)
  mean(y)
  
  propensity = response.prob(y, u1, u2)
  # propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)
  mean(r)
  
  mean(y[r == 1]); mean(y[r == 0]);
  
  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting2.A1 = function(n, response.rate){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)
  
  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)
  
  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2)  1/(1+exp(-0.1+0.5*y-0.5*u1-0.1*u2))
  
  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))
  mean(y)
  
  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)
  
  mean(y[r == 1]); mean(y[r == 0]);
  
  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting2.B1 = function(n){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)
  
  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)
  
  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2){
    propensity = 1/(1+exp(0.05+0.5*y-0.5*u1-0.1*u2))*exp(n^(-1/2)*y)
    # propensity[propensity > 1] = 0.95
    return(propensity)
  }
  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))
  mean(y)
  
  propensity = response.prob(y, u1, u2)
  # propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)
  mean(r)
  
  mean(y[r == 1]); mean(y[r == 0]);
  
  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting2.A2 = function(n){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)
  
  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)
  
  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2)  1/(1+exp(-0.8+0.5*y-0.5*u1-0.1*u2))
  
  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))
  mean(y)
  
  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)
  
  mean(y[r == 1]); mean(y[r == 0]);
  
  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}
setting2.B2 = function(n, response.rate){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)
  
  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)
  
  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2){
    propensity = 1/(1+exp(-0.2-0.5*y-0.5*u1-0.1*u2))*exp(n^(-1/2)*(y))
    # propensity[propensity > 1] = 0.95
    return(propensity)
  }
  
  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))
  mean(y)
  
  propensity = response.prob(y, u1, u2)
  # propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)
  mean(r)
  
  mean(y[r == 1]); mean(y[r == 0]);
  
  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting3.A1 = function(n){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)
  
  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)
  
  m = function(z1, z2, u1, u2) 1+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2)  1/(1+exp(-0.2+0.1*y+0.1*u1+0.1*u2))
  
  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)
  mean(y)
  
  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)
  
  mean(y[r == 1]); mean(y[r == 0]);
  
  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting3.B1 = function(n){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)
  
  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)
  
  m = function(z1, z2, u1, u2) 1+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2){
    propensity = 1/(1+exp(-0.15+0.1*y+0.1*u1+0.1*u2))*exp(n^(-1/2)*y)
    # propensity[propensity > 1] = 0.95
    return(propensity)
  }
  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)
  mean(y)
  
  propensity = response.prob(y, u1, u2)
  # propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)
  mean(r)
  
  mean(y[r == 1]); mean(y[r == 0]);
  
  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting3.A2 = function(n){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)
  
  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)
  
  m = function(z1, z2, u1, u2) 1+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2)  1/(1+exp(-1.1+0.1*y+0.1*u1+0.1*u2))
  
  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)
  mean(y)
  
  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)
  
  mean(y[r == 1]); mean(y[r == 0]);
  
  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting3.B2 = function(n){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)
  
  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)
  
  m = function(z1, z2, u1, u2) 1+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2){
    propensity = 1/(1+exp(-0.9+0.1*y+0.1*u1+0.1*u2))*exp(n^(-1/2)*y)
    # propensity[propensity > 1] = 0.95
    return(propensity)
  }
  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)
  mean(y)
  
  propensity = response.prob(y, u1, u2)
  # propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)
  mean(r)
  
  mean(y[r == 1]); mean(y[r == 0]);
  
  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}
