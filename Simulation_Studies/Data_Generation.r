# A: include correct model
# B: all models are misspecifed

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

setting5.A1 = function(n){
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

  v1 = rnorm(n, 2+r, sd = 0.5)
  v2 = rnorm(n, 1-r, sd = 0.5)
  v3 = rnorm(n, 2, sd = 0.5)
  v4 = rnorm(n, 1, sd = 0.5)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, v1 = v1, v2 = v2, v3 = v3, v4 = v4, y = y, r = r)
  return(dat)
}

setting5.B1 = function(n){
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

  v1 = rnorm(n, 2+r, sd = 0.5)
  v2 = rnorm(n, 1-r, sd = 0.5)
  v3 = rnorm(n, 2, sd = 0.5)
  v4 = rnorm(n, 1, sd = 0.5)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, v1 = v1, v2 = v2, v3 = v3, v4 = v4, y = y, r = r)
  return(dat)
}

setting5.A2 = function(n){
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

  v1 = rnorm(n, 2+r, sd = 0.5)
  v2 = rnorm(n, 1-r, sd = 0.5)
  v3 = rnorm(n, 2, sd = 0.5)
  v4 = rnorm(n, 1, sd = 0.5)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, v1 = v1, v2 = v2, v3 = v3, v4 = v4, y = y, r = r)
  return(dat)
}

setting5.B2 = function(n){
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

  v1 = rnorm(n, 2+r, sd = 0.5)
  v2 = rnorm(n, 1-r, sd = 0.5)
  v3 = rnorm(n, 2, sd = 0.5)
  v4 = rnorm(n, 1, sd = 0.5)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, v1 = v1, v2 = v2, v3 = v3, v4 = v4, y = y, r = r)
  return(dat)
}

setting6.A1 = function(n, response.rate){
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

  v1 = rnorm(n, 2+r, sd = 0.5)
  v2 = rnorm(n, 1-r, sd = 0.5)
  v3 = rnorm(n, 2, sd = 0.5)
  v4 = rnorm(n, 1, sd = 0.5)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, v1 = v1, v2 = v2, y = y, r = r)
  return(dat)
}

setting6.B1 = function(n){
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

  v1 = rnorm(n, 2+r, sd = 0.5)
  v2 = rnorm(n, 1-r, sd = 0.5)
  v3 = rnorm(n, 2, sd = 0.5)
  v4 = rnorm(n, 1, sd = 0.5)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, v1 = v1, v2 = v2, y = y, r = r)
  return(dat)
}

setting6.A2 = function(n){
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

  v1 = rnorm(n, 2+r, sd = 0.5)
  v2 = rnorm(n, 1-r, sd = 0.5)
  v3 = rnorm(n, 2, sd = 0.5)
  v4 = rnorm(n, 1, sd = 0.5)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, v1 = v1, v2 = v2, y = y, r = r)
  return(dat)
}

setting6.B2 = function(n, response.rate){
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

  v1 = rnorm(n, 2+r, sd = 0.5)
  v2 = rnorm(n, 1-r, sd = 0.5)
  v3 = rnorm(n, 2, sd = 0.5)
  v4 = rnorm(n, 1, sd = 0.5)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, v1 = v1, v2 = v2, y = y, r = r)
  return(dat)
}

setting7.A1 = function(n){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 3)

  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2) 1/(1+exp(0.2+0.1*y-0.5*u1-0.5*u2))

  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting7.B1 = function(n){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 3)

  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2

  response.prob = function(y, u1, u2){
    propensity = 1/(1+exp(0.2+0.1*y-0.5*u1-0.5*u2))*exp(n^(-1/2)*(-y+u1-u2))
    propensity[propensity > 1] = 0.95
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

setting7.A2 = function(n){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 3)

  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2)  1/(1+exp(-1+0.1*y-0.5*u1-0.5*u2))

  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting7.B2 = function(n){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 3)

  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2){
    propensity = 1/(1+exp(-1+0.1*y-0.5*u1-0.5*u2))*exp(n^(-1/2)*(-y+u1-u2))
    propensity[propensity > 1] = 0.95
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

setting8.A1 = function(n, response.rate){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 3)

  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2)  1/(1+exp(-0.04+0.2*y-0.2*u1-0.2*u2))

  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting8.B1 = function(n){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 3)

  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2){
    propensity = 1/(1+exp(-0.04+0.2*y-0.2*u1-0.2*u2))*exp(n^(-1/2)*(-y+u1-u2))
    propensity[propensity > 1] = 0.95
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

setting8.A2 = function(n){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 3)

  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2)  1/(1+exp(-0.9+0.2*y-0.2*u1-0.2*u2))

  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting8.B2 = function(n, response.rate){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 3)

  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2){
    propensity = 1/(1+exp(-0.9+0.2*y-0.2*u1-0.2*u2))*exp(n^(-1/2)*(-y+u1-u2))
    propensity[propensity > 1] = 0.95
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

setting9.A1 = function(n){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2) 1/(1+exp(0.6+0.2*y-u1-0.5*u2))

  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting9.B1 = function(n){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2

  response.prob = function(y, u1, u2){
    propensity = 1/(1+exp(0.6+0.2*y-u1-0.5*u2))*exp(n^(-1/2)*(-y+u1-u2))
    propensity[propensity > 1] = 0.95
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

setting9.A2 = function(n){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2)  1/(1+exp(-0.5+0.2*y-u1-0.5*u2))

  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting9.B2 = function(n){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2){
    propensity = 1/(1+exp(-0.5+0.2*y-u1-0.5*u2))*exp(n^(-1/2)*(-y+u1-u2))
    propensity[propensity > 1] = 0.95
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

setting10.A1 = function(n, response.rate){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+z1+z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2)  1/(1+exp(0.4+0.4*y-u1-0.5*u2))

  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting10.B1 = function(n){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+z1+z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2){
    propensity =  1/(1+exp(0.4+0.4*y-u1-0.5*u2))*exp(n^(-1/2)*(-y+u1-u2))
    propensity[propensity > 1] = 0.95
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

setting10.A2 = function(n){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+z1+z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2)  1/(1+exp(-0.6+0.4*y-u1-0.5*u2))

  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting10.B2 = function(n, response.rate){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+z1+z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2){
    propensity = 1/(1+exp(-0.6+0.4*y-u1-0.5*u2))*exp(n^(-1/2)*(-y+u1-u2))
    propensity[propensity > 1] = 0.95
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

setting11.A1 = function(n){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2) 1/(1+exp(1.3+0.2*y-2*u1-0.5*u2))

  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting11.B1 = function(n){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2

  response.prob = function(y, u1, u2){
    propensity = 1/(1+exp(1.3+0.2*y-2*u1-0.5*u2))*exp(n^(-1/2)*(-y+u1-u2))
    propensity[propensity > 1] = 0.95
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

setting11.A2 = function(n){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2)  1/(1+exp(0.2+0.2*y-2*u1-0.5*u2))

  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting11.B2 = function(n){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2){
    propensity = 1/(1+exp(0.2+0.2*y-2*u1-0.5*u2))*exp(n^(-1/2)*(-y+u1-u2))
    propensity[propensity > 1] = 0.95
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

setting12.A1 = function(n, response.rate){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+z1+z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2)  1/(1+exp(1.1+0.4*y-2*u1-0.5*u2))

  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting12.B1 = function(n){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+z1+z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2){
    propensity =   1/(1+exp(1.1+0.4*y-2*u1-0.5*u2))*exp(n^(-1/2)*(-y+u1-u2))
    propensity[propensity > 1] = 0.95
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

setting12.A2 = function(n){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+z1+z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2)   1/(1+exp(0.05+0.4*y-2*u1-0.5*u2))

  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting12.B2 = function(n, response.rate){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+z1+z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2){
    propensity = 1/(1+exp(0.05+0.4*y-2*u1-0.5*u2))*exp(n^(-1/2)*(-y+u1-u2))
    propensity[propensity > 1] = 0.95
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

Cho_RM2.A1 = function(n, response.rate){
  x1 = rnorm(n, mean = 1, sd = sqrt(1/3))
  x2 = rnorm(n, mean = 1, sd = sqrt(1/3))
  x3 = rnorm(n, mean = 1, sd = sqrt(1/3))

  m = function(x1, x2, x3) 0.5+0.5*x1+0.5*x2+0.5*x3
  y = rnorm(n, mean = m(x1, x2, x3), sd = sqrt(1/3))
  mean(y)

  response.prob = function(y, x1, x2){
    eta = -0.98 + 0.5*x1 + 0.25*y
    propensity = exp(eta)/(1+exp(eta))
    return(propensity)
  }

  propensity = response.prob(y, x1, x2)
  # propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(x1 = x1, x2 = x2, x12 = x1^2, x22 = x2^2, y = y, r = r)
  return(dat)
}

Cho_RM2.A2 = function(n, response.rate){
  x1 = rnorm(n, mean = 1, sd = sqrt(1/3))
  x2 = rnorm(n, mean = 1, sd = sqrt(1/3))
  # x2 = rbinom(n, size = 1, prob = 0.35)
  x3 = rnorm(n, mean = 1, sd = sqrt(1/3))
  # x4 = rnorm(n, mean = 1, sd = sqrt(1/3))

  m = function(x1, x2, x3, x4) 0.5+x1+0.5*x2+0.5*x3
  y = rnorm(n, mean = m(x1, x2, x3, x4), sd = sqrt(1/3))
  # y = rbinom(n, size = 1, prob = exp(m(x1, x2, x3))/(1+exp(m(x1, x2, x3))))
  mean(y)

  response.prob = function(y, x1, x2){
    eta = -0.114 + 0.5*x1 + 0.25*y
    propensity = exp(eta)/(1+exp(eta))
    return(propensity)
  }

  propensity = response.prob(y, x1, x2)
  # propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4, y = y, r = r)
  return(dat)
}

Cho_RM3.A1 = function(n, response.rate){
  x1 = rnorm(n, mean = 1, sd = sqrt(1/3))
  x2 = rnorm(n, mean = 1, sd = sqrt(1/3))

  m = function(x1, x2) 0.5+x1+0.5*x2
  y = rnorm(n, mean = m(x1, x2), sd = sqrt(1/3))
  mean(y)

  response.prob = function(y, x1, x2){
    eta = 0.02 + 0.5*x2 - 0.25*y
    propensity = exp(eta)/(1+exp(eta))
    return(propensity)
  }

  propensity = response.prob(y, x1, x2)
  # propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(x1 = x1, x2 = x2, y = y, r = r)
  return(dat)
}

Cho_RM3.A2 = function(n, response.rate){
  x1 = rnorm(n, mean = 1, sd = sqrt(1/3))
  x2 = rnorm(n, mean = 1, sd = sqrt(1/3))

  m = function(x1, x2) 0.5+x1+0.5*x2
  y = rnorm(n, mean = m(x1, x2), sd = sqrt(1/3))
  mean(y)

  response.prob = function(y, x1, x2){
    eta = 0.865 + 0.5*x2 - 0.25*y
    propensity = exp(eta)/(1+exp(eta))
    return(propensity)
  }

  propensity = response.prob(y, x1, x2)
  # propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(x1 = x1, x2 = x2, y = y, r = r)
  return(dat)
}

Cho_RM4 = function(n, response.rate){
  x1 = rnorm(n, mean = 1, sd = sqrt(1/3))
  x2 = rnorm(n, mean = 1, sd = sqrt(1/3))

  m = function(x1, x2) 0.5+x1+0.5*x2
  y = rnorm(n, mean = m(x1, x2), sd = sqrt(1/3))
  mean(y)

  response.prob = function(y, x1, x2){
    eta = ifelse(y <= 2, 1, 0)*(0.857  + 0.5*x1 - 0.25*y) + ifelse(y > 2, 1, 0)*(0.865  + 0.5*x2 - 0.25*y)
    propensity = exp(eta)/(1+exp(eta))
    return(propensity)
  }

  propensity = response.prob(y, x1, x2)
  # propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(x1 = x1, x2 = x2, y = y, r = r)
  return(dat)
}
