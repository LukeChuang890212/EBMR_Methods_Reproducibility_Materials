Chuang2023.3.1 = function(n, response.rate){
  z1 = rmultinom(n, size = 1, prob = c(0.3, 0.3, 0.4))
  z1 = apply(z1, 2, function(x) which(x == 1))
  z2 = rnorm(n, mean = 0, sd = 2)

  u1 = rmultinom(n, size = 1, prob = c(0.5, 0.25, 0.25))
  u1 = apply(u1, 2, function(x) which(x == 1))
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2)  1/(1+exp(-0.3+0.05*y+0.1*u1+0.1*u2))

  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

Chuang2023.3.1.mild2 = function(n, response.rate){
  z1 = rmultinom(n, size = 1, prob = c(0.3, 0.3, 0.4))
  z1 = apply(z1, 2, function(x) which(x == 1))
  z2 = rnorm(n, mean = 0, sd = 2)

  u1 = rmultinom(n, size = 1, prob = c(0.5, 0.25, 0.25))
  u1 = apply(u1, 2, function(x) which(x == 1))
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2)  1/(1+exp(-0.1+0.05*y+0.1*u1+0.1*u2))*exp(n^(-1/2)*y)

  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

Chuang2023.3.2 = function(n, response.rate){
  z1 = rmultinom(n, size = 1, prob = c(0.3, 0.3, 0.4))
  z1 = apply(z1, 2, function(x) which(x == 1))
  z2 = rnorm(n, mean = 0, sd = 2)

  u1 = rmultinom(n, size = 1, prob = c(0.5, 0.25, 0.25))
  u1 = apply(u1, 2, function(x) which(x == 1))
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2)  1/(1+exp(-1.2+0.05*y+0.1*u1+0.1*u2))

  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

Chuang2023.3.2.mild2 = function(n, response.rate){
  z1 = rmultinom(n, size = 1, prob = c(0.3, 0.3, 0.4))
  z1 = apply(z1, 2, function(x) which(x == 1))
  z2 = rnorm(n, mean = 0, sd = 2)

  u1 = rmultinom(n, size = 1, prob = c(0.5, 0.25, 0.25))
  u1 = apply(u1, 2, function(x) which(x == 1))
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2)  1/(1+exp(-0.95+0.05*y+0.1*u1+0.1*u2))*exp(n^(-1/2)*y)

  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

Chuang2023.1.3 = function(n, response.rate){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 1+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2)  1/(1+exp(-0.2+0.04*y+0.2*u1+0.2*u2))

  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

Chuang2023.1.3.mild = function(n, response.rate){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 1+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2){
    propensity = 1/(1+exp(-0.1-0.1*y+0.1*u1+0.1*u2))*exp(n^(-1/2)*(-y-u1-u2))
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

Chuang2023.1.3.mild2 = function(n, response.rate){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 1+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2){
    propensity = 1/(1+exp(0.1-0.1*y+0.1*u1+0.1*u2))*exp(n^(-1/2)*y)
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

Chuang2023.1.4 = function(n, response.rate){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 1+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2)  1/(1+exp(-0.5-0.2*y-0.5*u1-0.5*u2))

  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

Chuang2023.1.4.mild = function(n, response.rate){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 1+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2){
    propensity = 1/(1+exp(-0.7-0.5*y-0.5*u1-0.5*u2))*exp(n^(-1/2)*(-y-u1-u2))
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

Chuang2023.1.4.mild2 = function(n, response.rate){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 1+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2){
    propensity = 1/(1+exp(-0.5-0.1*y-0.1*u1+0.1*u2))*exp(n^(-1/2)*(y))
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

Chuang2023.1.5 = function(n, response.rate){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 1+2*z1+2*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2)  1/(1+exp(0.2-0.05*y-0.2*u1+0.2*u2))

  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

Chuang2023.1.5.mild = function(n, response.rate){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 1+2*z1+2*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2){
    propensity = 1/(1+exp(0.1-0.05*y-0.2*u1+0.2*u2))*exp(n^(-1/2)*(-y-u1-u2))
    print(sum(propensity > 1))
    propensity[propensity > 1] = 0.95
    return(propensity)
  }

  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

Chuang2023.1.6 = function(n, response.rate){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 1+2*z1+2*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2)  1/(1+exp(-0.5-0.05*y-0.4*u1+0.2*u2))

  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

Chuang2023.1.6.mild = function(n, response.rate){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 1+2*z1+2*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2){
    propensity = 1/(1+exp(-1.3+0.1*y-0.5*u1-0.5*u2))*exp(n^(-1/2)*(-y-u1-u2))
    print(sum(propensity > 1))
    propensity[propensity > 1] = 0.95
    return(propensity)
  }

  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

Chuang2023.4.1 = function(n, response.rate){
  z1 = rmultinom(n, size = 1, prob = c(0.3, 0.3, 0.4))
  z1 = apply(z1, 2, function(x) which(x == 1))
  z2 = rnorm(n, mean = 0, sd = 2)

  u1 = rmultinom(n, size = 1, prob = c(0.5, 0.25, 0.25))
  u1 = apply(u1, 2, function(x) which(x == 1))
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2)  1/(1+exp(0.4+0.5*y-0.5*u1-0.1*u2))

  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

Chuang2023.4.1.mild2 = function(n, response.rate){
  z1 = rmultinom(n, size = 1, prob = c(0.3, 0.3, 0.4))
  z1 = apply(z1, 2, function(x) which(x == 1))
  z2 = rnorm(n, mean = 0, sd = 2)

  u1 = rmultinom(n, size = 1, prob = c(0.5, 0.25, 0.25))
  u1 = apply(u1, 2, function(x) which(x == 1))
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2)  1/(1+exp(0.5+0.5*y-0.5*u1-0.1*u2))*exp(N^(-1/2)*y)

  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

Chuang2023.4.2 = function(n, response.rate){
  z1 = rmultinom(n, size = 1, prob = c(0.3, 0.3, 0.4))
  z1 = apply(z1, 2, function(x) which(x == 1))
  z2 = rnorm(n, mean = 0, sd = 2)

  u1 = rmultinom(n, size = 1, prob = c(0.5, 0.25, 0.25))
  u1 = apply(u1, 2, function(x) which(x == 1))
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2)  1/(1+exp(-0.4+0.5*y-0.5*u1-0.1*u2))

  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

Chuang2023.4.2.mild2 = function(n, response.rate){
  z1 = rmultinom(n, size = 1, prob = c(0.3, 0.3, 0.4))
  z1 = apply(z1, 2, function(x) which(x == 1))
  z2 = rnorm(n, mean = 0, sd = 2)

  u1 = rmultinom(n, size = 1, prob = c(0.5, 0.25, 0.25))
  u1 = apply(u1, 2, function(x) which(x == 1))
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2)  1/(1+exp(-0.4+0.5*y-0.5*u1-0.1*u2))*exp(n^(-1/2)*y)

  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

Chuang2023.2.3 = function(n, response.rate){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) -0.2-0.5*z1-0.5*z2-0.5*u1-0.5*u2
  response.prob = function(y, u1, u2)  1/(1+exp(-0.1+0.8*y-0.3*u1-0.3*u2))

  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

Chuang2023.2.3.mild = function(n, response.rate){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) -0.2-0.5*z1-0.5*z2-0.5*u1-0.5*u2
  response.prob = function(y, u1, u2){
    propensity = 1/(1+exp(-0.2+0.8*y-0.3*u1-0.3*u2))*exp(n^(-1/2)*(-y-u1-u2))
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

Chuang2023.2.3.mild2 = function(n, response.rate){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) -0.2-0.5*z1-0.5*z2-0.5*u1-0.5*u2
  response.prob = function(y, u1, u2){
    propensity = 1/(1+exp(-0.1+0.8*y-0.3*u1-0.3*u2))*exp(n^(-1/2)*(y))
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

Chuang2023.2.4 = function(n, response.rate){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) -0.2-0.5*z1-0.5*z2-0.5*u1-0.5*u2
  response.prob = function(y, u1, u2)  1/(1+exp(-0.3-0.8*y-0.5*u1-0.3*u2))

  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

Chuang2023.2.4.mild = function(n, response.rate){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) -0.2-0.5*z1-0.5*z2-0.5*u1-0.5*u2
  response.prob = function(y, u1, u2){
    propensity = 1/(1+exp(-0.6-0.8*y-0.5*u1-0.5*u2))*exp(n^(-1/2)*(-y-u1-u2))
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

Chuang2023.2.4.mild2 = function(n, response.rate){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) -0.2-0.5*z1-0.5*z2-0.5*u1-0.5*u2
  response.prob = function(y, u1, u2){
    propensity = 1/(1+exp(-0.2-0.8*y-0.5*u1-0.1*u2))*exp(n^(-1/2)*(y))
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

Chuang2023.2.5 = function(n, response.rate){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+z1+z2+0.1*u1+0.1*u2
  response.prob = function(y, u1, u2)  1/(1+exp(-0.1+0.1*y+0.1*u1-0.2*u2))

  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

Chuang2023.2.5.mild = function(n, response.rate){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+z1+z2+0.1*u1+0.1*u2
  response.prob = function(y, u1, u2){
    propensity = 1/(1+exp(-0.2+0.1*y+0.1*u1-0.2*u2))*exp(n^(-1/2)*(-y-u1-u2))
    propensity[propensity > 1] = 0.95
    return(propensity)
  }

  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

Chuang2023.2.6 = function(n, response.rate){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+z1+z2+0.1*u1+0.1*u2
  response.prob = function(y, u1, u2)  1/(1+exp(-0.5-0.2*y-0.5*u1-0.2*u2))

  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

Chuang2023.2.6.mild = function(n, response.rate){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 2)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+z1+z2+0.1*u1+0.1*u2
  response.prob = function(y, u1, u2){
    propensity = 1/(1+exp(-0.6-0.2*y-0.6*u1-0.2*u2))*exp(n^(-1/2)*(-y-u1-u2))
    propensity[propensity > 1] = 0.95
    return(propensity)
  }

  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))
  mean(y)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}
