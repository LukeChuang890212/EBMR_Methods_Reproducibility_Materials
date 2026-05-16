# A: include correct model
# B: all models are misspecifed

setting1.A1 = function(n){
  # Target: ~50% missing rate (response rate ~0.5)
  # Propensity = 1/(1+exp(a0 + 0.2*y - 0.8*u1 + 0.8*u2)), no tilt
  # m(x) = 1 + u1 + u2 + z1 + z2, z1~Bern(0.4), z2~N(0,1), u1~Bern(0.6), u2~N(0,1)
  alpha0 = 0.0803
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 1)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 0, sd = 1)

  m = function(z1, z2, u1, u2) 1+u1+u2+z1+z2
  response.prob = function(y, u1, u2) 1/(1+exp(alpha0+0.2*y-0.8*u1+0.8*u2))

  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)
  
  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  
  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting1.B1 = function(n){
  # Target: ~50% missing rate, with local misspecification via linear tilt
  # Propensity = 1/(1+exp(a0+0.2*y-0.8*u1+0.8*u2)) * exp(n^(-1/2)*(y+u1+u2))
  # m(x) = 1 + u1 + u2 + z1 + z2
  # z1~Bern(0.4), z2~N(0,1), u1~Bern(0.6), u2~N(0,1)
  # Intercept depends on n
  alpha0 = ifelse(n >= 1000, 0.1765, 0.2730)
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 1)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 0, sd = 1)
  
  m = function(z1, z2, u1, u2) 1+u1+u2+z1+z2
  
  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)
  
  propensity = 1/(1+exp(alpha0+0.2*y-0.8*u1+0.8*u2))*exp(n^(-1/2)*(y+u1+u2))
  propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)
  
  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting1.A2 = function(n){
  # Target: ~30% missing rate (response rate ~0.7)
  # Propensity = 1/(1+exp(a0 + 0.2*y - 0.8*u1 + 0.8*u2)), no tilt
  # m(x) = 1 + u1 + u2 + z1 + z2, z1~Bern(0.4), z2~N(0,1), u1~Bern(0.6), u2~N(0,1)
  alpha0 = -0.9646
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 1)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 0, sd = 1)
  
  m = function(z1, z2, u1, u2) 1+u1+u2+z1+z2
  response.prob = function(y, u1, u2) 1/(1+exp(alpha0+0.2*y-0.8*u1+0.8*u2))
  
  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)
  
  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  
  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting1.B2 = function(n){
  # Target: ~30% missing rate, with local misspecification via linear tilt
  # Propensity = 1/(1+exp(a0+0.2*y-0.8*u1+0.8*u2)) * exp(n^(-1/2)*(y+u1+u2))
  # m(x) = 1 + u1 + u2 + z1 + z2
  # z1~Bern(0.4), z2~N(0,1), u1~Bern(0.6), u2~N(0,1)
  # Intercept depends on n
  alpha0 = ifelse(n >= 1000, -0.7843, -0.6167)
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 1)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 0, sd = 1)
  
  m = function(z1, z2, u1, u2) 1+u1+u2+z1+z2
  
  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)
  
  propensity = 1/(1+exp(alpha0+0.2*y-0.8*u1+0.8*u2))*exp(n^(-1/2)*(y+u1+u2))
  propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)
  
  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting2.A1 = function(n){
  # Target: ~50% missing rate (response rate ~0.5), correctly specified PS model
  # Binary y with logistic link, Propensity = 1/(1+exp(a0 + 0.2*y - 0.8*u1 + 0.8*u2)), no tilt
  # m(x) = 1 + u1 + u2 + z1 + z2, z1~Bern(0.4), z2~N(0,1), u1~Bern(0.6), u2~N(0,1)
  alpha0 = 0.3193
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 1)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 0, sd = 1)
  
  m = function(z1, z2, u1, u2) 1+u1+u2+z1+z2
  # response.prob = function(y, u1, u2) 1/(1+exp(alpha0+0.15*y-0.30*u1+0.80*u2))
  response.prob = function(y, u1, u2) 1/(1+exp(alpha0+0.2*y-0.8*u1+0.8*u2))
  
  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))
  
  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  
  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting2.B1 = function(n){
  # Target: ~50% missing rate, with local misspecification via 10x linear tilt
  # Propensity = 1/(1+exp(a0+0.2*y-0.8*u1+0.8*u2)) * exp(10*n^(-1/2)*(y+u1+u2))
  # m(x) = 1 + u1 + u2 + z1 + z2
  # z1~Bern(0.4), z2~N(0,1), u1~Bern(0.6), u2~N(0,1)
  # Intercept depends on n
  alpha0 = ifelse(n >= 1000, 0.8794, 1.4141)
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 1)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 0, sd = 1)
  
  m = function(z1, z2, u1, u2) 1+u1+u2+z1+z2
  
  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))
  
  propensity = 1/(1+exp(alpha0+0.2*y-0.8*u1+0.8*u2))*exp(10*n^(-1/2)*(y+u1+u2))
  # propensity = 1/(1+exp(alpha0+y-0.5*u1+0.5*u2))*exp(10*n^(-1/2)*(y+u1+u2))
  propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)
  
  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting2.A2 = function(n){
  # Target: ~30% missing rate (response rate ~0.7), correctly specified PS model
  # Binary y with logistic link, Propensity = 1/(1+exp(a0 + 0.2*y - 0.8*u1 + 0.8*u2)), no tilt
  # m(x) = 1 + u1 + u2 + z1 + z2, z1~Bern(0.4), z2~N(0,1), u1~Bern(0.6), u2~N(0,1)
  alpha0 = -0.6734
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 1)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 0, sd = 1)
  
  m = function(z1, z2, u1, u2) 1+u1+u2+z1+z2
  # response.prob = function(y, u1, u2) 1/(1+exp(alpha0+0.15*y-0.30*u1+0.80*u2))
  response.prob = function(y, u1, u2) 1/(1+exp(alpha0+0.2*y-0.8*u1+0.8*u2))
  
  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))
  
  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))
  
  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting2.B2 = function(n){
  # Target: ~30% missing rate, with local misspecification via 10x linear tilt
  # Propensity = 1/(1+exp(a0+0.2*y-0.8*u1+0.8*u2)) * exp(10*n^(-1/2)*(y+u1+u2))
  # m(x) = 1 + u1 + u2 + z1 + z2
  # z1~Bern(0.4), z2~N(0,1), u1~Bern(0.6), u2~N(0,1)
  # Intercept depends on n
  alpha0 = ifelse(n >= 1000, 0.2012, 0.7290)
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 1)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 0, sd = 1)
  
  m = function(z1, z2, u1, u2) 1+u1+u2+z1+z2
  
  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))
  
  # propensity = 1/(1+exp(alpha0+0.15*y-0.30*u1+0.80*u2))*exp(10*n^(-1/2)*(y+u1+u2))
  propensity = 1/(1+exp(alpha0+0.2*y-0.8*u1+0.8*u2))*exp(10*n^(-1/2)*(y+u1+u2))
  propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)
  
  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting3.A1 = function(n){
  # Target: ~50% missing rate (response rate ~0.5)
  # Propensity = 1/(1+exp(a0 + 0.2*y - 0.4*u1 + 1.0*u2)), no tilt
  # m(x) = 1 + u1 + u2 + z1 + z2, z1~Bern(0.4), u1~Bern(0.6), u2~N(0,1)
  alpha0 = -0.1581
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 2)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 0, sd = 1)

  m = function(z1, z2, u1, u2) 1+u1+u2+z1+z2
  response.prob = function(y, u1, u2) 1/(1+exp(alpha0+0.2*y-0.4*u1+1.0*u2))

  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting3.B1 = function(n){
  # Target: ~50% response rate, with local misspecification via linear tilt
  # Propensity = 1/(1+exp(a0+0.2*y-0.4*u1+1.0*u2)) * exp(n^(-1/2)*(y+u1+u2))
  # m(x) = 1 + u1 + u2 + z1 + z2
  # z1~Bern(0.4), u1~Bern(0.6), u2~N(0,1)
  # Intercept depends on n
  alpha0 = ifelse(n >= 1000, -0.0774, 0.0070)
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 2)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 0, sd = 1)

  m = function(z1, z2, u1, u2) 1+u1+u2+z1+z2

  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)

  propensity = 1/(1+exp(alpha0+0.2*y-0.4*u1+1.0*u2))*exp(n^(-1/2)*(y+u1+u2))
  propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting3.A2 = function(n){
  # Target: ~30% missing rate (response rate ~0.7)
  # Propensity = 1/(1+exp(a0 + 0.2*y - 0.4*u1 + 1.0*u2)), no tilt
  # m(x) = 1 + u1 + u2 + z1 + z2, z1~Bern(0.4), u1~Bern(0.6), u2~N(0,1)
  alpha0 = -1.2697
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 2)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 0, sd = 1)

  m = function(z1, z2, u1, u2) 1+u1+u2+z1+z2
  response.prob = function(y, u1, u2) 1/(1+exp(alpha0+0.2*y-0.4*u1+1.0*u2))

  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting3.B2 = function(n){
  # Target: ~70% response rate, with local misspecification via linear tilt
  # Propensity = 1/(1+exp(a0+0.2*y-0.4*u1+1.0*u2)) * exp(n^(-1/2)*(y+u1+u2))
  # m(x) = 1 + u1 + u2 + z1 + z2
  # z1~Bern(0.4), u1~Bern(0.6), u2~N(0,1)
  # Intercept depends on n
  alpha0 = ifelse(n >= 1000, -1.0998, -0.9435)
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 2)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 0, sd = 1)

  m = function(z1, z2, u1, u2) 1+u1+u2+z1+z2

  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)

  propensity = 1/(1+exp(alpha0+0.2*y-0.4*u1+1.0*u2))*exp(n^(-1/2)*(y+u1+u2))
  propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting4.A1 = function(n){
  # Target: ~50% missing rate (response rate ~0.5), correctly specified PS model
  # Binary y with logistic link, Propensity = 1/(1+exp(a0 + y - 0.5*u1 + 0.5*u2)), no tilt
  # m(x) = 1 + u1 + u2 + z1 + z2, z1~Bern(0.4), z2~N(0,1), u1~Bern(0.6), u2~N(0,1)
  alpha0 = -0.5173
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 2)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 0, sd = 1)

  m = function(z1, z2, u1, u2) 1+u1+u2+z1+z2
  # response.prob = function(y, u1, u2) 1/(1+exp(alpha0+0.15*y-0.30*u1+0.80*u2))
  response.prob = function(y, u1, u2) 1/(1+exp(alpha0+0.15*y-0.30*u1+0.80*u2))

  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting4.B1 = function(n){
  # Target: ~50% missing rate, with local misspecification via 10x linear tilt
  # Propensity = 1/(1+exp(a0+y-0.5*u1+0.5*u2)) * exp(10*n^(-1/2)*(y+u1+u2))
  # m(x) = 1 + u1 + u2 + z1 + z2
  # z1~Bern(0.4), z2~N(0,1), u1~Bern(0.6), u2~N(0,1)
  # Intercept depends on n
  alpha0 = ifelse(n >= 1000, 0.0153, 0.5347)
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 2)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 0, sd = 1)

  m = function(z1, z2, u1, u2) 1+u1+u2+z1+z2

  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))

  propensity = 1/(1+exp(alpha0+0.15*y-0.30*u1+0.80*u2))*exp(10*n^(-1/2)*(y+u1+u2))
  # propensity = 1/(1+exp(alpha0+y-0.5*u1+0.5*u2))*exp(10*n^(-1/2)*(y+u1+u2))
  propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting4.A2 = function(n){
  # Target: ~30% missing rate (response rate ~0.7), correctly specified PS model
  # Binary y with logistic link, Propensity = 1/(1+exp(a0 + y - 0.5*u1 + 0.5*u2)), no tilt
  # m(x) = 1 + u1 + u2 + z1 + z2, z1~Bern(0.4), z2~N(0,1), u1~Bern(0.6), u2~N(0,1)
  alpha0 = -1.4592
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 2)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 0, sd = 1)

  m = function(z1, z2, u1, u2) 1+u1+u2+z1+z2
  # response.prob = function(y, u1, u2) 1/(1+exp(alpha0+0.15*y-0.30*u1+0.80*u2))
  response.prob = function(y, u1, u2) 1/(1+exp(alpha0+0.15*y-0.30*u1+0.80*u2))

  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting4.B2 = function(n){
  # Target: ~30% missing rate, with local misspecification via 10x linear tilt
  # Propensity = 1/(1+exp(a0+y-0.5*u1+0.5*u2)) * exp(10*n^(-1/2)*(y+u1+u2))
  # m(x) = 1 + u1 + u2 + z1 + z2
  # z1~Bern(0.4), z2~N(0,1), u1~Bern(0.6), u2~N(0,1)
  # Intercept depends on n
  alpha0 = ifelse(n >= 1000, -0.6297, -0.0120)
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 2)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 0, sd = 1)

  m = function(z1, z2, u1, u2) 1+u1+u2+z1+z2

  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))

  # propensity = 1/(1+exp(alpha0+0.15*y-0.30*u1+0.80*u2))*exp(10*n^(-1/2)*(y+u1+u2))
  propensity = 1/(1+exp(alpha0+0.15*y-0.30*u1+0.80*u2))*exp(10*n^(-1/2)*(y+u1+u2))
  propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting8.A1 = function(n, response.rate){
  z1 = rbinom(n, size = 1, prob = 0.3)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 3)

  m = function(z1, z2, u1, u2) 0.2+0.5*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2)  1/(1+exp(-0.04+0.2*y+0.2*u1))

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
    propensity = 1/(1+exp(-0.04+0.2*y+0.2*u1))*exp(n^(-1/2)*(-y+u1-u2))
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
  response.prob = function(y, u1, u2)  1/(1+exp(-0.9+0.2*y+0.2*u1))

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
    propensity = 1/(1+exp(-0.9+0.2*y+0.2*u1))*exp(n^(-1/2)*(-y+u1-u2))
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
  z1 = rnorm(n, mean = 0, sd = 1.5)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.3*z1+0.5*z2+0.5*u1+0.5*u2
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
  z1 = rnorm(n, mean = 0, sd = 1.5)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.3*z1+0.5*z2+0.5*u1+0.5*u2

  response.prob = function(y, u1, u2){
    propensity = 1/(1+exp(1.3+0.2*y-2*u1-0.5*u2))*exp(n^(-1/2)*(-0.5*y+u1-u2))
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
  z1 = rnorm(n, mean = 0, sd = 1.5)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.3*z1+0.5*z2+0.5*u1+0.5*u2
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
  z1 = rnorm(n, mean = 0, sd = 1.5)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.3*z1+0.5*z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2){
    propensity = 1/(1+exp(0.2+0.2*y-2*u1-0.5*u2))*exp(n^(-1/2)*(-0.5*y+u1-u2))
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
  z1 = rnorm(n, mean = 0, sd = 1.5)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.5*z1+z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2)  1/(1+exp(1.2+0.4*y-2*u1-0.5*u2))

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
  z1 = rnorm(n, mean = 0, sd = 1.5)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.5*z1+z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2){
    propensity =   1/(1+exp(1.2+0.4*y-2*u1-0.5*u2))*exp(n^(-1/2)*(-0.5*y+u1-u2))
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
  z1 = rnorm(n, mean = 0, sd = 1.5)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.5*z1+z2+0.5*u1+0.5*u2
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
  z1 = rnorm(n, mean = 0, sd = 1.5)
  z2 = rnorm(n, mean = 0, sd = 1)

  u1 = rbinom(n, size = 1, prob = 0.7)
  u2 = rnorm(n, mean = 0, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.5*z1+z2+0.5*u1+0.5*u2
  response.prob = function(y, u1, u2){
    propensity = 1/(1+exp(0.05+0.4*y-2*u1-0.5*u2))*exp(n^(-1/2)*(-0.5*y+u1-u2))
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

setting13.A1 = function(n){
  # Target: ~50% missing rate (response rate ~0.5)
  # alpha = (0.1161, 0.2, -0.4, -0.4)
  # m(x) = 0.2 + 0.3*u1 + 0.3*u2 + 0.6*z1 + 0.6*z2
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 2)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 1, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.3*u1+0.3*u2+0.6*z1+0.6*z2
  response.prob = function(y, u1, u2) 1/(1+exp(0.1161+0.2*y-0.4*u1-0.4*u2))

  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting13.B1 = function(n){
  # Target: ~50% response rate, with local misspecification via exp tilt
  # Propensity = 1/(1+exp(a0+0.2*y-0.4*u1-0.4*u2)) * exp(n^(-1/2)*(y+u1+u2))
  # Intercept depends on n: n=2000 -> 0.1764, n=500 -> 0.2396
  alpha0 = ifelse(n >= 1000, 0.1764, 0.2396)
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 2)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 1, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.3*u1+0.3*u2+0.6*z1+0.6*z2

  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)

  propensity = 1/(1+exp(alpha0+0.2*y-0.4*u1-0.4*u2))*exp(n^(-1/2)*(y+u1+u2))
  propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting13.A2 = function(n){
  # Target: ~30% missing rate (response rate ~0.7)
  # alpha = (-0.7780, 0.2, -0.4, -0.4)
  # m(x) = 0.2 + 0.3*u1 + 0.3*u2 + 0.6*z1 + 0.6*z2

  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 2)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 1, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.3*u1+0.3*u2+0.6*z1+0.6*z2
  response.prob = function(y, u1, u2)  1/(1+exp(-0.7780+0.2*y-0.4*u1-0.4*u2))

  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting13.B2 = function(n){
  # Target: ~70% response rate, with local misspecification via exp tilt
  # Propensity = 1/(1+exp(a0+0.2*y-0.4*u1-0.4*u2)) * exp(n^(-1/2)*(y+u1+u2))
  # Intercept depends on n: n=2000 -> -0.6807, n=500 -> -0.5863
  alpha0 = ifelse(n >= 1000, -0.6807, -0.5863)
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 2)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 1, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.3*u1+0.3*u2+0.6*z1+0.6*z2

  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)

  propensity = 1/(1+exp(alpha0+0.2*y-0.4*u1-0.4*u2))*exp(n^(-1/2)*(y+u1+u2))
  propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting14.A1 = function(n){
  # Target: ~50% missing rate (response rate ~0.5), correctly specified PS model
  # Binary y with logistic link, alpha = (0.1154, 0.2, -0.4, -0.4)
  # m(x) = 0.2 + 0.8*z1 + 0.8*z2 + 0.4*u1 + 0.4*u2
  alpha0 = ifelse(n >= 1000, -0.6653, -0.5688)
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 2)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 1, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.3*u1+0.3*u2+0.6*z1+0.6*z2
  response.prob = function(y, u1, u2)  1/(1+exp(0.1154+0.2*y-0.4*u1-0.4*u2))

  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting14.B1 = function(n){
  # Target: ~50% response rate, with local misspecification via exp tilt
  # Propensity = 1/(1+exp(a0+0.2*y-0.4*u1-0.4*u2)) * exp(n^(-1/2)*(y+u1+u2))
  # Intercept depends on n: n=2000 -> 0.1844, n=500 -> 0.2517
  alpha0 = ifelse(n >= 1000, -0.6653, -0.5688)
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 2)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 1, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.3*u1+0.3*u2+0.6*z1+0.6*z2

  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))

  propensity = 1/(1+exp(alpha0+0.2*y-0.4*u1-0.4*u2))*exp(n^(-1/2)*(y+u1+u2))
  propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting14.A2 = function(n){
  # Target: ~30% missing rate (response rate ~0.7), correctly specified PS model
  # Binary y with logistic link, alpha = (-0.7696, 0.2, -0.4, -0.4)
  # m(x) = 0.2 + 0.8*z1 + 0.8*z2 + 0.4*u1 + 0.4*u2
  alpha0 = ifelse(n >= 1000, -0.6653, -0.5688)
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 2)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 1, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.3*u1+0.3*u2+0.6*z1+0.6*z2
  response.prob = function(y, u1, u2)  1/(1+exp(-0.7696+0.2*y-0.4*u1-0.4*u2))

  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting14.B2 = function(n){
  # Target: ~70% response rate, with local misspecification via exp tilt
  # Propensity = 1/(1+exp(a0+0.2*y-0.4*u1-0.4*u2)) * exp(n^(-1/2)*(y+u1+u2))
  # Intercept depends on n: n=2000 -> -0.6653, n=500 -> -0.5688
  alpha0 = ifelse(n >= 1000, -0.6653, -0.5688)
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 2)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 1, sd = 2)
  # z1 = rbinom(n, size = 1, prob = 0.4)
  # z2 = rnorm(n, mean = 0, sd = 2)
  # u1 = rbinom(n, size = 1, prob = 0.6)
  # u2 = rnorm(n, mean = 0, sd = 1)

  m = function(z1, z2, u1, u2) 0.2+0.3*u1+0.3*u2+0.6*z1+0.6*z2

  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))

  propensity = 1/(1+exp(alpha0+0.2*y-0.4*u1-0.4*u2))*exp(n^(-1/2)*(y+u1+u2))
  propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

#------------------------------------------------------------------------------#
# Setting 15: Binary Y with nonlinear m(x)
# m(x) = 0.1 + 0.1*u1 + 0.1*u2 + 0.1*z1 + 0.1*z2 + 0.05*u2^2 + 0.05*z2^2
# B settings: exp tilt = exp(n^(-1/2)*(y+u1+u2))
#------------------------------------------------------------------------------#

setting15.A1 = function(n){
  # Target: ~50% missing rate (response rate ~0.5), correctly specified PS model
  # Binary y with logistic link, alpha = (0.0917, 0.2, -0.4, -0.4)
  # m(x) = 0.2 + 0.2*u1 + 0.2*u2 + 0.4*z1 + 0.4*z2 + 0.4*u2^2 + 0.4*z2^2
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 1)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 3, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.2*u1+0.2*u2+0.4*z1+0.4*z2+0.4*u2^2+0.4*z2^2
  response.prob = function(y, u1, u2)  1/(1+exp(0.0917+0.2*y-0.4*u1-0.4*u2))

  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting15.B1 = function(n){
  # Target: ~50% response rate, with local misspecification via exp tilt
  # Propensity = 1/(1+exp(a0+0.6*y+0.6*u1+0.6*u2)) * exp(10*n^(-1/2)*(y+u1+u2))
  # Intercept depends on n: n=2000 -> -0.3480, n=500 -> 0.0973
  alpha0 = ifelse(n >= 1000, -0.3480, 0.0973)
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 1)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 3, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.2*u1+0.2*u2+0.4*z1+0.4*z2+0.4*u2^2+0.4*z2^2

  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))

  propensity = 1/(1+exp(alpha0+0.6*y+0.6*u1+0.6*u2))*exp(10*n^(-1/2)*(y+u1+u2))
  propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting15.A2 = function(n){
  # Target: ~30% missing rate (response rate ~0.7), correctly specified PS model
  # Binary y with logistic link, alpha = (-0.7946, 0.2, -0.4, -0.4)
  # m(x) = 0.2 + 0.2*u1 + 0.2*u2 + 0.4*z1 + 0.4*z2 + 0.4*u2^2 + 0.4*z2^2
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 1)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 3, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.2*u1+0.2*u2+0.4*z1+0.4*z2+0.4*u2^2+0.4*z2^2
  response.prob = function(y, u1, u2)  1/(1+exp(-0.7946+0.2*y-0.4*u1-0.4*u2))

  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))

  r = rbinom(n, size = 1, prob = response.prob(y, u1, u2))

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting15.B2 = function(n){
  # Target: ~70% response rate, with local misspecification via exp tilt
  # Propensity = 1/(1+exp(a0+0.6*y+0.6*u1+0.6*u2)) * exp(10*n^(-1/2)*(y+u1+u2))
  # Intercept depends on n: n=2000 -> -1.0052, n=500 -> -0.4191
  alpha0 = ifelse(n >= 1000, -1.0052, -0.4191)
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 1)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 3, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.2*u1+0.2*u2+0.4*z1+0.4*z2+0.4*u2^2+0.4*z2^2

  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))

  propensity = 1/(1+exp(alpha0+0.6*y+0.6*u1+0.6*u2))*exp(10*n^(-1/2)*(y+u1+u2))
  propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

# Setting 16: Same as setting13 but with exp tilt = exp(10*n^(-1/2)*(y+u1+u2))
# Continuous Y, m(x) = 0.2 + 0.6*z1 + 0.6*z2 + 0.3*u1 + 0.3*u2

setting16.B1 = function(n){
  # Target: ~50% missing rate (response rate ~0.5), with local misspecification
  # Perturbed by exp(10*n^(-1/2)*(y+u1+u2))
  # m(x) = 0.2 + 0.6*z1 + 0.6*z2 + 0.3*u1 + 0.3*u2
  # Intercept depends on n: n=2000 -> 0.7340, n=500 -> 1.2000
  alpha0 = ifelse(n >= 1000, 0.7340, 1.2000)
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 2)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 3, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.6*z1+0.6*z2+0.3*u1+0.3*u2

  response.prob = function(y, z1, z2, u1, u2){
    propensity = 1/(1+exp(alpha0+0.2*y-0.4*u1-0.4*u2))*exp(10*n^(-1/2)*(y+u1+u2))
    propensity[propensity > 1] = 0.95
    return(propensity)
  }
  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)

  propensity = response.prob(y, z1, z2, u1, u2)
  propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting16.B2 = function(n){
  # Target: ~30% missing rate (response rate ~0.7), with local misspecification
  # Perturbed by exp(10*n^(-1/2)*(y+u1+u2))
  # m(x) = 0.2 + 0.6*z1 + 0.6*z2 + 0.3*u1 + 0.3*u2
  # Intercept depends on n: n=2000 -> -0.1460, n=500 -> 0.0560
  alpha0 = ifelse(n >= 1000, -0.1460, 0.0560)
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 2)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 3, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.6*z1+0.6*z2+0.3*u1+0.3*u2
  response.prob = function(y, z1, z2, u1, u2){
    propensity = 1/(1+exp(alpha0+0.2*y-0.4*u1-0.4*u2))*exp(10*n^(-1/2)*(y+u1+u2))
    propensity[propensity > 1] = 0.95
    return(propensity)
  }
  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)

  propensity = response.prob(y, z1, z2, u1, u2)
  propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

# Setting 17: Same as setting15 but with exp tilt = exp(10*n^(-1/2)*(y+u1+u2))
# Binary Y, m(x) = 0.2 + 0.6*z1 + 0.6*z2 + 0.3*u1 + 0.3*u2

setting17.B1 = function(n){
  # Target: ~50% missing rate (response rate ~0.5), with local misspecification
  # Perturbed by exp(10*n^(-1/2)*(y+u1+u2))
  # m(x) = 0.2 + 0.6*z1 + 0.6*z2 + 0.3*u1 + 0.3*u2
  # Intercept depends on n: n=2000 -> 0.7370, n=500 -> 1.1980
  alpha0 = ifelse(n >= 1000, 0.7370, 1.1980)
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 2)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 3, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.6*z1+0.6*z2+0.3*u1+0.3*u2
  response.prob = function(y, z1, z2, u1, u2){
    propensity = 1/(1+exp(alpha0+0.2*y-0.4*u1-0.4*u2))*exp(10*n^(-1/2)*(y+u1+u2))
    propensity[propensity > 1] = 0.95
    return(propensity)
  }
  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))

  propensity = response.prob(y, z1, z2, u1, u2)
  propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting17.B2 = function(n){
  # Target: ~30% missing rate (response rate ~0.7), with local misspecification
  # Perturbed by exp(10*n^(-1/2)*(y+u1+u2))
  # m(x) = 0.2 + 0.6*z1 + 0.6*z2 + 0.3*u1 + 0.3*u2
  # Intercept depends on n: n=2000 -> -0.0730, n=500 -> 0.2730
  alpha0 = ifelse(n >= 1000, -0.0730, 0.2730)
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 2)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 3, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.6*z1+0.6*z2+0.3*u1+0.3*u2
  response.prob = function(y, z1, z2, u1, u2){
    propensity = 1/(1+exp(alpha0+0.2*y-0.4*u1-0.4*u2))*exp(10*n^(-1/2)*(y+u1+u2))
    propensity[propensity > 1] = 0.95
    return(propensity)
  }

  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))

  propensity = response.prob(y, z1, z2, u1, u2)
  propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

# Setting 18: Same as setting16 but with exp tilt = exp(n^(-1/2)*(u1+u2+z1+z2))
# Continuous Y, m(x) = 0.2 + 0.6*z1 + 0.6*z2 + 0.3*u1 + 0.3*u2

setting18.B1 = function(n){
  # Target: ~50% missing rate (response rate ~0.5), with local misspecification
  # Base alpha = (0.1900, 0.2, -0.4, -0.4), perturbed by exp(n^(-1/2)*(u1+u2+z1+z2))
  # m(x) = 0.2 + 0.6*z1 + 0.6*z2 + 0.3*u1 + 0.3*u2
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 2)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 3, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.6*z1+0.6*z2+0.3*u1+0.3*u2

  response.prob = function(y, z1, z2, u1, u2){
    propensity = 1/(1+exp(0.1900+0.2*y-0.4*u1-0.4*u2))*exp(n^(-1/2)*(u1+u2+z1+z2))
    propensity[propensity > 1] = 0.95
    return(propensity)
  }
  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)

  propensity = response.prob(y, z1, z2, u1, u2)
  propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting18.B2 = function(n){
  # Target: ~30% missing rate (response rate ~0.7), with local misspecification
  # Base alpha = (-0.6600, 0.2, -0.4, -0.4), perturbed by exp(n^(-1/2)*(u1+u2+z1+z2))
  # m(x) = 0.2 + 0.6*z1 + 0.6*z2 + 0.3*u1 + 0.3*u2
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 2)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 3, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.6*z1+0.6*z2+0.3*u1+0.3*u2
  response.prob = function(y, z1, z2, u1, u2){
    propensity = 1/(1+exp(-0.6600+0.2*y-0.4*u1-0.4*u2))*exp(n^(-1/2)*(u1+u2+z1+z2))
    propensity[propensity > 1] = 0.95
    return(propensity)
  }
  y = rnorm(n, mean = m(z1, z2, u1, u2), sd = 1)

  propensity = response.prob(y, z1, z2, u1, u2)
  propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

# Setting 19: Same as setting17 but with exp tilt = exp(n^(-1/2)*(u1+u2+z1+z2))
# Binary Y, m(x) = 0.2 + 0.6*z1 + 0.6*z2 + 0.3*u1 + 0.3*u2

setting19.B1 = function(n){
  # Target: ~50% missing rate (response rate ~0.5), with local misspecification
  # Base alpha = (0.2100, 0.2, -0.4, -0.4), perturbed by exp(n^(-1/2)*(u1+u2+z1+z2))
  # m(x) = 0.2 + 0.6*z1 + 0.6*z2 + 0.3*u1 + 0.3*u2
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 2)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 3, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.6*z1+0.6*z2+0.3*u1+0.3*u2
  response.prob = function(y, z1, z2, u1, u2){
    propensity = 1/(1+exp(0.2100+0.2*y-0.4*u1-0.4*u2))*exp(n^(-1/2)*(u1+u2+z1+z2))
    propensity[propensity > 1] = 0.95
    return(propensity)
  }
  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))

  propensity = response.prob(y, z1, z2, u1, u2)
  propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

setting19.B2 = function(n){
  # Target: ~30% missing rate (response rate ~0.7), with local misspecification
  # Base alpha = (-0.6400, 0.2, -0.4, -0.4), perturbed by exp(n^(-1/2)*(u1+u2+z1+z2))
  # m(x) = 0.2 + 0.6*z1 + 0.6*z2 + 0.3*u1 + 0.3*u2
  z1 = rbinom(n, size = 1, prob = 0.4)
  z2 = rnorm(n, mean = 0, sd = 2)
  u1 = rbinom(n, size = 1, prob = 0.6)
  u2 = rnorm(n, mean = 3, sd = 2)

  m = function(z1, z2, u1, u2) 0.2+0.6*z1+0.6*z2+0.3*u1+0.3*u2
  response.prob = function(y, z1, z2, u1, u2){
    propensity = 1/(1+exp(-0.6400+0.2*y-0.4*u1-0.4*u2))*exp(n^(-1/2)*(u1+u2+z1+z2))
    propensity[propensity > 1] = 0.95
    return(propensity)
  }

  eta = m(z1, z2, u1, u2)
  y = rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))

  propensity = response.prob(y, z1, z2, u1, u2)
  propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)

  dat = data.frame(z1 = z1, z2 = z2, u1 = u1, u2 = u2, y = y, r = r)
  return(dat)
}

Cho_M1.A1 = function(n, response.rate){
  x1 = rnorm(n, mean = 1, sd = sqrt(1/3))
  x2 = rnorm(n, mean = 1, sd = sqrt(1/3))
  x3 = rnorm(n, mean = 1, sd = sqrt(1/3))

  m = function(x1, x2, x3) 0.5+x1+0.5*x2
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

Cho_M1.A2 = function(n, response.rate){
  x1 = rnorm(n, mean = 1, sd = sqrt(1/3))
  x2 = rnorm(n, mean = 1, sd = sqrt(1/3))
  # x2 = rbinom(n, size = 1, prob = 0.35)
  x3 = rnorm(n, mean = 1, sd = sqrt(1/3))
  x4 = rnorm(n, mean = 1, sd = sqrt(1/3))

  m = function(x1, x2, x3, x4) 0.5+x1+0.5*x2
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

  dat = data.frame(x1 = x1, x2 = x2, x3 = x1^2, x4 = x2^2, y = y, r = r)
  return(dat)
}

Cho_M2.A1 = function(n, response.rate){
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

Cho_M2.A2 = function(n, response.rate){
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

Cho_M1_gamma005.A1 = function(n, response.rate){
  x1 = rnorm(n, mean = 0, sd = sqrt(1/3))
  x2 = rnorm(n, mean = 0, sd = sqrt(1/3))
  x3 = rnorm(n, mean = 0, sd = 1)

  m = function(x1, x2, x3) 0.5+0.5*x1+x2+0.8*x3+0.05*(x1^2-1/3)
  y = rnorm(n, mean = m(x1, x2, x3), sd = sqrt(1/3))
  mean(y)

  response.prob = function(y, x1, x2){
    eta = 0.1244 - 0.5*x1 - 0.25*y
    propensity = exp(eta)/(1+exp(eta))
    return(propensity)
  }

  propensity = response.prob(y, x1, x2)
  # propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(x1 = x1, x2 = x2, x3 = x3, y = y, r = r)
  return(dat)
}

Cho_M1_gamma005.A2 = function(n, response.rate){
  x1 = rnorm(n, mean = 0, sd = sqrt(1/3))
  x2 = rnorm(n, mean = 0, sd = sqrt(1/3))
  x3 = rnorm(n, mean = 0, sd = 1)

  m = function(x1, x2, x3, x4) 0.5+0.5*x1+x2+0.8*x3+0.05*(x1^2-1/3)
  y = rnorm(n, mean = m(x1, x2, x3, x4), sd = sqrt(1/3))
  # y = rbinom(n, size = 1, prob = exp(m(x1, x2, x3))/(1+exp(m(x1, x2, x3))))
  mean(y)

  response.prob = function(y, x1, x2){
    eta = 1.0098 - 0.5*x1 - 0.25*y
    propensity = exp(eta)/(1+exp(eta))
    return(propensity)
  }

  propensity = response.prob(y, x1, x2)
  # propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x2^2, y = y, r = r)
  return(dat)
}

Cho_M2_gamma005.A1 = function(n, response.rate){
  x1 = rnorm(n, mean = 0, sd = sqrt(1/3))
  x2 = rnorm(n, mean = 0, sd = sqrt(1/3))
  x3 = rnorm(n, mean = 0, sd = 1)

  m = function(x1, x2) 0.5+0.5*x1+x2+0.8*x3+0.05*(x1^2-1/3)
  y = rnorm(n, mean = m(x1, x2), sd = sqrt(1/3))
  mean(y)

  response.prob = function(y, x1, x2){
    eta = -0.125 + 0.5*x2 + 0.25*y
    propensity = exp(eta)/(1+exp(eta))
    return(propensity)
  }

  propensity = response.prob(y, x1, x2)
  # propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(x1 = x1, x2 = x2, x3 = x3, y = y, r = r)
  return(dat)
}

Cho_M2_gamma005.A2 = function(n, response.rate){
  x1 = rnorm(n, mean = 0, sd = sqrt(1/3))
  x2 = rnorm(n, mean = 0, sd = sqrt(1/3))
  x3 = rnorm(n, mean = 0, sd = 1)

  m = function(x1, x2) 0.5+0.5*x1+x2+0.8*x3+0.05*(x1^2-1/3)
  y = rnorm(n, mean = m(x1, x2), sd = sqrt(1/3))
  mean(y)

  response.prob = function(y, x1, x2){
    eta = 0.7679 + 0.5*x2 + 0.25*y
    propensity = exp(eta)/(1+exp(eta))
    return(propensity)
  }

  propensity = response.prob(y, x1, x2)
  # propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(x1 = x1, x2 = x2, x3 = x3, y = y, r = r)
  return(dat)
}

Cho_M1_gamma000.A1 = function(n, response.rate){
  x1 = rnorm(n, mean = 0, sd = sqrt(1/3))
  x2 = rnorm(n, mean = 0, sd = sqrt(1/3))
  x3 = rnorm(n, mean = 0, sd = 1)

  m = function(x1, x2, x3) 0.5+0.5*x1+x2+0.8*x3
  y = rnorm(n, mean = m(x1, x2, x3), sd = sqrt(1/3))
  mean(y)

  response.prob = function(y, x1, x2){
    eta = 0.1246 - 0.5*x1 - 0.25*y
    propensity = exp(eta)/(1+exp(eta))
    return(propensity)
  }

  propensity = response.prob(y, x1, x2)
  # propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(x1 = x1, x2 = x2, x3 = x3, y = y, r = r)
  return(dat)
}

Cho_M1_gamma000.A2 = function(n, response.rate){
  x1 = rnorm(n, mean = 0, sd = sqrt(1/3))
  x2 = rnorm(n, mean = 0, sd = sqrt(1/3))
  x3 = rnorm(n, mean = 0, sd = 1)

  m = function(x1, x2, x3, x4) 0.5+0.5*x1+x2+0.8*x3
  y = rnorm(n, mean = m(x1, x2, x3, x4), sd = sqrt(1/3))
  # y = rbinom(n, size = 1, prob = exp(m(x1, x2, x3))/(1+exp(m(x1, x2, x3))))
  mean(y)

  response.prob = function(y, x1, x2){
    eta = 1.0099 - 0.5*x1 - 0.25*y
    propensity = exp(eta)/(1+exp(eta))
    return(propensity)
  }

  propensity = response.prob(y, x1, x2)
  # propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(x1 = x1, x2 = x2, x3 = x3, y = y, r = r)
  return(dat)
}

Cho_M2_gamma000.A1 = function(n, response.rate){
  x1 = rnorm(n, mean = 0, sd = sqrt(1/3))
  x2 = rnorm(n, mean = 0, sd = sqrt(1/3))
  x3 = rnorm(n, mean = 0, sd = 1)

  m = function(x1, x2) 0.5+0.5*x1+x2+0.8*x3
  y = rnorm(n, mean = m(x1, x2), sd = sqrt(1/3))
  mean(y)

  response.prob = function(y, x1, x2){
    eta = -0.125 + 0.5*x2 + 0.25*y
    propensity = exp(eta)/(1+exp(eta))
    return(propensity)
  }

  propensity = response.prob(y, x1, x2)
  # propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(x1 = x1, x2 = x2, x3 = x3, y = y, r = r)
  return(dat)
}

Cho_M2_gamma000.A2 = function(n, response.rate){
  x1 = rnorm(n, mean = 0, sd = sqrt(1/3))
  x2 = rnorm(n, mean = 0, sd = sqrt(1/3))
  x3 = rnorm(n, mean = 0, sd = 1)

  m = function(x1, x2) 0.5+0.5*x1+x2+0.8*x3
  y = rnorm(n, mean = m(x1, x2), sd = sqrt(1/3))
  mean(y)

  response.prob = function(y, x1, x2){
    eta = 0.7679 + 0.5*x2 + 0.25*y
    propensity = exp(eta)/(1+exp(eta))
    return(propensity)
  }

  propensity = response.prob(y, x1, x2)
  # propensity[propensity > 1] = 0.95
  r = rbinom(n, size = 1, prob = propensity)
  mean(r)

  mean(y[r == 1]); mean(y[r == 0]);

  dat = data.frame(x1 = x1, x2 = x2, x3 = x3, y = y, r = r)
  return(dat)
}

# Cho_M1_gamma005.A1 = function(n, response.rate){
#   x1 = rnorm(n, mean = 0, sd = sqrt(1/3))
#   x2 = rnorm(n, mean = 0, sd = sqrt(1/3))
#   
#   m = function(x1, x2, x3) 0.5+0.5*x1+x2+0.05*(x1^2-1/3)
#   y = rnorm(n, mean = m(x1, x2, x3), sd = sqrt(1/3))
#   mean(y)
#   
#   response.prob = function(y, x1, x2){
#     eta = 0.1244 - 0.5*x1 - 0.25*y
#     propensity = exp(eta)/(1+exp(eta))
#     return(propensity)
#   }
#   
#   propensity = response.prob(y, x1, x2)
#   # propensity[propensity > 1] = 0.95
#   r = rbinom(n, size = 1, prob = propensity)
#   mean(r)
#   
#   mean(y[r == 1]); mean(y[r == 0]);
#   
#   dat = data.frame(x1 = x1, x2 = x2, y = y, r = r)
#   return(dat)
# }
# 
# Cho_M1_gamma005.A2 = function(n, response.rate){
#   x1 = rnorm(n, mean = 0, sd = sqrt(1/3))
#   x2 = rnorm(n, mean = 0, sd = sqrt(1/3))
#   
#   m = function(x1, x2, x3, x4) 0.5+0.5*x1+x2+0.05*(x1^2-1/3)
#   y = rnorm(n, mean = m(x1, x2, x3, x4), sd = sqrt(1/3))
#   # y = rbinom(n, size = 1, prob = exp(m(x1, x2, x3))/(1+exp(m(x1, x2, x3))))
#   mean(y)
#   
#   response.prob = function(y, x1, x2){
#     eta = 1.0098 - 0.5*x1 - 0.25*y
#     propensity = exp(eta)/(1+exp(eta))
#     return(propensity)
#   }
#   
#   propensity = response.prob(y, x1, x2)
#   # propensity[propensity > 1] = 0.95
#   r = rbinom(n, size = 1, prob = propensity)
#   mean(r)
#   
#   mean(y[r == 1]); mean(y[r == 0]);
#   
#   dat = data.frame(x1 = x1, x2 = x2, y = y, r = r)
#   return(dat)
# }
# 
# Cho_M2_gamma005.A1 = function(n, response.rate){
#   x1 = rnorm(n, mean = 0, sd = sqrt(1/3))
#   x2 = rnorm(n, mean = 0, sd = sqrt(1/3))
#   
#   m = function(x1, x2) 0.5+0.5*x1+x2+0.05*(x1^2-1/3)
#   y = rnorm(n, mean = m(x1, x2), sd = sqrt(1/3))
#   mean(y)
#   
#   response.prob = function(y, x1, x2){
#     eta = -0.125 + 0.5*x2 + 0.25*y
#     propensity = exp(eta)/(1+exp(eta))
#     return(propensity)
#   }
#   
#   propensity = response.prob(y, x1, x2)
#   # propensity[propensity > 1] = 0.95
#   r = rbinom(n, size = 1, prob = propensity)
#   mean(r)
#   
#   mean(y[r == 1]); mean(y[r == 0]);
#   
#   dat = data.frame(x1 = x1, x2 = x2, y = y, r = r)
#   return(dat)
# }
# 
# Cho_M2_gamma005.A2 = function(n, response.rate){
#   x1 = rnorm(n, mean = 0, sd = sqrt(1/3))
#   x2 = rnorm(n, mean = 0, sd = sqrt(1/3))
#   
#   m = function(x1, x2) 0.5+0.5*x1+x2+0.05*(x1^2-1/3)
#   y = rnorm(n, mean = m(x1, x2), sd = sqrt(1/3))
#   mean(y)
#   
#   response.prob = function(y, x1, x2){
#     eta = 0.7679 + 0.5*x2 + 0.25*y
#     propensity = exp(eta)/(1+exp(eta))
#     return(propensity)
#   }
#   
#   propensity = response.prob(y, x1, x2)
#   # propensity[propensity > 1] = 0.95
#   r = rbinom(n, size = 1, prob = propensity)
#   mean(r)
#   
#   mean(y[r == 1]); mean(y[r == 0]);
#   
#   dat = data.frame(x1 = x1, x2 = x2, y = y, r = r)
#   return(dat)
# }
# 
# Cho_M1_gamma000.A1 = function(n, response.rate){
#   x1 = rnorm(n, mean = 0, sd = sqrt(1/3))
#   x2 = rnorm(n, mean = 0, sd = sqrt(1/3))
#   
#   m = function(x1, x2, x3) 0.5+0.5*x1+x2
#   y = rnorm(n, mean = m(x1, x2, x3), sd = sqrt(1/3))
#   mean(y)
#   
#   response.prob = function(y, x1, x2){
#     eta = 0.1246 - 0.5*x1 - 0.25*y
#     propensity = exp(eta)/(1+exp(eta))
#     return(propensity)
#   }
#   
#   propensity = response.prob(y, x1, x2)
#   # propensity[propensity > 1] = 0.95
#   r = rbinom(n, size = 1, prob = propensity)
#   mean(r)
#   
#   mean(y[r == 1]); mean(y[r == 0]);
#   
#   dat = data.frame(x1 = x1, x2 = x2, y = y, r = r)
#   return(dat)
# }
# 
# Cho_M1_gamma000.A2 = function(n, response.rate){
#   x1 = rnorm(n, mean = 0, sd = sqrt(1/3))
#   x2 = rnorm(n, mean = 0, sd = sqrt(1/3))
#   
#   m = function(x1, x2, x3, x4) 0.5+0.5*x1+x2
#   y = rnorm(n, mean = m(x1, x2, x3, x4), sd = sqrt(1/3))
#   # y = rbinom(n, size = 1, prob = exp(m(x1, x2, x3))/(1+exp(m(x1, x2, x3))))
#   mean(y)
#   
#   response.prob = function(y, x1, x2){
#     eta = 1.0099 - 0.5*x1 - 0.25*y
#     propensity = exp(eta)/(1+exp(eta))
#     return(propensity)
#   }
#   
#   propensity = response.prob(y, x1, x2)
#   # propensity[propensity > 1] = 0.95
#   r = rbinom(n, size = 1, prob = propensity)
#   mean(r)
#   
#   mean(y[r == 1]); mean(y[r == 0]);
#   
#   dat = data.frame(x1 = x1, x2 = x2, y = y, r = r)
#   return(dat)
# }
# 
# Cho_M2_gamma000.A1 = function(n, response.rate){
#   x1 = rnorm(n, mean = 0, sd = sqrt(1/3))
#   x2 = rnorm(n, mean = 0, sd = sqrt(1/3))
#   
#   m = function(x1, x2) 0.5+0.5*x1+x2
#   y = rnorm(n, mean = m(x1, x2), sd = sqrt(1/3))
#   mean(y)
#   
#   response.prob = function(y, x1, x2){
#     eta = -0.125 + 0.5*x2 + 0.25*y
#     propensity = exp(eta)/(1+exp(eta))
#     return(propensity)
#   }
#   
#   propensity = response.prob(y, x1, x2)
#   # propensity[propensity > 1] = 0.95
#   r = rbinom(n, size = 1, prob = propensity)
#   mean(r)
#   
#   mean(y[r == 1]); mean(y[r == 0]);
#   
#   dat = data.frame(x1 = x1, x2 = x2, y = y, r = r)
#   return(dat)
# }
# 
# Cho_M2_gamma000.A2 = function(n, response.rate){
#   x1 = rnorm(n, mean = 0, sd = sqrt(1/3))
#   x2 = rnorm(n, mean = 0, sd = sqrt(1/3))
#   
#   m = function(x1, x2) 0.5+0.5*x1+x2
#   y = rnorm(n, mean = m(x1, x2), sd = sqrt(1/3))
#   mean(y)
#   
#   response.prob = function(y, x1, x2){
#     eta = 0.7679 + 0.5*x2 + 0.25*y
#     propensity = exp(eta)/(1+exp(eta))
#     return(propensity)
#   }
#   
#   propensity = response.prob(y, x1, x2)
#   # propensity[propensity > 1] = 0.95
#   r = rbinom(n, size = 1, prob = propensity)
#   mean(r)
#   
#   mean(y[r == 1]); mean(y[r == 0]);
#   
#   dat = data.frame(x1 = x1, x2 = x2, y = y, r = r)
#   return(dat)
# }
