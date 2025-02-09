KimYu2012 = function(gamma, y, u1, u2, z1, z2, r){
  N = length(r)
  h = apply(cbind(u2, z2), 2, sd)*(N^(-1/5))
  m0 = rep(NA, sum(r == 0))
  id = which(r == 0)
  for(i in 1:length(id)){
    z1.value = as.numeric(z1[id[i]])
    z2.value = ifelse(is.null(z2[id[i]]), NA, z2[id[i]])
    u1.value = ifelse(is.null(u1[id[i]]), NA, u1[id[i]])
    u2.value = ifelse(is.null(u2[id[i]]), NA, u2[id[i]])
    
    w0.value = cbind(ifelse(u1 == u1.value, 1, 0),
                     dnorm((u2-u2.value)/h[1])/h[1],
                     ifelse(z1 == z1.value, 1, 0),
                     dnorm((z2-z2.value)/h[2])/h[2])
    w0.value = r*apply(w0.value, 1, prod)*exp(gamma*y)
    w0.value = w0.value/sum(w0.value)
    
    m0[i] = sum(w0.value*y)
  }
  return(mean(c(y[r == 1], m0)))
}