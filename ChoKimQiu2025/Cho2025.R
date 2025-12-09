Cho2025 = function(dat, or_x_names, s, ps_fit.list, Aux_names, inv_link){
    source("../ChoKimQiu2025/functions_common.R")
    source("../ChoKimQiu2025/functions_sim2.R")
    
    family = ifelse(length(unique(dat$y)) > 2, "gaussian", "binomial")
    
    deltavec = dat$r
    ysvec = dat$y[deltavec == 1]
    
    Xtmp <- as.matrix(cbind(1, dat[deltavec == 1, or_x_names]))
    # nuis_init <- as.vector(solve(t(Xtmp)%*%Xtmp)%*%t(Xtmp)%*%ysvec)
    nuis_init = glm(ysvec~Xtmp-1, family = family)$coef
    
    J = length(ps_fit.list)
    n1 = sum(deltavec)
    n = nrow(dat)
    pismat = matrix(NA, n1, J)
    pimat = matrix(NA, n, J)
    pitil = matrix(NA, n1, J)
    pvec = NULL
    is.mnar = unlist(lapply(ps_fit.list, function(ps_fit) ps_fit$is.mnar))
    
    ##############################
    # Invariant terms
    for(j in 1:J){
      pismat[, j] <- ps_fit.list[[j]]$fitted.values[deltavec == 1]
    }
    if(any(!is.mnar)){
      pimat[, which(!is.mnar)] <- ps_fit.list[[which(!is.mnar)]]$fitted.values
      pitil[, which(!is.mnar)] <- pismat[, which(!is.mnar)] - mean(pimat[, which(!is.mnar)])
    }
    H <- 1
    ##############################
    # Iteration
    if(J > 1){
      nuis_prev <- nuis_init
      is_conv <- T
      count <- 0
      or_x = as.matrix(dat[or_x_names])
      while(TRUE){
        count <- count+1
        for(j in which(is.mnar)){
          ##############################
          ### Step 1. EL-step: Update pvec
          eta <- c(cbind(1, or_x)%*%nuis_prev)
          ## PS
          phihat = ps_fit.list[[j]]$coefficients
          Xmat = as.matrix(dat[ps_fit.list[[j]]$model_x_names])
          q <- length(phihat)
          eta_tmp <- c(cbind(1, eta, Xmat)%*%phihat)
          # eta_x <- c(cbind(1, Xmat)%*%phihat[-q])
          res_int = NULL
          if(family == "gaussian"){
            res_int = integrate(function(uvec){
              value = sapply(uvec,function(u){
                pivec <- inv_link(eta_tmp+phihat[2]*s*u)
                return( mean( (1-pivec)^H*pivec )*dnorm(x=u,mean=0,sd=1) )
              })
              # sapply(uvec,function(u){
              #   pivec <- inv_link(eta_x+phihat[q]*u)
              #   return( mean( (1-pivec)^H*pivec )*dnorm(x=(u-eta)/s,mean=0,sd=1) )
              # })
              return(ifelse(is.finite(value), value, 0))
            },-Inf,Inf)
          }else if(family == "binomial"){
            res_int = list(value = sum(sapply(0:1, function(u){
              pivec <- inv_link(c(cbind(1, Xmat)%*%phihat[-2])+phihat[2]*u)
              return( mean((1-pivec)^H*pivec*exp(eta)^u/(1+exp(eta))))
            })))
          }
          pitil[, j] <- pismat[, j] - res_int$value - sum(1-(1-pismat[, j])^H)/n
        }
        
        ## Combine the results.
        Aux = as.matrix(dat[Aux_names])
        Aux = t(t(Aux[deltavec == 1, ]) - apply(Aux, 2, mean))
        res_el <- EstEL1(pitil, Aux, n-n1)
        pvec <- res_el$pvec
        
        ##############################
        ### Step 2. M-step
        nuis_new = glm(ysvec~Xtmp-1, family = family, weights = pvec)$coef
        # nuis_new <- c(solve(t(Xtmp)%*%diag(pvec)%*%Xtmp)%*%t(Xtmp)%*%diag(pvec)%*%ysvec)
        diff <- sqrt(sum((nuis_new-nuis_prev)^2))
        if(diff < 1e-6 ){
          break
        }else{
          nuis_prev <- nuis_new
        }
        if(count >= 1000){
          is_conv <- F
          break
        }
      }
    
      # for(j in which(is.mnar)){
      #   # Impose MBC assumption
      #   phihat = ps_fit.list[[j]]$coefficients
      #   q <- length(phihat)
      #   Xmat = as.matrix(dat[ps_fit.list[[j]]$model_x_names])
      #   eta <- c(cbind(1, or_x)%*%nuis_prev)
      #   est_adj = NULL
      #   if(family == "gaussian"){
      #     est_adj <- function(adj){
      #       res_int <- integrate(function(uvec){
      #         value = sapply(uvec,function(u){
      #           sum(inv_link(adj+c(cbind(1, Xmat)%*%phihat[-2])+
      #                          phihat[2]*( eta+s*u )))*
      #             dnorm(x=u,mean=0,sd=1)
      #         })
      #         return(ifelse(is.finite(value), value, 0))
      #       },-Inf,Inf)
      #       return(res_int$value-n1)
      #     }
      #   }else if(family == "binomial"){
      #     est_adj <- function(adj){
      #       res_int = list(value = sum(sapply(0:1, function(u){
      #         sum(inv_link(c(cbind(1, Xmat)%*%phihat[-2])+phihat[2]*u)*
      #                       exp(eta)^u/(1+exp(eta)))
      #       })))
      #       return(res_int$value-n1)
      #     }
      #   }
      #   M <- 1
      #   while(TRUE){
      #     if( est_adj(M)*est_adj(-M)<0 ){
      #       break
      #     }else{
      #       M <- M+1
      #     }
      #   }
      #   res_adj <- uniroot(est_adj,interval=c(-M,M))
      #   phihatadj <- phihat
      #   phihatadj[1] <- phihatadj[1] + res_adj$root
      #   pismat[, j] <- inv_link(cbind(1, ysvec, Xmat[deltavec == 1,])%*%phihatadj)
      # }
      # 
      # Aux = as.matrix(dat[Aux_names])
      # Aux = t(t(Aux[deltavec == 1, ]) - apply(Aux, 2, mean))
      # res_el <- EstEL0(pismat, Aux, n-n1)
      # pvec <- res_el$pvec
    }else{
      Aux = as.matrix(dat[Aux_names])
      Aux = t(t(Aux[deltavec == 1, ]) - apply(Aux, 2, mean))
      res_el <- EstEL0(as.matrix(ps_fit.list[[1]]$fitted.values[deltavec == 1]), Aux, n-n1)
      pvec <- res_el$pvec
    }
    
    ##########################################
    muhat <- sum(pvec*ysvec)/sum(pvec)
    # res_quan <- uniroot(Estquan,interval=range(ysvec)+c(-3,3),ysvec=ysvec,pvec=pvec,alpha=alpha_quan)
    # quanhat <- res_quan$root
    # thetahat <- c(muhat,quanhat)
    # cbind(theta0,thetahat)
    return(muhat)
}
