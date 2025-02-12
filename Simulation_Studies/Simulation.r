simulate1 = function(alldata, alpha.true, true.pi.model, family, h.list, propensity.list, N, save.file = NULL){
  source("DataGeneration_SM.r", local = TRUE)
  source("Wang2014_12.r", local = TRUE)
  source("ChuangChao2023_SM.r", local = TRUE)
  # source("Chuang2024_SM.r", local = TRUE)
  # source("KimYu2012.r", local = TRUE)
  # source("Morikawa2021.r", local = TRUE)
  # source("Han2018_2.r", local = TRUE)
  # source("ZhaoMa2022.r", local = TRUE)
  
  true.mu = mean(alldata$y)
  
  cores = detectCores()
  cl = makeCluster(cores - 4) ## number of clusters
  registerDoSNOW(cl)
  
  pb = txtProgressBar(max = n.sim, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts = list(progress = progress)
  parallel_packages = c("pracma", "tidyverse", "Matrix", "MASS", "sandwich", "randomForest")
  
  start = Sys.time()
  sim.res = foreach(i= 1:n.sim, .combine = 'cbind', .options.snow = opts, .packages = parallel_packages) %dopar% {
    tryCatch({
      dat = alldata[((i-1)*N+1):(i*N),]
      # dat = Chuang2023.1.1(n = N)
      
      y = dat$y
      u1 = dat$u1
      u2 = dat$u2
      z1 = dat$z1
      z2 = dat$z2
      r = dat$r
      
      n = sum(r)
      J = length(propensity.list)
      
      true.pi = true.pi.model(y, u1, u2, r, N, alpha.true)
      
      
      
      c(est1, 
        unlist(lapply(pi.fit.list, function(pi.fit) pi.fit$theta.hat)),
        unlist(lapply(pi.fit.list, function(pi.fit) pi.fit$se)), 
        est.res1$w.pi,
        est.res1$nu,
        est.res1$imbalance)    
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  close(pb)
  print(Sys.time()-start)
  
  if(!is.null(save.file)){
    if(file.exists(save.file)){
      sim.res.tmp = sim.res
      source(save.file)
      sim.res = cbind(sim.res, sim.res.tmp)
    }
    
    dump("sim.res", save.file)
  }
  
  cat("\n", "Before removing the outliers", "\n")
  print(rbind(apply(sim.res, 1, mean, na.rm = TRUE), apply(sim.res, 1, sd, na.rm = TRUE))); cat("\n");
  cat("After removing the outliers", "\n")
  print(rbind(apply(sim.res, 1, function(v) if(!all(is.na(v))) mean(rm.extreme(na.omit(v))) else NA),
              apply(sim.res, 1, function(v) if(!all(is.na(v))) sd(rm.extreme(na.omit(v))) else NA)))
  cat("\n", "Number of replicates:", ncol(sim.res), "/", 
      "Number of outliers removed:",
      apply(sim.res, 1, function(v) if(!any(is.na(v))) ncol(sim.res)-length(rm.extreme(na.omit(v))) else 0),
      "\n\n")
  
  for(j in 1){
    cp = rep(NA, 2)
    for(i in 1:2){
      ci = cbind(sim.res[4*(j-1)+i,]-1.96*sim.res[4*(j-1)+2+i,],
                 sim.res[4*(j-1)+i,]+1.96*sim.res[4*(j-1)+2+i,])
      cp[i] = mean(na.omit(apply(ci, 1, function(v) ifelse(v[1] < true.mu & v[2] > true.mu, 1, 0))))
    }
    names(cp) = c("IPW.CP", "IPW.true.CP")
    print(cp)
  }
  
  gc()
}

simulate_all = function(read.file, save.root, vers, alpha.true, N.v, miss.true.pi.model, true.pi.model, h.list, propensity.list){
  dat = readRDS(read.file)
  dat$r[is.na(dat$r)] = 1
  plot.misspecification(dat, alpha.true, miss.true.pi.model, propensity.list,
                        N = 10^4,
                        save.file = paste0(c("ChuangResults_SM_Graphs/ChuangChao2023_", save.root, "_", N.v[1], "_", vers, ".png"), collapse = ""))
  for(N in N.v){
    for(model.num in 1:length(propensity.list)){
      model.sets = combn(length(propensity.list), model.num)
      for(i in 1:ncol(model.sets)){
        model.set = model.sets[, i]
        save.file = paste0(c("ChuangResults_SM/ChuangChao2023_", save.root, "_", model.set, "_", N, "_", vers, ".RData"), collapse = "")
        print(paste("model set:", paste0(model.set, collapse = ""), "/", "N:", N, "/", "file:", save.file))
        simulate1(dat, alpha.true, true.pi.model, "gaussian", h.list[model.set], propensity.list[model.set], N, save.file)
      }
    }
  }
  # for(N in N.v){
  #   for(model.num in 3){
  #     model.sets = combn(length(propensity.list), model.num)
  #     for(i in 1){
  #       model.set = model.sets[, i]
  #       save.file = paste0(c("ChuangResults_SM/ChuangChao2023_", save.root, "_", model.set, "_", N, "_", vers, ".RData"), collapse = "")
  #       print(paste("model set:", paste0(model.set, collapse = ""), "/", "N:", N, "/", "file:", save.file))
  #       simulate1(dat, alpha.true, true.pi.model, "gaussian", propensity.list[model.set], N)
  #     }
  #   }
  # }
}

N.v = c(1000, 300)
true.pi.model = function(y, u1, u2, r, n, alpha.true) 1/(1+exp(cbind(rep(1, n), y, u1, u2)%*%alpha.true))
miss.true.pi.model = function(y, u1, u2, r, n, alpha.true) 1/(1+exp(cbind(rep(1, n), y, u1, u2)%*%alpha.true))*exp(n^(-1/2)*y)

# version: 21, 29, 34, 43
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
                            model.x2.names = NULL))

# version: 30, 31, 32, 38
propensity.list = list(list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = c("u1"),
                            model.x2.names =c("u2")),
                       list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = NULL,
                            model.x2.names = c("u2")),
                       list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = NULL,
                            model.x2.names = NULL))

# version: 45
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


# version: 22, 25
propensity.list = list(list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = c("u1"),
                            model.x2.names =c("u2")),
                       list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = NULL,
                            model.x2.names = c("u2")),
                       list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), x)%*%theta),
                            model.y = function(y) NULL,
                            model.x1.names = NULL,
                            model.x2.names = c("u2")))

# version: 23, 26
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
                       list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), x)%*%theta),
                            model.y = function(y) NULL,
                            model.x1.names = c("u1"),
                            model.x2.names = NULL))

# version: 24
propensity.list = list(list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = c("u1"),
                            model.x2.names = c("u2")),
                       list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) NULL,
                            model.x1.names = c("u1"),
                            model.x2.names = c("u2")),
                       list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), x)%*%theta),
                            model.y = function(y) NULL,
                            model.x1.names = c("u1"),
                            model.x2.names = NULL))

# version: 40
propensity.list = list(list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = c("u1"),
                            model.x2.names = c("u2")),
                       list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = c("z1"),
                            model.x2.names = NULL),
                       list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = NULL,
                            model.x2.names = c("z2")))

# version: 50
propensity.list = list(list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = c("u1"),
                            model.x2.names = c("u2")),
                       list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = c("z1"),
                            model.x2.names = c("u2")),
                       list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = c("u1"),
                            model.x2.names = c("z2")))

# version: 51
h.list = list(function(u1, u2, z1, z2) list(cbind(as.factor(u1), as.factor(z1)), cbind(u2, z2)),
              function(u1, u2, z1, z2) list(cbind(as.factor(z1)), cbind(z2)),
              function(u1, u2, z1, z2) list(cbind(as.factor(u1), as.factor(z1)), cbind(z2)))

propensity.list = list(list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = c("u1"),
                            model.x2.names = c("u2")),
                       list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = c("z1"),
                            model.x2.names = NULL),
                       list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = c("u1"),
                            model.x2.names = NULL))

# version: 52, 62
h.list = list(function(u1, u2, z1, z2) list(cbind(as.factor(u1), as.factor(z1)), cbind(u2, z2)),
              function(u1, u2, z1, z2) list(cbind(as.factor(z1)), cbind(z2)),
              function(u1, u2, z1, z2) list(cbind(as.factor(z1)), cbind(u2, z2)))

propensity.list = list(list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = c("u1"),
                            model.x2.names = c("u2")),
                       list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = NULL,
                            model.x2.names = c("z2")),
                       list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = NULL,
                            model.x2.names = c("u2")))

# version: 55
h.list = list(function(u1, u2, z1, z2) list(cbind(as.factor(u1), as.factor(z1)), cbind(u2, z2)),
              function(u1, u2, z1, z2) list(cbind(as.factor(u1), as.factor(z1)), cbind(u2)),
              function(u1, u2, z1, z2) list(cbind(as.factor(z1)), cbind(u2, z2)))

propensity.list = list(list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = c("u1"),
                            model.x2.names = c("u2")),
                       list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = c("u1"),
                            model.x2.names = c("u2")),
                       list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = c("u2"),
                            model.x2.names = NULL))

#===============================================================================

simulate_all(read.file = "ChuangData_SM/ChuangData_SM1.1_2000.RDS",
             save.root = "ChuangData1.1-nested",
             vers = "62",
             alpha.true = c(-0.2, 0.1, 0.1, 0.1),
             N.v = c(1000, 300), miss.true.pi.model, true.pi.model, h.list, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM1.1.mild2_1000.RDS",
             save.root = "ChuangData1.1.mild2-nested",
             vers = "62",
             alpha.true = c(0.1, -0.1, -0.1, 0.1),
             N.v = c(1000), miss.true.pi.model, true.pi.model, h.list, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM1.1.mild2_300.RDS",
             save.root = "ChuangData1.1.mild2-nested",
             vers = "62",
             alpha.true = c(0.1, -0.1, -0.1, 0.1),
             N.v = c(300), miss.true.pi.model, true.pi.model, h.list, propensity.list)

#===============================================================================

simulate_all(read.file = "ChuangData_SM/ChuangData_SM1.2_2000.RDS",
             save.root = "ChuangData1.2-nested",
             vers = "62",
             alpha.true = c(-1.1, 0.1, 0.1, 0.1),
             N.v = c(1000, 300), miss.true.pi.model, true.pi.model, h.list, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM1.2.mild2_1000.RDS",
             save.root = "ChuangData1.2.mild2-nested",
             vers = "62",
             alpha.true = c(-1, 0.1, 0.1, 0.1),
             N.v = c(1000), miss.true.pi.model, true.pi.model, h.list, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM1.2.mild2_300.RDS",
             save.root = "ChuangData1.2.mild2-nested",
             vers = "62",
             alpha.true = c(-1, 0.1, 0.1, 0.1),
             N.v = c(300), miss.true.pi.model, true.pi.model, h.list, propensity.list)

#===============================================================================

simulate_all(read.file = "ChuangData_SM/ChuangData_SM1.3_2000.RDS",
             save.root = "ChuangData1.3-nested",
             vers = "62",
             alpha.true = c(-0.2, 0.04, 0.2, 0.2),
             N.v = c(1000, 300), miss.true.pi.model, true.pi.model, h.list, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM1.3.mild2_1000.RDS",
             save.root = "ChuangData1.3.mild2-nested",
             vers = "62",
             alpha.true = c(0.1, -0.1, 0.1, 0.1),
             N.v = c(1000), miss.true.pi.model, true.pi.model, h.list, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM1.3.mild2_300.RDS",
             save.root = "ChuangData1.3.mild2-nested",
             vers = "62",
             alpha.true = c(0.1, -0.1, 0.1, 0.1),
             N.v = c(300), miss.true.pi.model, true.pi.model, h.list, propensity.list)

#===============================================================================

simulate_all(read.file = "ChuangData_SM/ChuangData_SM1.4_2000.RDS",
             save.root = "ChuangData1.4-nested",
             vers = "62",
             alpha.true = c(-0.5, -0.2, -0.5, -0.5),
             N.v = c(1000, 300), miss.true.pi.model, true.pi.model, h.list, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM1.4.mild2_1000.RDS",
             save.root = "ChuangData1.4.mild2-nested",
             vers = "62",
             alpha.true = c(-0.5, -0.1, -0.1, 0.1),
             N.v = c(1000), miss.true.pi.model, true.pi.model, h.list, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM1.4.mild2_300.RDS",
             save.root = "ChuangData1.4.mild2-nested",
             vers = "62",
             alpha.true = c(-0.5, -0.1, -0.1, 0.1),
             N.v = c(300), miss.true.pi.model, true.pi.model, h.list, propensity.list)

#===============================================================================

simulate_all(read.file = "ChuangData_SM/ChuangData_SM2.1_2000.RDS",
             save.root = "ChuangData2.1-nested",
             vers = "62",
             alpha.true = c(-0.1, 0.5, -0.5, -0.1),
             N.v = c(1000, 300), miss.true.pi.model, true.pi.model, h.list, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM2.1.mild2_1000.RDS",
             save.root = "ChuangData2.1.mild2-nested",
             vers = "62",
             alpha.true = c(0.05, 0.5, -0.5, -0.1),
             N.v = c(1000), miss.true.pi.model, true.pi.model, h.list, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM2.1.mild2_300.RDS",
             save.root = "ChuangData2.1.mild2-nested",
             vers = "62",
             alpha.true = c(0.05, 0.5, -0.5, -0.1),
             N.v = c(300), miss.true.pi.model, true.pi.model, h.list, propensity.list)

#===============================================================================

simulate_all(read.file = "ChuangData_SM/ChuangData_SM2.2_2000.RDS",
             save.root = "ChuangData2.2-nested",
             vers = "62",
             alpha.true = c(-0.8, 0.5, -0.5, -0.1),
             N.v = c(1000, 300), miss.true.pi.model, true.pi.model, h.list, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM2.2.mild2_1000.RDS",
             save.root = "ChuangData2.2.mild2-nested",
             vers = "62",
             alpha.true = c(-0.2, -0.5, -0.5, -0.1),
             N.v = c(1000), miss.true.pi.model, true.pi.model, h.list, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM2.2.mild2_300.RDS",
             save.root = "ChuangData2.2.mild2-nested",
             vers = "62",
             alpha.true = c(-0.2, -0.5, -0.5, -0.1),
             N.v = c(300), miss.true.pi.model, true.pi.model, h.list, propensity.list)

#===============================================================================

simulate_all(read.file = "ChuangData_SM/ChuangData_SM2.3_2000.RDS",
             save.root = "ChuangData2.3-nested",
             vers = "62",
             alpha.true = c(-0.2, 0.1, 0.1, 0.1),
             N.v = c(1000, 300), miss.true.pi.model, true.pi.model, h.list, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM2.3.mild2_1000.RDS",
             save.root = "ChuangData2.3.mild2-nested",
             vers = "62",
             alpha.true = c(-0.1, 0.8, -0.3, -0.3),
             N.v = c(1000), miss.true.pi.model, true.pi.model, h.list, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM2.3.mild2_300.RDS",
             save.root = "ChuangData2.3.mild2-nested",
             vers = "62",
             alpha.true = c(-0.1, 0.8, -0.3, -0.3),
             N.v = c(300), miss.true.pi.model, true.pi.model, h.list, propensity.list)

#===============================================================================

simulate_all(read.file = "ChuangData_SM/ChuangData_SM2.4_2000.RDS",
             save.root = "ChuangData2.4-nested",
             vers = "62",
             alpha.true = c(-0.3, -0.8, -0.5, -0.3),
             N.v = c(1000, 300), miss.true.pi.model, true.pi.model, h.list, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM2.4.mild2_1000.RDS",
             save.root = "ChuangData2.4.mild2-nested",
             vers = "62",
             alpha.true = c(-0.2, -0.8, -0.5, -0.1),
             N.v = c(1000), miss.true.pi.model, true.pi.model, h.list, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM2.4.mild2_300.RDS",
             save.root = "ChuangData2.4.mild2-nested",
             vers = "62",
             alpha.true = c(-0.2, -0.8, -0.5, -0.1),
             N.v = c(300), miss.true.pi.model, true.pi.model, h.list, propensity.list)

#===============================================================================

propensity.list = list(list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = c("u1"),
                            model.x2.names =c("u2")),
                       list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = c("z1"),
                            model.x2.names = c("u2")),
                       list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = NULL,
                            model.x2.names = NULL))

simulate_all(read.file = "ChuangData_SM/ChuangData_SM1.5_1000.RDS",
             save.root = "ChuangData1.5-nested",
             vers = "31",
             alpha.true = c(0.2, -0.05, -0.2, 0.2),
             N.v = c(1000, 300),
             true.pi.model, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM1.5.mild2_1000.RDS",
             save.root = "ChuangData1.5.mild2-nested",
             vers = "31",
             alpha.true = c(0.1, -0.05, -0.2, 0.2),
             N.v = c(1000),
             true.pi.model, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM1.5.mild2_300.RDS",
             save.root = "ChuangData1.5.mild2-nested",
             vers = "31",
             alpha.true = c(0.1, -0.05, -0.2, 0.2),
             N.v = c(300),
             true.pi.model, propensity.list)

#===============================================================================

propensity.list = list(list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = c("u1"),
                            model.x2.names =c("u2")),
                       list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = c("z1"),
                            model.x2.names = c("u2")),
                       list(w = function(theta, y, x, L) 1+exp(cbind(rep(1, L), y, x)%*%theta),
                            w.prime = function(theta, y, x, L) exp(cbind(rep(1, L), y, x)%*%theta),
                            model.y = function(y) y,
                            model.x1.names = NULL,
                            model.x2.names = NULL))

simulate_all(read.file = "ChuangData_SM/ChuangData_SM1.6_1000.RDS",
             save.root = "ChuangData1.6-nested",
             vers = "31",
             alpha.true = c(-0.5, -0.05, -0.4, 0.2),
             N.v = c(1000, 300),
             true.pi.model, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM1.6.mild2_1000.RDS",
             save.root = "ChuangData1.6.mild2-nested",
             vers = "31",
             alpha.true = c(-1, -0.5, -1, -0.5),
             N.v = c(1000),
             true.pi.model, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM1.6.mild2_300.RDS",
             save.root = "ChuangData1.6.mild2-nested",
             vers = "31",
             alpha.true = c(-1, -0.5, -1, -0.5),
             N.v = c(300),
             true.pi.model, propensity.list)

#===============================================================================

simulate_all(read.file = "ChuangData_SM/ChuangData_SM2.5_1000.RDS",
             save.root = "ChuangData2.5-nested",
             vers = "31",
             alpha.true = c(-0.1, 0.1, 0.1, -0.2),
             N.v = c(1000, 300),
             true.pi.model, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM2.5.mild2_1000.RDS",
             save.root = "ChuangData2.5.mild2-nested",
             vers = "31",
             alpha.true = c(-0.2, 0.1, 0.1, -0.2),
             N.v = c(1000),
             true.pi.model, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM2.5.mild2_300.RDS",
             save.root = "ChuangData2.5.mild2-nested",
             vers = "31",
             alpha.true = c(-0.2, 0.1, 0.1, -0.2),
             N.v = c(300),
             true.pi.model, propensity.list)

#===============================================================================

simulate_all(read.file = "ChuangData_SM/ChuangData_SM2.6_1000.RDS",
             save.root = "ChuangData2.6-nested",
             vers = "31",
             alpha.true = c(-0.5, -0.2, -0.5, -0.2),
             N.v = c(1000, 300),
             true.pi.model, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM2.6.mild2_1000.RDS",
             save.root = "ChuangData2.6.mild2-nested",
             vers = "31",
             alpha.true = c(-0.6, -0.2, -0.6, -0.2),
             N.v = c(1000),
             true.pi.model, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM2.6.mild2_300.RDS",
             save.root = "ChuangData2.6.mild2-nested",
             vers = "31",
             alpha.true = c(-0.6, -0.2, -0.6, -0.2),
             N.v = c(300),
             true.pi.model, propensity.list)

#===============================================================================

simulate_all(read.file = "ChuangData_SM/ChuangData_SM3.1_1000.RDS",
             save.root = "ChuangData3.1-nested",
             vers = "43",
             alpha.true = c(-0.3, 0.05, 0.1, 0.1),
             N.v = c(1000, 300), miss.true.pi.model, true.pi.model, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM3.1.mild2_1000.RDS",
             save.root = "ChuangData3.1.mild2-nested",
             vers = "43",
             alpha.true = c(-0.1, 0.05, 0.1, 0.1),
             N.v = c(1000), miss.true.pi.model, true.pi.model, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM3.1.mild2_300.RDS",
             save.root = "ChuangData3.1.mild2-nested",
             vers = "43",
             alpha.true = c(-0.1, 0.05, 0.1, 0.1),
             N.v = c(300), miss.true.pi.model, true.pi.model, propensity.list)

#===============================================================================

simulate_all(read.file = "ChuangData_SM/ChuangData_SM3.2_1000.RDS",
             save.root = "ChuangData3.2-nested",
             vers = "43",
             alpha.true = c(-1.2, 0.05, 0.1, 0.1),
             N.v = c(1000, 300), miss.true.pi.model, true.pi.model, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM3.2.mild2_1000.RDS",
             save.root = "ChuangData3.2.mild2-nested",
             vers = "43",
             alpha.true = c(-0.95, 0.05, 0.1, 0.1),
             N.v = c(1000), miss.true.pi.model, true.pi.model, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM3.2.mild2_300.RDS",
             save.root = "ChuangData3.2.mild2-nested",
             vers = "43",
             alpha.true = c(-0.95, 0.05, 0.1, 0.1),
             N.v = c(300), miss.true.pi.model, true.pi.model, propensity.list)

#===============================================================================

simulate_all(read.file = "ChuangData_SM/ChuangData_SM4.1_1000.RDS",
             save.root = "ChuangData4.1-nested",
             vers = "43",
             alpha.true = c(0.4, 0.5, -0.5, -0.1),
             N.v = c(1000, 300), miss.true.pi.model, true.pi.model, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM4.1.mild2_1000.RDS",
             save.root = "ChuangData4.1.mild2-nested",
             vers = "43",
             alpha.true = c(0.5, 0.5, -0.5, -0.1),
             N.v = c(1000), miss.true.pi.model, true.pi.model, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM4.1.mild2_300.RDS",
             save.root = "ChuangData4.1.mild2-nested",
             vers = "43",
             alpha.true = c(0.5, 0.5, -0.5, -0.1),
             N.v = c(300), miss.true.pi.model, true.pi.model, propensity.list)

#===============================================================================

simulate_all(read.file = "ChuangData_SM/ChuangData_SM4.2_1000.RDS",
             save.root = "ChuangData4.2-nested",
             vers = "43",
             alpha.true = c(-0.4, 0.5, -0.5, -0.1),
             N.v = c(1000, 300), miss.true.pi.model, true.pi.model, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM4.2.mild2_1000.RDS",
             save.root = "ChuangData4.2.mild2-nested",
             vers = "43",
             alpha.true = c(-1, 0.1, 0.1, 0.1),
             N.v = c(1000), miss.true.pi.model, true.pi.model, propensity.list)

simulate_all(read.file = "ChuangData_SM/ChuangData_SM4.2.mild2_300.RDS",
             save.root = "ChuangData4.2.mild2-nested",
             vers = "43",
             alpha.true = c(-1, 0.1, 0.1, 0.1),
             N.v = c(300), miss.true.pi.model, true.pi.model, propensity.list)

#===============================================================================



