simulate = function(all_data, ps_model.true, alpha.true, ps_specifications, n, replicate_num, save.file = NULL){
  # source("Data_Generation.r", local = TRUE)
  library(EBMRalgorithm)
  library(parallel)
  library(foreach)
  library(doSNOW)

  mu.true = mean(all_data$y)

  cores = detectCores()
  cl = makeCluster(cores - 6) ## number of clusters
  registerDoSNOW(cl)

  pb = txtProgressBar(max = replicate_num, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts = list(progress = progress)
  parallel_packages = c("EBMRalgorithm")

  start = Sys.time()
  sim_result = foreach(i = 1:replicate_num, .combine = 'cbind', .options.snow = opts, .packages = parallel_packages) %dopar% {
    tryCatch({
      dat = all_data[((i-1)*n+1):(i*n), ]

      ebmr = EBMRAlgorithm$new(y_names = "y",
                               ps_specifications = ps_specifications,
                               data = dat)
      result = ebmr$EBMR_IPW(h_x_names = c("u1", "u2"),
                             true_ps = ps_model.true(dat$y, dat$u1, dat$u2, dat$r, alpha.true))
      estimates = unlist(result[1:4])

      c(estimates,
        alpha.hat = unlist(lapply(ebmr$ps_fit.list, function(ps_fit) ps_fit$coefficients)),
        se_alpha.hat = unlist(lapply(ebmr$ps_fit.list, function(ps_fit) ps_fit$se)),
        nu.hat = result$nu.hat,
        w.hat = result$w.hat
      )
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  close(pb)
  print(Sys.time()-start)

  if(!is.null(save.file)){
    if(file.exists(save.file)){
      sim_result_temp = sim_result
      sim_result = readRDS(save.file)
      sim_result = cbind(sim_result, sim_result_temp)
    }
    saveRDS("sim_result", save.file)
  }

  cat("\n", "Before removing the outliers", "\n")
  print(rbind(apply(sim_result, 1, mean, na.rm = TRUE), apply(sim_result, 1, sd, na.rm = TRUE))); cat("\n");
  cat("After removing the outliers", "\n")
  print(rbind(apply(sim_result, 1, function(v) if(!all(is.na(v))) mean(rm.extreme(na.omit(v))) else NA),
              apply(sim_result, 1, function(v) if(!all(is.na(v))) sd(rm.extreme(na.omit(v))) else NA)))
  cat("\n", "Number of replicates:", ncol(sim_result), "/",
      "Number of outliers removed:",
      apply(sim_result, 1, function(v) if(!any(is.na(v))) ncol(sim_result)-length(rm.extreme(na.omit(v))) else 0),
      "\n\n")

  cp = rep(NA, 2)
  for(i in 1:2){
    ci = cbind(sim_result[i,] - 1.96*sim_result[2 + i,],
               sim_result[i,] + 1.96*sim_result[2 + i,])
    cp[i] = mean(na.omit(apply(ci, 1, function(v) ifelse(v[1] < mu.true & v[2] > mu.true, 1, 0))))
  }
  names(cp) = c("ipw_cp", "ipw.true_cp")
  print(cp)

  gc()
}

simulate_all_model_combinations_and_sample_sizes = function(
                                                   setting,
                                                   scenario,
                                                   missing_rate,
                                                   full_ps_specifications,
                                                   n.vector,
                                                   all_data_file,
                                                   alpha.true,
                                                   ps_model.true,
                                                   replicate_num,
                                                   version){
  all_data = readRDS(all_data_file)
  all_data$r[is.na(all_data$r)] = 1
  # plot.misspecification(dat, alpha.true, miss.ps_model.true, propensity.list,
  #                       N = 10^4,
  #                       save.file = paste0(c("ChuangResults_SM_Graphs/ChuangChao2023_",save_root, "_", n.vector[1], "_", version, ".png"), collapse = ""))

  J = length(full_ps_specifications)
  for(n in n.vector){
    for(model_num in 1:J){
      model_combinations = combn(J, model_num)
      for(i in 1:ncol(model_combinations)){
        model_set = model_combinations[, i]
        save.file = paste0(c("Simulation_Results/EBMR_IPW_", setting, "-", missing_rate, "-scenario", scenario, "_", model_set, "_n", n, "_replicate", replicate_num, "_", version, ".RDS"), collapse = "")
        print(paste("model set:", paste0(model_set, collapse = ""), "/", "file:", save.file))
        ps_specifications <- list(
          formula.list = full_ps_specifications$formula.list[model_set],
          h_x_names.list = full_ps_specifications$h_x_names.list[model_set],
          inv_link = full_ps_specifications$inv_link
        )
        simulate(all_data, ps_model.true, alpha.true, ps_specifications, n, replicate_num, save.file)
      }
    }
  }
}

sim_all_settings_with_all_missing_mechanisms = function(scenario,
                                                        full_ps_specifications,
                                                        n.vector.list,
                                                        all_data_file.list,
                                                        alpha.true.list,
                                                        version){
  for(setting in c("setting1", "setting2")){
      for(missing_rate in c("miss50", "miss30")){
        for(i in 1:length(all_data_file.list[[setting]][[missing_rate]])){
        simulate_all_model_combinations_and_sample_sizes(
          setting = setting,
          scenario = scenario,
          missing_rate = missing_rate,
          full_ps_specifications = full_ps_specifications,
          n.vector = n.vector.list[[i]],
          all_data_file = all_data_file.list[[setting]][[missing_rate]][[i]],
          alpha.true = alpha.true.list[[setting]][[missing_rate]][[i]],
          ps_model.true = ps_model.true,
          replicate_num,
          version = version
        )
      }
    }
  }
}

