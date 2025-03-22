#------------------------------------------------------------------------------#
# Functions for running the simulation results ----
#------------------------------------------------------------------------------#

simulate = function(all_data, ps_model.true, alpha.true, ps_specifications, n, replicate_num, save_file = NULL){
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
  parallel_packages = c("EBMRalgorithm", "stringr", "Matrix", "numDeriv")

  start = Sys.time()
  sim_result = foreach(i= 1:replicate_num, .combine = 'cbind', .options.snow = opts, .packages = parallel_packages, .export = c("WangShaoKim2014", "EBMR_IPW", "estimate_nu", "ensemble", "parse_formula", "separate_variable_types")) %dopar% {
    tryCatch({
      dat = all_data[((i-1)*n+1):(i*n), ]

      # ebmr = EBMRAlgorithm$new(y_names = "y",
      #                          ps_specifications = ps_specifications,
      #                          data = dat)
      # result = ebmr$EBMR_IPW(h_x_names = c("u1", "u2"),
      #                        true_ps = ps_model.true(dat$y, dat$u1, dat$u2, dat$r, alpha.true))

      W = function(g.matrix){
        return(solve(t(g.matrix)%*%g.matrix/n))
      }

      ps_fit.list = list()
      J = length(ps_specifications$formula.list)
      for(j in 1:J){
        formula = ps_specifications$formula.list[[j]]
        h_x_names = ps_specifications$h_x_names.list[[j]]
        inv_link = ps_specifications$inv_link
        ps_fit.list[[j]] = WangShaoKim2014(formula, h_x_names, inv_link, W, data = dat)
      }

      # j = 1
      # formula = ps_specifications$formula.list[[j]]
      # h_x_names = ps_specifications$h_x_names.list[[j]]
      # inv_link = ps_specifications$inv_link
      # ps_fit.list[[j]] = WangShaoKim2014(formula, h_x_names, inv_link, data = dat)
      # ps_fit.list[[j]]$coefficients
      # ps_fit.list[[j]]$se

      result = EBMR_IPW(h_x_names = c("u1", "u2"), W, dat, se.fit = TRUE,
                        true_ps = ps_model.true(dat$y, dat$u1, dat$u2, dat$r, alpha.true))
      estimates = unlist(result[1:4])

      c(estimates,
        alpha.hat = unlist(lapply(ps_fit.list, function(ps_fit) ps_fit$coefficients)),
        se_alpha.hat = unlist(lapply(ps_fit.list, function(ps_fit) ps_fit$se)),
        nu.hat = result$nu.hat,
        w.hat = result$w.hat
      )
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  close(pb)
  print(Sys.time()-start)

  if(!is.null(save_file)){
    if(file.exists(save_file)){
      sim_result_temp = sim_result
      sim_result = readRDS(save_file)
      sim_result = cbind(sim_result, sim_result_temp)
    }
    saveRDS(sim_result, save_file)
  }

  cat("\n", "Before removing the outliers", "\n")
  print(rbind(apply(sim_result, 1, mean, na.rm = TRUE), apply(sim_result, 1, sd, na.rm = TRUE))); cat("\n");
  # cat("After removing the outliers", "\n")
  # print(rbind(apply(sim_result, 1, function(v) if(!all(is.na(v))) mean(rm.extreme(na.omit(v))) else NA),
  #             apply(sim_result, 1, function(v) if(!all(is.na(v))) sd(rm.extreme(na.omit(v))) else NA)))
  # cat("\n", "Number of replicates:", ncol(sim_result), "/",
  #     "Number of outliers removed:",
  #     apply(sim_result, 1, function(v) if(!any(is.na(v))) ncol(sim_result)-length(rm.extreme(na.omit(v))) else 0),
  #     "\n\n")

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
  #                       save_file = paste0(c("ChuangResults_SM_Graphs/ChuangChao2023_",save_root, "_", n.vector[1], "_", version, ".png"), collapse = ""))

  J = length(full_ps_specifications)
  for(n in n.vector){
    for(model_num in 1:J){
      model_combinations = combn(J, model_num)
      for(i in 1:ncol(model_combinations)){
        model_set = model_combinations[, i]
        save_file = paste0(c("Simulation_Results/EBMR_IPW_", setting, "-", missing_rate, "-scenario", scenario, "_", model_set, "_n", n, "_replicate", replicate_num, "_", version, ".RDS"), collapse = "")
        print(paste("model set:", paste0(model_set, collapse = ""), "/", "file:", save_file))
        ps_specifications <- list(
          formula.list = full_ps_specifications$formula.list[model_set],
          h_x_names.list = full_ps_specifications$h_x_names.list[model_set],
          inv_link = full_ps_specifications$inv_link
        )
        simulate(all_data, ps_model.true, alpha.true, ps_specifications, n, replicate_num, save_file)
      }
    }
  }
}

simulate_all_settings_with_all_missing_rates = function(settings,
                                                        missing_rates,
                                                        scenario,
                                                        full_ps_specifications,
                                                        n.vector.list,
                                                        all_data_file.list,
                                                        alpha.true.list,
                                                        version){
  for(setting in settings){
    for(missing_rate in missing_rates){
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

#------------------------------------------------------------------------------#
# Functions for summarizing the simulation results ----
#------------------------------------------------------------------------------#

# which.extreme = function(v){
#   v = na.omit(v)
#   lower = quantile(v, 0.25) - 1.5*IQR(v)
#   upper = quantile(v, 0.75) + 1.5*IQR(v)
#   print(paste(sum(v < lower | v > upper), length(v)))
#   return(which(v < lower | v > upper))
# }

which.not.extreme = function(v){
  z.score = scale(v)
  extreme_num = sum(abs(z.score) > 3)
  print(paste(extreme_num, length(v)))
  if(extreme_num == 0){
    return(1:length(v))
  }else{
    return((1:length(v))[-order(abs(z.score), decreasing = TRUE)[1:min(extreme_num, extreme_num)]])
  }
}

summarize_results = function(sim_result, pe_index, ese_index, mu.true, is.original){
  # process_sim_replicates = switch(is.original,
  #                                 `TRUE` = function(v) v,
  #                                 `FALSE` = rm.extreme)

  if(!is.original){
    print(ncol(sim_result)-length(which.not.extreme(sim_result[ese_index,])))
    sim_result = sim_result[, which.not.extreme(sim_result[ese_index,])]
  }
  # pe = rep(NA, length(pe_index))
  # esd = rep(NA, length(pe_index))
  # for(i in 1:length(pe_index)){
  #   pe[i] = mean(process_sim_replicates(sim_result[pe_index[i],]))
  #   esd[i] = sd(process_sim_replicates(sim_result[pe_index[i],]))
  # }

  pe = mean(sim_result[pe_index,])
  esd = sd(sim_result[pe_index,])

  bias = pe - mu.true
  mse = bias^2+esd^2
  bias = round(bias, 3)
  mse = round(mse, 3)
  esd = round(esd, 3)

  ese = sim_result[ese_index,]
  cp = rep(NA, length(ese_index))
  for(k in 1:length(ese_index)){
    ci = cbind(sim_result[pe_index[k],]-1.96*ese, sim_result[pe_index[k],]+1.96*ese)
    coverage = apply(ci, 1,
                     function(interval) ifelse(mu.true >= interval[1] & mu.true <= interval[2], 1, 0))
    cp[k] = mean(na.omit(coverage))
    # se.cp[k] = sd(na.omit(coverage))
  }
  ese = round(mean(ese), 3)
  cp = round(cp, 3)

  return(c(bias, esd, ese, mse, cp))
}

summarize_all_model_combinations_and_sample_sizes = function(setting,
                                                             scenario,
                                                             J,
                                                             missing_rate,
                                                             n.vector,
                                                             replicate_num,
                                                             mu.true,
                                                             version,
                                                             is.original){

  result = matrix(NA, 8*length(n.vector),  5)
  j = 1
  for(n in n.vector){
    for(model_num in 1:J){
      model_combinations = combn(J, model_num)
      for(i in 1:ncol(model_combinations)){
        # print(c(model_num, i))
        model_set = model_combinations[, i]
        read_file = paste0(c("Simulation_Results/EBMR_IPW_", setting, "-", missing_rate, "-scenario", scenario, "_", model_set, "_n", n, "_replicate", replicate_num, "_", version, ".RDS"), collapse = "")
        sim_result = readRDS(read_file)

        if(model_num == 1 & i == 1){
          result[j, ] = summarize_results(sim_result, pe_index = 2, ese_index = 4, mu.true, is.original)
          j = j + 1
        }

        result[j, ] = summarize_results(sim_result, pe_index = 1, ese_index = 3, mu.true, is.original)
        j = j + 1
      }
    }
  }
  return(result)
}

summarize_all_settings_with_all_missing_rates = function(settings,
                                                         missing_rates,
                                                         scenario,
                                                         J,
                                                         n.vector,
                                                         all_data_file.list,
                                                         version,
                                                         is.original = TRUE){
  for(setting in settings){
    results_with_all_missing_rates = matrix(NA, 8*length(n.vector), 5*length(missing_rates))
    mu.true = mean(readRDS(all_data_file.list[[setting]][[1]][[1]])$y)
    for(i in 1:length((missing_rates))){
      results_with_all_missing_rates[, ((i-1)*5+1):((i-1)*5+5)] = summarize_all_model_combinations_and_sample_sizes(
        setting = setting,
        scenario = scenario,
        J = J,
        missing_rate = missing_rates[i],
        n.vector = n.vector,
        replicate_num = replicate_num,
        mu.true = mu.true,
        is.original,
        version = version
      )
    }
    estimator_names= rep(c("$\\hat{\\mu}_\\text{IPW}$",
                           "$\\hat{\\mu}_{100}$", "$\\hat{\\mu}_{010}$", "$\\hat{\\mu}_{001}$",
                           "$\\hat{\\mu}_{110}$", "$\\hat{\\mu}_{101}$", "$\\hat{\\mu}_{011}$",
                           "$\\hat{\\mu}_{111}$"), length(n.vector))

    results_with_all_missing_rates = cbind(estimator_names, as.data.frame(results_with_all_missing_rates)) %>%
      as.data.frame
    colnames(results_with_all_missing_rates) = c("", rep(c("Bias", "ESD", "ESE", "MSE", "CP"), length(missing_rates)))

    # print(results_with_all_missing_rates)
    print(
      kable(results_with_all_missing_rates, align = "c", booktabs = TRUE, escape = FALSE, linesep = "") %>%
        kable_styling(full_width = FALSE, latex_options = c("hold_position")) %>%
        add_header_above(c("", "$50\\%$ missing" = 5, "$30\\%$ missing" = 5))
    )
    # return(
    #   kable(results_with_all_missing_rates, align = "c", booktabs = TRUE, escape = FALSE, linesep = "") %>%
    #     kable_styling(full_width = FALSE, latex_options = c("hold_position")) %>%
    #     add_header_above(c("", "$50\\%$ missing" = 5, "$30\\%$ missing" = 5))
    # )
  }
}
