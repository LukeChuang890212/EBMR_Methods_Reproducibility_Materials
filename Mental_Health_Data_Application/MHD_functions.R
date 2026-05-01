#------------------------------------------------------------------------------#
# Mental Health Data Application - Utility Functions
# This file contains all helper functions used in Mental_Health_Data_Analysis.R
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Data Generation
#------------------------------------------------------------------------------#

#' Generate dataset from aggregated data
#'
#' @param original_data Data frame with aggregated data and percentages
#' @param class_n Vector of counts for each class
#' @param n Total sample size
#' @return Data frame with individual-level data
gen_data <- function(original_data, class_n, n = 2486) {
  dat <- matrix(NA, n, 5)
  for (i in 1:nrow(original_data)) {
    dat[((sum(class_n[0:(i-1)])) + 1):sum(class_n[1:i]), ] <-
      matrix(rep(unlist(original_data[i, ]), each = class_n[i]), class_n[i], 5)
  }
  dat <- as.data.frame(dat[, -5])
  colnames(dat) <- c("father", "health", "teacher_report", "parent_report")

  dat$father[dat$father == "Yes"] <- 0
  dat$father[dat$father == "No"] <- 1
  dat$health[dat$health == "No"] <- 0
  dat$health[dat$health == "Yes"] <- 1
  dat[dat == "Abnormal"] <- 1
  dat[dat == "Normal"] <- 0
  dat[dat == "Missing"] <- NA
  dat$r <- ifelse(is.na(dat$teacher_report), 0, 1)
  dat$father <- as.numeric(dat$father)
  dat$health <- as.numeric(dat$health)
  dat$parent_report <- as.numeric(dat$parent_report)
  dat$teacher_report <- as.numeric(dat$teacher_report)
  dat$fp <- dat$father * dat$parent_report
  dat$fh <- dat$father * dat$health
  dat$hp <- dat$health * dat$parent_report
  dat$fhp <- dat$father * dat$health * dat$parent_report

  dat$teacher_report[dat$r == 0] <- -1
  dat$y <- dat$teacher_report
  return(dat)
}

#------------------------------------------------------------------------------#
# Effect Estimation Functions
#------------------------------------------------------------------------------#

#' Log odds ratio effect function
#'
#' @param x Probability for group 0 (unexposed)
#' @param y Probability for group 1 (exposed)
#' @return Log odds ratio
effect.f <- function(x, y) {
  log((y / (1 - y)) / (x / (1 - x)))
}

#' Delta method gradient for log odds ratio
#'
#' @param x Probability for group 0
#' @param y Probability for group 1
#' @return Gradient vector for delta method SE computation
delta_method <- function(x, y) {
  c(-(y / (1 - y)) / ((x / (1 - x))^2) * ((1 / x^2) / ((1 / x - 1)^2)),
    ((1 / y^2) / ((1 / y - 1)^2)) / (x / (1 - x))) / exp(effect.f(x, y))
}

#------------------------------------------------------------------------------#
# Summary and Utility Functions
#------------------------------------------------------------------------------#

#' Summarize GMM propensity score model fits
#'
#' @param ps_fit.list List of propensity score model fits
#' @return List of summary tables with estimates, SE, and p-values
summary.gmm <- function(ps_fit.list) {
  alpha.list <- lapply(ps_fit.list, function(ps_fit) ps_fit$coef)
  alpha_se.list <- lapply(ps_fit.list, function(ps_fit) ps_fit$se)
  for (i in 1:length(alpha.list)) {
    alpha.list[[i]] <- rbind(alpha.list[[i]], alpha_se.list[[i]])
    p.value <- (1 - pnorm(abs(alpha.list[[i]][1, ] / alpha.list[[i]][2, ]))) * 2
    alpha.list[[i]] <- rbind(alpha.list[[i]], p.value)
    rownames(alpha.list[[i]]) <- c("Estimate", "SE", "p-value")
    # Check if model is MNAR (has teacher_report) or MAR
    if (ps_fit.list[[i]]$is.mnar) {
      colnames(alpha.list[[i]]) <- c("Intercept", "teacher_report", ps_fit.list[[i]]$model_x_names)
    } else {
      colnames(alpha.list[[i]]) <- c("Intercept", ps_fit.list[[i]]$model_x_names)
    }
  }
  return(alpha.list)
}

#' Print propensity score model summary in a beautiful table format
#'
#' @param ps_fit.list List of propensity score model fits
#' @param health_label Label for the health group (e.g., "Health = 0")
print_ps_summary <- function(ps_fit.list, health_label = "") {
  alpha.list <- summary.gmm(ps_fit.list)

  # Model descriptions
  model_names <- c(
    "Model 1: MNAR (Y, Father)",
    "Model 2: MNAR (Y, Parent Report)",
    "Model 3: MAR (Father, Parent Report)"
  )

  # Print header
  if (health_label != "") {
    cat("\n")
    cat("  ", paste(rep("=", 76), collapse = ""), "\n")
    cat("  ", health_label, "\n")
    cat("  ", paste(rep("=", 76), collapse = ""), "\n")
  }

  for (i in 1:length(alpha.list)) {
    mat <- alpha.list[[i]]
    col_names <- colnames(mat)
    n_cols <- ncol(mat)

    # Model header
    cat("\n  ", paste(rep("-", 76), collapse = ""), "\n")
    cat("  ", model_names[i], "\n")
    cat("  ", paste(rep("-", 76), collapse = ""), "\n\n")

    # Table header with variable names
    cat("  ", sprintf("%-14s", ""), sep = "")
    for (j in 1:n_cols) {
      cat(sprintf("%14s", col_names[j]))
    }
    cat("\n")
    cat("  ", paste(rep("-", 14 + 14 * n_cols), collapse = ""), "\n")

    # Estimate row
    cat("  ", sprintf("%-14s", "Estimate"), sep = "")
    for (j in 1:n_cols) {
      cat(sprintf("%14.4f", mat["Estimate", j]))
    }
    cat("\n")

    # SE row
    cat("  ", sprintf("%-14s", "Std. Error"), sep = "")
    for (j in 1:n_cols) {
      cat(sprintf("%14.4f", mat["SE", j]))
    }
    cat("\n")

    # z-value row
    cat("  ", sprintf("%-14s", "z value"), sep = "")
    for (j in 1:n_cols) {
      z_val <- mat["Estimate", j] / mat["SE", j]
      cat(sprintf("%14.3f", z_val))
    }
    cat("\n")

    # p-value row with significance stars
    cat("  ", sprintf("%-14s", "Pr(>|z|)"), sep = "")
    for (j in 1:n_cols) {
      p_val <- mat["p-value", j]
      stars <- ""
      if (p_val < 0.001) stars <- "***"
      else if (p_val < 0.01) stars <- "**"
      else if (p_val < 0.05) stars <- "*"
      else if (p_val < 0.1) stars <- "."
      cat(sprintf("%10.4f%-4s", p_val, stars))
    }
    cat("\n")
    cat("  ", paste(rep("-", 14 + 14 * n_cols), collapse = ""), "\n")
  }

  # Significance codes legend
  cat("\n  Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
}

#' Remove extreme values using IQR method
#'
#' @param v Numeric vector
#' @param multiplier IQR multiplier (default 3 for extreme outliers)
#' @return Vector with extreme values removed
rm.extreme <- function(v, multiplier = 3) {
  iqr <- IQR(v)
  lower_bound <- quantile(v, 0.25) - multiplier * iqr
  upper_bound <- quantile(v, 0.75) + multiplier * iqr
  return(v[v >= lower_bound & v <= upper_bound])
}

#------------------------------------------------------------------------------#
# Pretty Printing Functions
#------------------------------------------------------------------------------#

#' Print section header
#'
#' @param title Section title
#' @param width Total width (default 80)
print_header <- function(title, width = 80) {
  cat("\n")
  cat(paste(rep("=", width), collapse = ""), "\n")
  padding <- (width - nchar(title)) %/% 2
  cat(paste(rep(" ", padding), collapse = ""), title, "\n")
  cat(paste(rep("=", width), collapse = ""), "\n\n")
}

#' Print subsection header
#'
#' @param title Subsection title
#' @param width Total width (default 80)
print_section <- function(title, width = 80) {
  cat("\n")
  cat(paste(rep("-", width), collapse = ""), "\n")
  cat(" ", title, "\n")
  cat(paste(rep("-", width), collapse = ""), "\n")
}

#' Format confidence interval for display
#'
#' @param lower Lower bound
#' @param upper Upper bound
#' @param digits Number of decimal places (default 4)
#' @return Formatted string "[lower, upper]"
format_ci <- function(lower, upper, digits = 4) {
  sprintf("[%.*f, %.*f]", digits, lower, digits, upper)
}

#------------------------------------------------------------------------------#
# Bootstrap Functions
#------------------------------------------------------------------------------#

#' Standard Bootstrap for EBMR estimation
#'
#' @param full_ps_specifications List of propensity score model specifications
#' @param W Weighting function for GMM
#' @param data Full dataset
#' @param B Number of bootstrap samples
#' @param init.list Initial values for optimization by health group
#' @param nu_init.list Initial nu values for ensemble step
#' @return Matrix of bootstrap results (rows = estimands, cols = bootstrap samples)
bootstrap <- function(full_ps_specifications, W, data, B, init.list, nu_init.list) {
  n <- nrow(data)

  # Extract nu estimates for each health group (from main analysis)
  nu_hat0 <- nu_init.list[[1]]
  nu_hat1 <- nu_init.list[[2]]

  # Pre-compute model combinations using nu estimates as initial values
  model_configs <- list(
    list(model_set = 1, nu_init0 = 1, nu_init1 = 1),
    list(model_set = 2, nu_init0 = 1, nu_init1 = 1),
    list(model_set = 3, nu_init0 = 1, nu_init1 = 1),
    list(model_set = c(1, 2), nu_init0 = nu_hat0[1:2], nu_init1 = nu_hat1[1:2]),
    list(model_set = c(1, 3), nu_init0 = nu_hat0[c(1,3)], nu_init1 = nu_hat1[c(1,3)]),
    list(model_set = c(2, 3), nu_init0 = nu_hat0[2:3], nu_init1 = nu_hat1[2:3]),
    list(model_set = c(1, 2, 3), nu_init0 = nu_hat0, nu_init1 = nu_hat1)
  )

  # Define h_nu function once
  h_nu <- function(data) cbind(f = data$father, p = data$parent_report, fp = data$fp)

  # Setup parallel cluster
  cores <- parallel::detectCores()
  n_cores <- max(1, cores - 2)
  cl <- parallel::makeCluster(n_cores)
  doSNOW::registerDoSNOW(cl)

  pb <- txtProgressBar(max = B, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  parallel_packages <- c("EBMRalgorithmFast2", "stringr", "Matrix", "numDeriv")

  boot_result <- foreach::foreach(
    i = 1:B,
    .combine = 'cbind',
    .options.snow = opts,
    .packages = parallel_packages,
    .export = c("effect.f", "delta_method", "model_configs", "h_nu")
  ) %dopar% {

    tryCatch({
      # Resample indices
      indices <- sample.int(n, replace = TRUE)
      boot_dat <- data[indices, ]

      boot_health0 <- boot_dat$health == 0
      boot_health1 <- !boot_health0

      theta.hat <- numeric(7)
      se_theta.hat <- numeric(7)
      w.hat.list <- vector("list", 7)
      nu.hat.list <- vector("list", 7)
      alpha_hat <- numeric(16)
      alpha_se_hat <- numeric(16)

      sub_boot_dat0 <- boot_dat[boot_health0, ]
      sub_boot_dat1 <- boot_dat[boot_health1, ]
      sub_boot_list <- list(sub_boot_dat0, sub_boot_dat1)

      r0 <- sub_boot_dat0$r == 1
      r1 <- sub_boot_dat1$r == 1
      p_cc <- c(mean(sub_boot_dat0$teacher_report[r0]),
                mean(sub_boot_dat1$teacher_report[r1]))
      se_p_cc <- c(sd(sub_boot_dat0$teacher_report[r0]) / sqrt(sum(r0)),
                   sd(sub_boot_dat1$teacher_report[r1]) / sqrt(sum(r1)))

      for (k in 1:7) {
        config <- model_configs[[k]]
        model_set <- config$model_set
        nu_init_by_health <- list(config$nu_init0, config$nu_init1)

        p.hat <- numeric(2)
        se_p.hat <- numeric(2)
        w.hat <- c()
        nu.hat <- c()

        for (l in 1:2) {
          sub_boot_dat <- sub_boot_list[[l]]
          nu_init <- nu_init_by_health[[l]]

          ps_specifications <- list(
            formula.list = full_ps_specifications$formula.list[model_set],
            h_x_names.list = full_ps_specifications$h_x_names.list[model_set],
            alpha_init.list = init.list[[l]][model_set],
            outcome = full_ps_specifications$outcome,
            inv_link = full_ps_specifications$inv_link
          )

          ebmr <- EBMRAlgorithmFast2$new("teacher_report", ps_specifications,
                                         sub_boot_dat, W)
          result_ebmr <- ebmr$EBMR_IPW(h_nu = h_nu, nu_init = nu_init, true_ps = NULL, se.fit = FALSE)

          p.hat[l] <- result_ebmr$mu_ipw
          se_p.hat[l] <- result_ebmr$se_ipw
          w.hat <- c(w.hat, result_ebmr$w.hat)
          nu.hat <- c(nu.hat, result_ebmr$nu.hat)

          if (k == 7) {
            idx_start <- (l - 1) * 8 + 1
            idx_end <- l * 8
            alpha_hat[idx_start:idx_end] <- unlist(lapply(ebmr$ps_fit.list,
                                                           function(ps_fit) ps_fit$coefficients))
            alpha_se_hat[idx_start:idx_end] <- unlist(lapply(ebmr$ps_fit.list,
                                                              function(ps_fit) ps_fit$se))
          }
        }

        theta.hat[k] <- effect.f(p.hat[1], p.hat[2])
        se_theta.hat[k] <- sqrt(sum((delta_method(p.hat[1], p.hat[2])^2) * (se_p.hat^2)))
        w.hat.list[[k]] <- w.hat
        nu.hat.list[[k]] <- nu.hat
      }

      theta_cc <- effect.f(p_cc[1], p_cc[2])
      se_theta_cc <- sqrt(sum((delta_method(p_cc[1], p_cc[2])^2) * (se_p_cc^2)))

      c(theta.hat = theta.hat, theta_cc = theta_cc,
        se_theta.hat = se_theta.hat, se_theta_cc = se_theta_cc,
        unlist(w.hat.list), unlist(nu.hat.list),
        alpha_hat0 = alpha_hat[1:8], alpha_hat1 = alpha_hat[9:16],
        alpha_se_hat0 = alpha_se_hat[1:8], alpha_se_hat1 = alpha_se_hat[9:16])

    }, error = function(e) {
      rep(NA, 8 + 8 + sum(sapply(1:7, function(k) 2 * length(model_configs[[k]]$model_set))) * 2 + 32)
    })
  }

  close(pb)
  parallel::stopCluster(cl)
  gc()
  return(boot_result)
}

#' Perturbation Bootstrap for EBMR estimation
#'
#' Uses random weights instead of resampling - more stable for singular cases
#'
#' @param full_ps_specifications List of propensity score model specifications
#' @param W Weighting function for GMM
#' @param data Full dataset
#' @param B Number of bootstrap samples
#' @param init.list Initial values for optimization by health group
#' @param nu_init.list Initial nu values for ensemble step
#' @return Matrix of bootstrap results
perturbation_bootstrap <- function(full_ps_specifications, W, data, B, init.list, nu_init.list) {
  n <- nrow(data)

  nu_hat0 <- nu_init.list[[1]]
  nu_hat1 <- nu_init.list[[2]]

  model_configs <- list(
    list(model_set = 1, nu_init0 = 1, nu_init1 = 1),
    list(model_set = 2, nu_init0 = 1, nu_init1 = 1),
    list(model_set = 3, nu_init0 = 1, nu_init1 = 1),
    list(model_set = c(1, 2), nu_init0 = nu_hat0[1:2], nu_init1 = nu_hat1[1:2]),
    list(model_set = c(1, 3), nu_init0 = nu_hat0[c(1,3)], nu_init1 = nu_hat1[c(1,3)]),
    list(model_set = c(2, 3), nu_init0 = nu_hat0[2:3], nu_init1 = nu_hat1[2:3]),
    list(model_set = c(1, 2, 3), nu_init0 = nu_hat0, nu_init1 = nu_hat1)
  )

  h_nu <- function(data) cbind(f = data$father, p = data$parent_report, fp = data$fp)

  cores <- parallel::detectCores()
  n_cores <- max(1, cores - 2)
  cl <- parallel::makeCluster(n_cores)
  doSNOW::registerDoSNOW(cl)

  pb <- txtProgressBar(max = B, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  parallel_packages <- c("EBMRalgorithmFast2", "stringr", "Matrix", "numDeriv")

  boot_result <- foreach::foreach(
    i = 1:B,
    .combine = 'cbind',
    .options.snow = opts,
    .packages = parallel_packages,
    .export = c("effect.f", "delta_method", "model_configs", "h_nu")
  ) %dopar% {

    tryCatch({
      wt <- rexp(n, rate = 1)

      theta.hat <- numeric(7)
      se_theta.hat <- numeric(7)
      w.hat.list <- vector("list", 7)
      nu.hat.list <- vector("list", 7)
      alpha_hat <- numeric(16)
      alpha_se_hat <- numeric(16)

      health0_idx <- data$health == 0
      health1_idx <- data$health == 1
      sub_dat0 <- data[health0_idx, ]
      sub_dat1 <- data[health1_idx, ]
      sub_dat_list <- list(sub_dat0, sub_dat1)

      wt0 <- wt[health0_idx]
      wt1 <- wt[health1_idx]
      wt_list <- list(wt0, wt1)

      r0 <- sub_dat0$r == 1
      r1 <- sub_dat1$r == 1
      p_cc <- c(
        sum(wt0[r0] * sub_dat0$teacher_report[r0]) / sum(wt0[r0]),
        sum(wt1[r1] * sub_dat1$teacher_report[r1]) / sum(wt1[r1])
      )
      se_p_cc <- c(
        sqrt(sum(wt0[r0] * (sub_dat0$teacher_report[r0] - p_cc[1])^2) / sum(wt0[r0])) / sqrt(sum(r0)),
        sqrt(sum(wt1[r1] * (sub_dat1$teacher_report[r1] - p_cc[2])^2) / sum(wt1[r1])) / sqrt(sum(r1))
      )

      for (k in 1:7) {
        config <- model_configs[[k]]
        model_set <- config$model_set
        nu_init_by_health <- list(config$nu_init0, config$nu_init1)

        p.hat <- numeric(2)
        se_p.hat <- numeric(2)
        w.hat <- c()
        nu.hat <- c()

        for (l in 1:2) {
          sub_dat <- sub_dat_list[[l]]
          sub_wt <- wt_list[[l]]
          nu_init <- nu_init_by_health[[l]]

          ps_specifications <- list(
            formula.list = full_ps_specifications$formula.list[model_set],
            h_x_names.list = full_ps_specifications$h_x_names.list[model_set],
            alpha_init.list = init.list[[l]][model_set],
            outcome = full_ps_specifications$outcome,
            inv_link = full_ps_specifications$inv_link
          )

          ebmr <- EBMRAlgorithmFast2$new("teacher_report", ps_specifications,
                                         sub_dat, W, wt = sub_wt)
          result_ebmr <- ebmr$EBMR_IPW(h_nu = h_nu, nu_init = nu_init, true_ps = NULL,
                                        se.fit = FALSE, wt = sub_wt)

          p.hat[l] <- result_ebmr$mu_ipw
          se_p.hat[l] <- result_ebmr$se_ipw
          w.hat <- c(w.hat, result_ebmr$w.hat)
          nu.hat <- c(nu.hat, result_ebmr$nu.hat)

          if (k == 7) {
            idx_start <- (l - 1) * 8 + 1
            idx_end <- l * 8
            alpha_hat[idx_start:idx_end] <- unlist(lapply(ebmr$ps_fit.list,
                                                           function(ps_fit) ps_fit$coefficients))
            alpha_se_hat[idx_start:idx_end] <- unlist(lapply(ebmr$ps_fit.list,
                                                              function(ps_fit) ps_fit$se))
          }
        }

        theta.hat[k] <- effect.f(p.hat[1], p.hat[2])
        se_theta.hat[k] <- sqrt(sum((delta_method(p.hat[1], p.hat[2])^2) * (se_p.hat^2)))
        w.hat.list[[k]] <- w.hat
        nu.hat.list[[k]] <- nu.hat
      }

      theta_cc <- effect.f(p_cc[1], p_cc[2])
      se_theta_cc <- sqrt(sum((delta_method(p_cc[1], p_cc[2])^2) * (se_p_cc^2)))

      c(theta.hat = theta.hat, theta_cc = theta_cc,
        se_theta.hat = se_theta.hat, se_theta_cc = se_theta_cc,
        unlist(w.hat.list), unlist(nu.hat.list),
        alpha_hat0 = alpha_hat[1:8], alpha_hat1 = alpha_hat[9:16],
        alpha_se_hat0 = alpha_se_hat[1:8], alpha_se_hat1 = alpha_se_hat[9:16])

    }, error = function(e) {
      rep(NA, 8 + 8 + sum(sapply(1:7, function(k) 2 * length(model_configs[[k]]$model_set))) * 2 + 32)
    })
  }

  close(pb)
  parallel::stopCluster(cl)
  gc()
  return(boot_result)
}

#' Standard Bootstrap for separate health groups and model sets
#'
#' This allows computing SE for all 49 theta combinations by running
#' separate bootstraps for each health group and model set.
#' Uses resampling (like standard bootstrap) instead of weights.
#'
#' @param full_ps_specifications List of propensity score model specifications
#' @param W Weighting function for GMM
#' @param data Full dataset
#' @param B Number of bootstrap samples
#' @param init.list Initial values for optimization by health group
#' @param nu_init.list Initial nu values for ensemble step
#' @param health_value Health group value (0 or 1)
#' @param model_set_idx Index of model set in all_model_sets
#' @param all_model_sets List of all model sets
#' @return Vector of B mu_ipw replicates
bootstrap_separate <- function(full_ps_specifications, W, data, B, init.list,
                               nu_init.list, health_value, model_set_idx,
                               all_model_sets) {
  n <- nrow(data)
  health_idx <- data$health == health_value
  sub_dat <- data[health_idx, ]
  n_sub <- nrow(sub_dat)

  model_set <- all_model_sets[[model_set_idx]]

  nu_hat_full <- nu_init.list[[health_value + 1]]
  if (length(model_set) == 1) {
    nu_init <- 1
  } else if (length(model_set) == 2) {
    nu_init <- nu_hat_full[model_set]
  } else {
    nu_init <- nu_hat_full
  }

  h_nu <- function(data) cbind(f = data$father, p = data$parent_report, fp = data$fp)

  cores <- parallel::detectCores()
  n_cores <- max(1, cores - 2)
  cl <- parallel::makeCluster(n_cores)
  doSNOW::registerDoSNOW(cl)

  pb <- txtProgressBar(max = B, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  parallel_packages <- c("EBMRalgorithmFast2", "stringr", "Matrix", "numDeriv")

  mu_replicates <- foreach::foreach(
    i = 1:B,
    .combine = 'c',
    .options.snow = opts,
    .packages = parallel_packages
  ) %dopar% {

    tryCatch({
      # Resample with replacement (standard bootstrap)
      boot_indices <- sample.int(n_sub, replace = TRUE)
      boot_dat <- sub_dat[boot_indices, ]

      ps_specifications <- list(
        formula.list = full_ps_specifications$formula.list[model_set],
        h_x_names.list = full_ps_specifications$h_x_names.list[model_set],
        alpha_init.list = init.list[[health_value + 1]][model_set],
        outcome = full_ps_specifications$outcome,
        inv_link = full_ps_specifications$inv_link
      )

      ebmr <- EBMRAlgorithmFast2$new("teacher_report", ps_specifications,
                                     boot_dat, W)
      result_ebmr <- ebmr$EBMR_IPW(h_nu = h_nu, nu_init = nu_init, true_ps = NULL,
                                    se.fit = FALSE)

      result_ebmr$mu_ipw

    }, error = function(e) {
      NA
    })
  }

  close(pb)
  parallel::stopCluster(cl)
  gc()

  return(mu_replicates)
}

#' Perturbation Bootstrap for separate health groups and model sets
#'
#' This allows computing SE for all 49 theta combinations by running
#' separate bootstraps for each health group and model set.
#'
#' @param full_ps_specifications List of propensity score model specifications
#' @param W Weighting function for GMM
#' @param data Full dataset
#' @param B Number of bootstrap samples
#' @param init.list Initial values for optimization by health group
#' @param nu_init.list Initial nu values for ensemble step
#' @param health_value Health group value (0 or 1)
#' @param model_set_idx Index of model set in all_model_sets
#' @param all_model_sets List of all model sets
#' @return Vector of B mu_ipw replicates
perturbation_bootstrap_separate <- function(full_ps_specifications, W, data, B, init.list,
                                            nu_init.list, health_value, model_set_idx,
                                            all_model_sets) {
  n <- nrow(data)
  health_idx <- data$health == health_value
  sub_dat <- data[health_idx, ]
  n_sub <- nrow(sub_dat)

  model_set <- all_model_sets[[model_set_idx]]

  nu_hat_full <- nu_init.list[[health_value + 1]]
  if (length(model_set) == 1) {
    nu_init <- 1
  } else if (length(model_set) == 2) {
    nu_init <- nu_hat_full[model_set]
  } else {
    nu_init <- nu_hat_full
  }

  h_nu <- function(data) cbind(f = data$father, p = data$parent_report, fp = data$fp)

  cores <- parallel::detectCores()
  n_cores <- max(1, cores - 2)
  cl <- parallel::makeCluster(n_cores)
  doSNOW::registerDoSNOW(cl)

  pb <- txtProgressBar(max = B, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  parallel_packages <- c("EBMRalgorithmFast2", "stringr", "Matrix", "numDeriv")

  mu_replicates <- foreach::foreach(
    i = 1:B,
    .combine = 'c',
    .options.snow = opts,
    .packages = parallel_packages
  ) %dopar% {

    tryCatch({
      wt <- rexp(n_sub, rate = 1)

      ps_specifications <- list(
        formula.list = full_ps_specifications$formula.list[model_set],
        h_x_names.list = full_ps_specifications$h_x_names.list[model_set],
        alpha_init.list = init.list[[health_value + 1]][model_set],
        outcome = full_ps_specifications$outcome,
        inv_link = full_ps_specifications$inv_link
      )

      ebmr <- EBMRAlgorithmFast2$new("teacher_report", ps_specifications,
                                     sub_dat, W, wt = wt)
      result_ebmr <- ebmr$EBMR_IPW(h_nu = h_nu, nu_init = nu_init, true_ps = NULL,
                                    se.fit = FALSE, wt = wt)

      result_ebmr$mu_ipw

    }, error = function(e) {
      NA
    })
  }

  close(pb)
  parallel::stopCluster(cl)
  gc()

  return(mu_replicates)
}

#------------------------------------------------------------------------------#
# Visualization Functions
#------------------------------------------------------------------------------#

#' Generate heatmap of theta estimates
#'
#' Creates a heatmap visualization highlighting specified model combinations
#'
#' @param result_matrix 7x7 matrix of theta estimates
#' @param model_set_labels Labels for model sets
#' @param model3_included Indices of rows with Model 3 included
#' @param model1_included Indices of columns with Model 1 included
#' @param output_file Path to save PNG file
generate_theta_heatmap <- function(result_matrix, model_set_labels,
                                    model3_included = c(3, 5, 6, 7),
                                    model1_included = c(1, 4, 5, 7),
                                    output_file = "MHD_results/theta_heatmap.png") {

  # Create highlight matrix
  highlight_matrix <- matrix(0, 7, 7)
  for (i in 1:7) {
    for (j in 1:7) {
      has_model3 <- i %in% model3_included
      has_model1 <- j %in% model1_included
      if (has_model3 && has_model1) {
        highlight_matrix[i, j] <- 2
      } else if (has_model3 || has_model1) {
        highlight_matrix[i, j] <- 1
      }
    }
  }

  # Save heatmap as PNG
  png(output_file, width = 10, height = 8, units = "in", res = 300)

  layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1))

  n_colors <- 100
  pe_range <- range(result_matrix)
  max_abs <- max(abs(pe_range))
  breaks <- seq(-max_abs, max_abs, length.out = n_colors + 1)

  # Muted color palette for non-highlighted cells
  gray_blue <- colorRampPalette(c("#8C96A0", "#A8B0B8", "#C4CCD4", "#E0E4E8", "#F0F0F0"))
  gray_red <- colorRampPalette(c("#F0F0F0", "#E8E0DC", "#D4C4BC", "#C0A8A0", "#A08C88"))
  muted_colors <- c(gray_blue(n_colors/2), gray_red(n_colors/2))

  # Vibrant color palette for highlighted cells
  highlight_blue <- colorRampPalette(c("#1A5276", "#2874A6", "#3498DB", "#85C1E9", "#FFFFFF"))
  highlight_orange <- colorRampPalette(c("#FFFFFF", "#FAD7A0", "#F5B041", "#E67E22", "#D35400"))
  highlight_colors <- c(highlight_blue(n_colors/2), highlight_orange(n_colors/2))

  get_color <- function(val, highlighted) {
    idx <- findInterval(val, breaks, all.inside = TRUE)
    if (highlighted) {
      return(highlight_colors[idx])
    } else {
      return(muted_colors[idx])
    }
  }

  row_labels <- c("{1}", "{2}", "{3}", "{1,2}", "{1,3}", "{2,3}", "{1,2,3}")
  col_labels <- row_labels

  par(mar = c(6, 7, 5, 1))
  plot(NULL, xlim = c(0.5, 7.5), ylim = c(0.5, 7.5), xlab = "", ylab = "",
       xaxt = "n", yaxt = "n", bty = "n", asp = 1)

  # Draw non-highlighted cells
  for (i in 1:7) {
    for (j in 1:7) {
      y_pos <- 8 - i
      is_highlighted <- highlight_matrix[i, j] == 2

      if (!is_highlighted) {
        cell_col <- get_color(result_matrix[i, j], FALSE)
        rect(j - 0.5, y_pos - 0.5, j + 0.5, y_pos + 0.5,
             col = cell_col, border = "#CCCCCC", lwd = 0.5)
        text(j, y_pos, sprintf("%.3f", result_matrix[i, j]), cex = 0.75, col = "#666666", font = 1)
      }
    }
  }

  # Draw highlighted cells
  for (i in 1:7) {
    for (j in 1:7) {
      y_pos <- 8 - i
      is_highlighted <- highlight_matrix[i, j] == 2

      if (is_highlighted) {
        cell_col <- get_color(result_matrix[i, j], TRUE)
        rect(j - 0.52, y_pos - 0.52, j + 0.52, y_pos + 0.52,
             col = NA, border = "#FFD700", lwd = 6)
        rect(j - 0.5, y_pos - 0.5, j + 0.5, y_pos + 0.5,
             col = cell_col, border = "#000000", lwd = 3)
        text_col <- ifelse(abs(result_matrix[i, j]) > max_abs * 0.5, "white", "black")
        text(j, y_pos, sprintf("%.3f", result_matrix[i, j]), cex = 1.0, col = text_col, font = 2)
      }
    }
  }

  # Add axes labels
  axis_col_default <- "gray40"
  axis_col_highlight <- "#D35400"

  for (j in 1:7) {
    col <- ifelse(j %in% model1_included, axis_col_highlight, axis_col_default)
    font_style <- ifelse(j %in% model1_included, 2, 1)
    mtext(col_labels[j], side = 1, at = j, line = 0.5, cex = 1.0, col = col, font = font_style)
  }

  for (i in 1:7) {
    y_pos <- 8 - i
    col <- ifelse(i %in% model3_included, axis_col_highlight, axis_col_default)
    font_style <- ifelse(i %in% model3_included, 2, 1)
    mtext(row_labels[i], side = 2, at = y_pos, line = 0.5, cex = 1.0, col = col, font = font_style, las = 1)
  }

  mtext("Model Set for Health = 1", side = 1, line = 4.5, cex = 1.2, font = 2)
  mtext("Model Set for Health = 0", side = 2, line = 5, cex = 1.2, font = 2)

  mtext(expression(paste("Log Odds Ratio Estimates (", hat(theta), ") by Model Combination")),
        side = 3, line = 2.5, cex = 1.4, font = 2)
  mtext("Highlighted: Model 3 (MAR) for Health=0, Model 1 (MNAR:father) for Health=1",
        side = 3, line = 0.8, cex = 0.95, font = 3)

  # Color legend
  par(mar = c(6, 1, 5, 3))
  legend_y <- seq(0, 1, length.out = n_colors)
  image(1, legend_y, matrix(legend_y, nrow = 1), col = highlight_colors,
        xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  box(lwd = 2)

  legend_ticks <- pretty(c(-max_abs, max_abs), n = 5)
  legend_ticks <- legend_ticks[legend_ticks >= -max_abs & legend_ticks <= max_abs]
  legend_pos <- (legend_ticks + max_abs) / (2 * max_abs)
  axis(4, at = legend_pos, labels = sprintf("%.2f", legend_ticks), las = 1, cex.axis = 0.9)
  mtext(expression(hat(theta)), side = 4, line = 2, cex = 1.1)
  mtext("(Highlighted cells)", side = 1, line = 1, cex = 0.8, font = 3)

  dev.off()

  cat("  Heatmap saved to:", output_file, "\n")
}

#------------------------------------------------------------------------------#
# Model Combination Utilities
#------------------------------------------------------------------------------#

#' Generate all model set combinations
#'
#' @return List of 7 model sets: {1}, {2}, {3}, {1,2}, {1,3}, {2,3}, {1,2,3}
generate_all_model_sets <- function() {
  all_model_sets <- list()
  for (model_num in 1:3) {
    model_combinations <- combn(3, model_num)
    for (i in 1:ncol(model_combinations)) {
      all_model_sets <- c(all_model_sets, list(model_combinations[, i]))
    }
  }
  return(all_model_sets)
}

#' Create binary label for model combination
#'
#' @param model_set Vector of model indices (e.g., c(1,3))
#' @return Binary string representation (e.g., "101")
model_set_to_binary <- function(model_set) {
  paste0(as.integer(1 %in% model_set),
         as.integer(2 %in% model_set),
         as.integer(3 %in% model_set))
}

#' Create theta label for combination
#'
#' @param model_set0 Model set for health=0
#' @param model_set1 Model set for health=1
#' @return Formatted label like "theta_{101}^{110}"
create_theta_label <- function(model_set0, model_set1) {
  binary0 <- model_set_to_binary(model_set0)
  binary1 <- model_set_to_binary(model_set1)
  paste0("theta_{", binary0, "}^{", binary1, "}")
}
