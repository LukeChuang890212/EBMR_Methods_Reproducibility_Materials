## Test multiple model/h_alpha combinations to find one with good ESE2/ESD
## For each config: point estimate + 500 bootstrap reps to check ESE vs Boot SE
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Mental_Health_Data_Application")
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)
library(Matrix); library(numDeriv)
library(parallel); library(foreach); library(doSNOW)
source("MHD_functions.R")

original_data <- read.csv("data_application.csv")
percent <- original_data$Percentage
n <- 2486
class_n <- round(n * percent / 100)
dat <- gen_data(original_data, class_n, n)
dat$y <- dat$teacher_report

W <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
mu_cc <- mean(dat$teacher_report[dat$r == 1])

## Define configs to test
configs <- list(
  ## Config A: Original models (with teacher_report interactions), standard h_alpha
  A = list(
    name = "Original (y-interactions, h_alpha=6)",
    formulas = list(
      r ~ teacher_report + health + father + teacher_report:health + teacher_report:father,
      r ~ teacher_report + health + parent_report + teacher_report:health + teacher_report:parent_report,
      r ~ teacher_report + parent_report + father + teacher_report:parent_report + teacher_report:father
    ),
    h_alpha = list(
      c("health", "father", "parent_report", "fp", "fh", "hp"),
      c("health", "father", "parent_report", "fp", "fh", "hp"),
      c("health", "father", "parent_report", "fp", "fh", "hp")
    ),
    h_nu = function(data) cbind(health=data$health, father=data$father, parent_report=data$parent_report,
                                fp=data$fp, fh=data$fh, hp=data$hp, fhp=data$fhp)
  ),

  ## Config B: Original models, saturated h_alpha (add fhp)
  B = list(
    name = "Original (y-interactions, h_alpha=7 saturated)",
    formulas = list(
      r ~ teacher_report + health + father + teacher_report:health + teacher_report:father,
      r ~ teacher_report + health + parent_report + teacher_report:health + teacher_report:parent_report,
      r ~ teacher_report + parent_report + father + teacher_report:parent_report + teacher_report:father
    ),
    h_alpha = list(
      c("health", "father", "parent_report", "fp", "fh", "hp", "fhp"),
      c("health", "father", "parent_report", "fp", "fh", "hp", "fhp"),
      c("health", "father", "parent_report", "fp", "fh", "hp", "fhp")
    ),
    h_nu = function(data) cbind(health=data$health, father=data$father, parent_report=data$parent_report,
                                fp=data$fp, fh=data$fh, hp=data$hp, fhp=data$fhp)
  ),

  ## Config C: Simple models (no interactions), k=3, overid=4
  C = list(
    name = "Simple (no interactions, k=3, overid=4)",
    formulas = list(
      r ~ teacher_report + health,
      r ~ teacher_report + father,
      r ~ teacher_report + parent_report
    ),
    h_alpha = list(
      c("health", "father", "parent_report", "fp", "fh", "hp"),
      c("health", "father", "parent_report", "fp", "fh", "hp"),
      c("health", "father", "parent_report", "fp", "fh", "hp")
    ),
    h_nu = function(data) cbind(health=data$health, father=data$father, parent_report=data$parent_report,
                                fp=data$fp, fh=data$fh, hp=data$hp, fhp=data$fhp)
  ),

  ## Config D: Two covariates, no interactions, k=4, overid=3
  D = list(
    name = "Two covariates (no interactions, k=4, overid=3)",
    formulas = list(
      r ~ teacher_report + health + father,
      r ~ teacher_report + health + parent_report,
      r ~ teacher_report + father + parent_report
    ),
    h_alpha = list(
      c("health", "father", "parent_report", "fp", "fh", "hp"),
      c("health", "father", "parent_report", "fp", "fh", "hp"),
      c("health", "father", "parent_report", "fp", "fh", "hp")
    ),
    h_nu = function(data) cbind(health=data$health, father=data$father, parent_report=data$parent_report,
                                fp=data$fp, fh=data$fh, hp=data$hp, fhp=data$fhp)
  ),

  ## Config E: Covariate interactions (no y-interactions), k=5, overid=2
  E = list(
    name = "Covariate interactions (k=5, overid=2)",
    formulas = list(
      r ~ teacher_report + health + father + health:father,
      r ~ teacher_report + health + parent_report + health:parent_report,
      r ~ teacher_report + father + parent_report + father:parent_report
    ),
    h_alpha = list(
      c("health", "father", "parent_report", "fp", "fh", "hp"),
      c("health", "father", "parent_report", "fp", "fh", "hp"),
      c("health", "father", "parent_report", "fp", "fh", "hp")
    ),
    h_nu = function(data) cbind(health=data$health, father=data$father, parent_report=data$parent_report,
                                fp=data$fp, fh=data$fh, hp=data$hp, fhp=data$fhp)
  ),

  ## Config F: One y-interaction each, k=4, overid=3
  F_cfg = list(
    name = "One y-interaction each (k=4, overid=3)",
    formulas = list(
      r ~ teacher_report + health + teacher_report:health,
      r ~ teacher_report + father + teacher_report:father,
      r ~ teacher_report + parent_report + teacher_report:parent_report
    ),
    h_alpha = list(
      c("health", "father", "parent_report", "fp", "fh", "hp"),
      c("health", "father", "parent_report", "fp", "fh", "hp"),
      c("health", "father", "parent_report", "fp", "fh", "hp")
    ),
    h_nu = function(data) cbind(health=data$health, father=data$father, parent_report=data$parent_report,
                                fp=data$fp, fh=data$fh, hp=data$hp, fhp=data$fhp)
  )
)

## Run each config: point estimates + 500 bootstrap
B <- 500
n_cores <- min(detectCores() - 1, 10)

for (cfg_name in names(configs)) {
  cfg <- configs[[cfg_name]]
  cat(sprintf("\n========== Config %s: %s ==========\n\n", cfg_name, cfg$name))

  ps_spec <- list(
    formula.list = cfg$formulas,
    h_alpha.list = cfg$h_alpha,
    outcome = "teacher_report",
    inv_link = function(eta) 1 / (1 + exp(eta))
  )

  # Point estimates
  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat, W)

    # Model diagnostics
    for (i in 1:3) {
      pf <- ebmr$ps_fit.list[[i]]
      ps <- pf$fitted.values
      cat(sprintf("  M%d: k=%d, h_dim=%d, PS=[%.4f,%.4f], >0.99:%d, <0.01:%d, mu=%.4f\n",
          i, ncol(pf$design_matrix), ncol(pf$h_x),
          min(ps), max(ps), sum(ps>0.99), sum(ps<0.01),
          mean(dat$r * dat$teacher_report / ps)))
    }

    # All 7 model combinations
    all_sets <- list(c(1), c(2), c(3), c(1,2), c(1,3), c(2,3), c(1,2,3))
    set_names <- c("1", "2", "3", "12", "13", "23", "123")
    point_mu <- point_se <- numeric(7)

    cat(sprintf("\n  %-8s %10s %10s\n", "Models", "Estimate", "SE"))
    cat("  ", paste(rep("-", 30), collapse=""), "\n")
    for (j in 1:7) {
      res <- ebmr$EBMR_IPW(h_nu = cfg$h_nu, model_indices = all_sets[[j]], true_ps = NULL)
      point_mu[j] <- res$mu_ipw
      point_se[j] <- res$se_ipw
      cat(sprintf("  %-8s %10.4f %10.4f\n", set_names[j], res$mu_ipw, res$se_ipw))
    }

    # Bootstrap
    cat(sprintf("\n  Running %d bootstrap reps...\n", B))
    cl <- makeCluster(n_cores)
    registerDoSNOW(cl)
    clusterExport(cl, c("dat", "n", "ps_spec", "W", "cfg"), envir = environment())
    clusterEvalQ(cl, devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE))

    pb <- txtProgressBar(max = B, style = 3)
    progress <- function(nn) setTxtProgressBar(pb, nn)
    opts <- list(progress = progress)

    boot_mat <- foreach(b = 1:B, .combine = 'rbind', .options.snow = opts,
                        .packages = c("stringr", "Matrix", "numDeriv")) %dopar% {
      set.seed(12345 + b)
      idx <- sample(1:n, n, replace = TRUE)
      dat_b <- dat[idx, ]
      all_sets <- list(c(1), c(2), c(3), c(1,2), c(1,3), c(2,3), c(1,2,3))
      mu_vec <- rep(NA, 7)
      tryCatch({
        ebmr_b <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat_b, W)
        for (j in 1:7) {
          tryCatch({
            res_b <- ebmr_b$EBMR_IPW(h_nu = cfg$h_nu, model_indices = all_sets[[j]],
                                       true_ps = NULL, se.fit = FALSE)
            mu_vec[j] <- res_b$mu_ipw
          }, error = function(e) {})
        }
      }, error = function(e) {})
      mu_vec
    }
    close(pb)
    stopCluster(cl)

    # Compare ESE vs Boot SE
    cat(sprintf("\n  %-8s %10s %10s %10s %10s\n", "Models", "Estimate", "ESE", "Boot SE", "ESE/Boot"))
    cat("  ", paste(rep("-", 52), collapse=""), "\n")
    for (j in 1:7) {
      boot_v <- boot_mat[, j]
      boot_se <- sd(boot_v, na.rm = TRUE)
      n_na <- sum(is.na(boot_v))
      ratio <- if (boot_se > 0) point_se[j] / boot_se else NA
      cat(sprintf("  %-8s %10.4f %10.4f %10.4f %10.3f  (NA:%d)\n",
          set_names[j], point_mu[j], point_se[j], boot_se, ratio, n_na))
    }

  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
  })
}
