#------------------------------------------------------------------------------#
# AIDS Clinical Trial Data Analysis (ACTG175)
# Estimate population mean of cd496 (CD4 count at 96 weeks)
# Using EPS
# Analysis by treatment arm
#------------------------------------------------------------------------------#

setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/AIDS")
devtools::load_all("../EPS", quiet = TRUE)
library(speff2trial)
library(Matrix)
library(numDeriv)
library(parallel)
library(foreach)
library(doSNOW)

# Load data
data(ACTG175)
dat <- ACTG175

cat("================================================================================\n")
cat("  ACTG175 AIDS Clinical Trial Data Analysis (EPS)\n")
cat("  Goal: Estimate population mean of CD4 count at 96 weeks (cd496)\n")
cat("================================================================================\n\n")

# Data summary
cat("=== Data Summary ===\n")
cat("Total observations:", nrow(dat), "\n")
cat("Respondents (r=1):", sum(dat$r), "\n")
cat("Non-respondents (r=0):", sum(dat$r == 0), "\n")
cat("Response rate:", round(mean(dat$r) * 100, 2), "%\n\n")

# Treatment arms
cat("=== Treatment Arms ===\n")
cat("arms = 0: zidovudine (ZDV) only\n")
cat("arms = 1: ZDV + didanosine (ddI)\n")
cat("arms = 2: ZDV + zalcitabine (ddC)\n")
cat("arms = 3: didanosine (ddI) only\n\n")

# Summary by arm
cat("=== Summary by Treatment Arm ===\n")
for (arm in 0:3) {
  arm_dat <- dat[dat$arms == arm, ]
  cat(sprintf("Arm %d: n=%d, respondents=%d (%.1f%%), mean cd496 (observed)=%.1f\n",
              arm, nrow(arm_dat), sum(arm_dat$r),
              mean(arm_dat$r) * 100,
              mean(arm_dat$cd496[arm_dat$r == 1], na.rm = TRUE)))
}
cat("\n")

# Prepare data: y = cd496, set to -1 for non-respondents
dat$y <- dat$cd496
dat$y[dat$r == 0] <- -1

# Weight matrix function
W <- function(g.matrix) {
  m <- t(g.matrix) %*% g.matrix / nrow(g.matrix)
  tryCatch(solve(m), error = function(e) MASS::ginv(m))
}

# h_alpha: main effects + all pairwise interactions + squared terms
h_alpha <- function(data) {
  d <- data
  # Main effects + pairwise interactions (21 terms)
  mm <- model.matrix(~ (cd40 + cd80 + cd420 + cd820 + age + wtkg)^2, data = d)
  mm <- mm[, -1, drop = FALSE]
  # Add squared terms (6 terms)
  sq <- cbind(cd40_sq = d$cd40^2, cd80_sq = d$cd80^2, cd420_sq = d$cd420^2,
              cd820_sq = d$cd820^2, age_sq = d$age^2, wtkg_sq = d$wtkg^2)
  cbind(mm, sq)  # total: 27 terms
}

# h_nu: h_alpha terms + squared*main (x_i^2 * x_j) + cubic terms (x_i^3)
h_nu <- function(data) {
  d <- data
  # Start with h_alpha terms (27)
  mm <- model.matrix(~ (cd40 + cd80 + cd420 + cd820 + age + wtkg)^2, data = d)
  mm <- mm[, -1, drop = FALSE]
  sq <- cbind(cd40_sq = d$cd40^2, cd80_sq = d$cd80^2, cd420_sq = d$cd420^2,
              cd820_sq = d$cd820^2, age_sq = d$age^2, wtkg_sq = d$wtkg^2)
  # 3-way interactions (20 terms)
  mm3 <- model.matrix(~ (cd40 + cd80 + cd420 + cd820 + age + wtkg)^3, data = d)
  mm3 <- mm3[, -1, drop = FALSE]
  # Keep only 3-way terms (columns not in mm)
  three_way <- mm3[, !(colnames(mm3) %in% colnames(mm)), drop = FALSE]
  # Squared * main effect: x_i^2 * x_j for all i != j (6*5 = 30 terms)
  v <- cbind(d$cd40, d$cd80, d$cd420, d$cd820, d$age, d$wtkg)
  vnames <- c("cd40", "cd80", "cd420", "cd820", "age", "wtkg")
  sq_main_list <- list()
  for (i in 1:6) {
    for (j in 1:6) {
      if (i != j) {
        sq_main_list[[paste0(vnames[i], "sq_", vnames[j])]] <- v[, i]^2 * v[, j]
      }
    }
  }
  sq_main <- do.call(cbind, sq_main_list)
  colnames(sq_main) <- names(sq_main_list)
  # Cubic terms: x_i^3 (6 terms)
  cu <- cbind(cd40_cu = d$cd40^3, cd80_cu = d$cd80^3, cd420_cu = d$cd420^3,
              cd820_cu = d$cd820^3, age_cu = d$age^3, wtkg_cu = d$wtkg^3)
  cbind(mm, sq, three_way, sq_main, cu)
}

# PS model specifications
# M1: r ~ y + cd40 + cd420 + cd80 + age + wtkg
# M2: r ~ y + cd40 + cd420 + cd820 + age + wtkg
# M3: r ~ y + cd40 + cd80 + cd820 + age + wtkg
# Note: y = cd496 with -1 for non-respondents, outcome = "y"
ps_specifications <- list(
  formula.list = list(
    r ~ y + cd40 + cd420 + cd80 + age + wtkg,
    r ~ y + cd40 + cd420 + cd820 + age + wtkg,
    r ~ y + cd40 + cd80 + cd820 + age + wtkg
  ),
  h_alpha.list = list(h_alpha, h_alpha, h_alpha),
  outcome = "y",
  inv_link = function(eta) 1 / (1 + exp(eta))
)

all_model_sets <- list(c(1), c(2), c(3), c(1,2), c(1,3), c(2,3), c(1,2,3))
model_set_labels <- c("100", "010", "001", "110", "101", "011", "111")
model_set_names <- c("M1 only", "M2 only", "M3 only",
                     "M1+M2", "M1+M3", "M2+M3", "M1+M2+M3")

#------------------------------------------------------------------------------#
# Outlier-trimmed Bootstrap SE (IQR x 3, 1% cap)
#------------------------------------------------------------------------------#
compute_trimmed_boot_se <- function(x) {
  x <- x[!is.na(x)]
  n_total <- length(x)
  if (n_total < 10) return(list(se_raw = NA, se_trim = NA, n_na = 0, n_outlier = NA))
  Q1 <- quantile(x, 0.25); Q3 <- quantile(x, 0.75)
  IQR_val <- Q3 - Q1
  is_outlier <- x < (Q1 - 3 * IQR_val) | x > (Q3 + 3 * IQR_val)
  max_remove <- floor(0.01 * n_total)
  if (sum(is_outlier) > max_remove && max_remove > 0) {
    dist_med <- abs(x - median(x))
    oi <- which(is_outlier)
    keep <- oi[order(dist_med[oi], decreasing = TRUE)[(max_remove + 1):length(oi)]]
    is_outlier[keep] <- FALSE
  }
  list(se_raw = sd(x), se_trim = sd(x[!is_outlier]), n_outlier = sum(is_outlier))
}

#------------------------------------------------------------------------------#
# Analysis function for a single arm
#------------------------------------------------------------------------------#
run_arm_analysis <- function(arm_data, arm_name) {
  cat(strrep("=", 80), "\n")
  cat("Analysis for", arm_name, "\n")
  cat(strrep("=", 80), "\n\n")

  n <- nrow(arm_data)
  cat("Sample size:", n, "\n")
  cat("Respondents:", sum(arm_data$r), "(", round(mean(arm_data$r)*100, 1), "%)\n")
  naive_est <- mean(arm_data$y[arm_data$r == 1])
  cat("Naive estimate (observed mean):", round(naive_est, 2), "\n\n")

  # Fit EBMR
  ebmr <- EPS$new("y", ps_specifications, arm_data, W)

  # Model summaries
  cat("=== Model Summaries ===\n\n")
  for (i in 1:3) {
    pf <- ebmr$ps_fit.list[[i]]
    ps <- pf$fitted.values
    formula_str <- deparse(ps_specifications$formula.list[[i]])
    cat(sprintf("Model %d: %s\n", i, formula_str))

    coef <- pf$coefficients
    se <- pf$se
    z_val <- coef / se
    p_val <- 2 * pnorm(-abs(z_val))
    var_names <- colnames(pf$design_matrix)
    max_w <- max(15, max(nchar(var_names)) + 2)

    cat(sprintf("  %*s %10s %10s %10s %10s\n", max_w, "Variable", "Estimate", "SE", "z", "p"))
    for (k in seq_along(coef)) {
      stars <- ""
      if (p_val[k] < 0.001) stars <- "***"
      else if (p_val[k] < 0.01) stars <- "**"
      else if (p_val[k] < 0.05) stars <- "*"
      else if (p_val[k] < 0.1) stars <- "."
      cat(sprintf("  %*s %10.4f %10.4f %10.3f %10.4f%s\n",
          max_w, var_names[k], coef[k], se[k], z_val[k], p_val[k], stars))
    }
    cat(sprintf("  k=%d, h_dim=%d, overid=%d\n", ncol(pf$design_matrix), ncol(pf$h_x),
        ncol(pf$h_x) - ncol(pf$design_matrix)))
    cat(sprintf("  PS range: [%.4f, %.4f], >0.99:%d, <0.01:%d\n",
        min(ps), max(ps), sum(ps > 0.99), sum(ps < 0.01)))
    cat(sprintf("  GMM obj=%.2e, grad=%.2e\n\n",
        pf$gmm_fit$opt$objective, pf$gmm_fit$opt$final_grad_norm))
  }

  # All model combinations
  results_list <- list()
  cat("=== Point Estimates ===\n\n")
  cat(sprintf("  %-15s %12s %12s %20s %30s\n",
              "Models", "Estimate", "SE", "95% CI", "w.hat"))
  cat("  ", paste(rep("-", 95), collapse = ""), "\n")

  cat(sprintf("  %-15s %12.2f %12s %20s\n", "Naive (CC)", naive_est, "-", "-"))

  for (j in seq_along(all_model_sets)) {
    res <- ebmr$EBMR_IPW(h_nu = h_nu, model_indices = all_model_sets[[j]], true_ps = NULL)
    results_list[[j]] <- list(mu = res$mu_ipw, se = res$se_ipw,
                               nu = res$nu.hat, w = res$w.hat)
    ci <- sprintf("[%.1f, %.1f]", res$mu_ipw - 1.96*res$se_ipw, res$mu_ipw + 1.96*res$se_ipw)
    w_str <- paste(round(res$w.hat, 4), collapse = ", ")
    cat(sprintf("  %-15s %12.2f %12.2f %20s %30s\n",
        model_set_names[j], res$mu_ipw, res$se_ipw, ci, w_str))
  }
  cat("  ", paste(rep("-", 95), collapse = ""), "\n")

  # Bootstrap
  cat("\n=== Standard Bootstrap (B=1000) ===\n\n")
  B <- 1000
  n_cores <- min(detectCores() - 1, 10)
  cat(sprintf("B=%d, cores=%d\n", B, n_cores))

  boot_file <- sprintf("Results/boot_Fast4_%s_B1000.RDS", gsub("[: +]", "_", arm_name))

  if (file.exists(boot_file)) {
    cat("Loading existing bootstrap results...\n")
    boot_mat <- readRDS(boot_file)
  } else {
    cl <- makeCluster(n_cores)
    registerDoSNOW(cl)
    clusterExport(cl, c("arm_data", "ps_specifications", "W", "h_nu",
                         "all_model_sets", "model_set_labels"), envir = environment())
    clusterEvalQ(cl, devtools::load_all("../EPS", quiet = TRUE))

    pb <- txtProgressBar(max = B, style = 3)
    progress <- function(nn) setTxtProgressBar(pb, nn)
    opts <- list(progress = progress)

    boot_mat <- foreach(
      b = 1:B, .combine = 'rbind', .options.snow = opts,
      .packages = c("stringr", "Matrix", "numDeriv")
    ) %dopar% {
      set.seed(12345 + b)
      n_arm <- nrow(arm_data)
      idx <- sample(1:n_arm, n_arm, replace = TRUE)
      dat_b <- arm_data[idx, ]
      mu_cc_b <- mean(dat_b$y[dat_b$r == 1])

      mu_vec <- setNames(rep(NA, 7), model_set_labels)
      tryCatch({
        ebmr_b <- EPS$new("y", ps_specifications, dat_b, W)
        for (jj in 1:7) {
          tryCatch({
            res_b <- ebmr_b$EBMR_IPW(h_nu = h_nu, model_indices = all_model_sets[[jj]],
                                       true_ps = NULL, se.fit = FALSE)
            mu_vec[jj] <- res_b$mu_ipw
          }, error = function(e) {})
        }
      }, error = function(e) {})
      c(CC = mu_cc_b, mu_vec)
    }

    close(pb)
    stopCluster(cl)

    dir.create("Results", showWarnings = FALSE)
    saveRDS(boot_mat, boot_file)
    cat("\nSaved:", boot_file, "\n")
  }

  # Summary table
  cat("\n=== Final Summary ===\n\n")
  cat(sprintf("  %-15s %12s %12s %12s %12s %8s\n",
              "Models", "Estimate", "Analytical", "Boot(raw)", "Boot(trim)", "Out"))
  cat("  ", paste(rep("-", 72), collapse = ""), "\n")

  # CC
  bt_cc <- compute_trimmed_boot_se(boot_mat[, "CC"])
  cat(sprintf("  %-15s %12.2f %12s %12.2f %12.2f %8d\n",
              "Naive (CC)", naive_est, "-", bt_cc$se_raw, bt_cc$se_trim, bt_cc$n_outlier))

  for (j in 1:7) {
    bt <- compute_trimmed_boot_se(boot_mat[, model_set_labels[j]])
    n_na <- sum(is.na(boot_mat[, model_set_labels[j]]))
    cat(sprintf("  %-15s %12.2f %12.2f %12.2f %12.2f %8d  (NA:%d)\n",
        model_set_names[j], results_list[[j]]$mu, results_list[[j]]$se,
        bt$se_raw, bt$se_trim, bt$n_outlier, n_na))
  }
  cat("  ", paste(rep("-", 72), collapse = ""), "\n")

  return(results_list)
}

#------------------------------------------------------------------------------#
# Run analysis for each treatment arm
#------------------------------------------------------------------------------#
arm_names <- c("Arm 0: ZDV only", "Arm 1: ZDV + ddI",
               "Arm 2: ZDV + ddC", "Arm 3: ddI only")

results_by_arm <- list()
for (arm in 0:3) {
  arm_data <- dat[dat$arms == arm, ]
  results_by_arm[[arm + 1]] <- run_arm_analysis(arm_data, arm_names[arm + 1])
}

#------------------------------------------------------------------------------#
# Summary comparison across arms
#------------------------------------------------------------------------------#
cat("\n")
cat(strrep("=", 80), "\n")
cat("SUMMARY: COMPARISON ACROSS TREATMENT ARMS (M1+M2+M3 ensemble)\n")
cat(strrep("=", 80), "\n\n")

cat(sprintf("  %-20s %8s %12s %12s %20s\n",
            "Treatment Arm", "n", "Estimate", "SE", "95% CI"))
cat("  ", paste(rep("-", 76), collapse = ""), "\n")

for (arm in 0:3) {
  arm_data <- dat[dat$arms == arm, ]
  res <- results_by_arm[[arm + 1]][[7]]  # Models 1+2+3
  if (!is.na(res$mu)) {
    ci <- sprintf("[%.1f, %.1f]", res$mu - 1.96*res$se, res$mu + 1.96*res$se)
    cat(sprintf("  %-20s %8d %12.2f %12.2f %20s\n",
        arm_names[arm + 1], nrow(arm_data), res$mu, res$se, ci))
  } else {
    cat(sprintf("  %-20s %8d %12s %12s %20s\n",
        arm_names[arm + 1], nrow(arm_data), "ERROR", "-", "-"))
  }
}
cat("  ", paste(rep("-", 76), collapse = ""), "\n")

# Naive comparison
cat("\nNaive estimates (observed mean):\n")
cat(sprintf("  %-20s %8s %12s\n", "Treatment Arm", "n(obs)", "Mean cd496"))
cat("  ", paste(rep("-", 44), collapse = ""), "\n")
for (arm in 0:3) {
  arm_data <- dat[dat$arms == arm, ]
  cat(sprintf("  %-20s %8d %12.2f\n",
      arm_names[arm + 1], sum(arm_data$r),
      mean(arm_data$y[arm_data$r == 1])))
}

cat("\nNote: cd496 = CD4 count at 96 weeks (cells/mm^3)\n")
cat("Higher values indicate better immune function.\n")
cat("\n================================================================================\n")
