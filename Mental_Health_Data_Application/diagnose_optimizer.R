## Diagnose optimizer instability across bootstrap reps
## Focus on M1 (100) since even single-model has SE discrepancy
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Mental_Health_Data_Application")
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)
library(Matrix); library(numDeriv)
source("MHD_functions.R")

original_data <- read.csv("data_application.csv")
percent <- original_data$Percentage
n <- 2486
class_n <- round(n * percent / 100)
dat <- gen_data(original_data, class_n, n)
dat$y <- dat$teacher_report

W <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

h_alpha <- c("health", "father", "parent_report", "fp", "fh", "hp")

# M1 only
ps_spec_m1 <- list(
  formula.list = list(
    r ~ teacher_report + health + father + teacher_report:health + teacher_report:father
  ),
  h_alpha.list = list(h_alpha),
  outcome = "teacher_report",
  inv_link = function(eta) 1 / (1 + exp(eta))
)

h_nu <- function(data) {
  cbind(health = data$health, father = data$father, parent_report = data$parent_report,
        fp = data$fp, fh = data$fh, hp = data$hp, fhp = data$fhp)
}

cat("=== Diagnosing optimizer instability ===\n\n")

# 1. Point estimate with both optimizers
cat("--- Point estimate comparison ---\n\n")

for (opt in c("L-BFGS-B", "constrained_nr")) {
  ps_spec <- ps_spec_m1
  ps_spec$optimizer <- opt
  ebmr <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat, W)
  pf <- ebmr$ps_fit.list[[1]]
  res <- ebmr$EBMR_IPW(h_nu = h_nu, true_ps = NULL)
  cat(sprintf("  %s:\n", opt))
  cat(sprintf("    alpha = (%s)\n", paste(round(pf$coefficients, 6), collapse=", ")))
  cat(sprintf("    GMM obj = %.6e, grad = %.6e\n", pf$gmm_fit$opt$objective, pf$gmm_fit$opt$final_grad_norm))
  cat(sprintf("    mu = %.6f, se = %.6f\n", res$mu_ipw, res$se_ipw))
  cat(sprintf("    PS range = [%.6f, %.6f]\n\n", min(pf$fitted.values), max(pf$fitted.values)))
}

# 2. Bootstrap: track alpha, objective, grad for each rep
cat("--- Bootstrap stability analysis (200 reps) ---\n\n")

B <- 200
boot_results <- list()

for (opt in c("L-BFGS-B", "constrained_nr")) {
  cat(sprintf("Optimizer: %s\n", opt))

  alphas <- matrix(NA, B, 6)
  mus <- ses <- objs <- grads <- rep(NA, B)

  for (b in 1:B) {
    set.seed(12345 + b)
    idx <- sample(1:n, n, replace = TRUE)
    dat_b <- dat[idx, ]

    tryCatch({
      ps_spec <- ps_spec_m1
      ps_spec$optimizer <- opt
      ebmr_b <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat_b, W)
      pf <- ebmr_b$ps_fit.list[[1]]
      res <- ebmr_b$EBMR_IPW(h_nu = h_nu, true_ps = NULL, se.fit = FALSE)

      alphas[b, ] <- pf$coefficients
      mus[b] <- res$mu_ipw
      objs[b] <- pf$gmm_fit$opt$objective
      grads[b] <- pf$gmm_fit$opt$final_grad_norm
    }, error = function(e) {})
  }

  valid <- !is.na(mus)
  cat(sprintf("  Valid: %d/%d, NA: %d\n", sum(valid), B, sum(!valid)))
  cat(sprintf("  mu: mean=%.4f, sd=%.4f\n", mean(mus, na.rm=TRUE), sd(mus, na.rm=TRUE)))
  cat(sprintf("  Bootstrap SE of mu: %.4f\n", sd(mus, na.rm=TRUE)))

  # Alpha variability
  cat("  Alpha variability (sd across bootstrap reps):\n")
  alpha_sds <- apply(alphas, 2, sd, na.rm=TRUE)
  cat(sprintf("    %s\n", paste(round(alpha_sds, 4), collapse=", ")))

  # Convergence quality
  cat(sprintf("  GMM obj: mean=%.2e, max=%.2e\n", mean(objs, na.rm=TRUE), max(objs, na.rm=TRUE)))
  cat(sprintf("  Grad norm: mean=%.2e, max=%.2e\n", mean(grads, na.rm=TRUE), max(grads, na.rm=TRUE)))
  cat(sprintf("  Grad > 1e-4: %d reps\n", sum(grads > 1e-4, na.rm=TRUE)))
  cat(sprintf("  Grad > 1e-3: %d reps\n", sum(grads > 1e-3, na.rm=TRUE)))
  cat(sprintf("  Grad > 1e-2: %d reps\n\n", sum(grads > 1e-2, na.rm=TRUE)))

  boot_results[[opt]] <- list(alphas=alphas, mus=mus, objs=objs, grads=grads)
}

# 3. Check if poorly converged reps drive the excess variability
cat("--- Effect of filtering by convergence quality ---\n\n")

for (opt in names(boot_results)) {
  br <- boot_results[[opt]]
  for (tol in c(1e-2, 1e-3, 1e-4, 1e-5, 1e-6)) {
    good <- !is.na(br$grads) & br$grads < tol
    if (sum(good) > 10) {
      cat(sprintf("  %s, grad < %.0e: n=%d, mu_sd=%.4f\n",
          opt, tol, sum(good), sd(br$mus[good])))
    }
  }
  cat("\n")
}

# 4. Check: are there multiple distinct solutions?
cat("--- Checking for distinct alpha clusters ---\n\n")
for (opt in names(boot_results)) {
  br <- boot_results[[opt]]
  valid_alpha <- br$alphas[!is.na(br$mus), ]

  # Check alpha[2] (teacher_report coefficient) distribution
  a2 <- valid_alpha[, 2]
  cat(sprintf("  %s: alpha_teacher_report: mean=%.4f, sd=%.4f, min=%.4f, max=%.4f\n",
      opt, mean(a2), sd(a2), min(a2), max(a2)))
  cat(sprintf("    Q5=%.4f, Q25=%.4f, med=%.4f, Q75=%.4f, Q95=%.4f\n",
      quantile(a2, 0.05), quantile(a2, 0.25), median(a2),
      quantile(a2, 0.75), quantile(a2, 0.95)))
}
