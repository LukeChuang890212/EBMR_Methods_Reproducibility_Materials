## Test: Recode binary variables from {0,1} to {-1,1}
## teacher_report: Normal(0)->-1, Abnormal(1)->1, Missing(-1)->0
## health, father, parent_report: 0->-1, 1->1
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

# Save original for comparison
dat_orig <- dat

# Recode: teacher_report {-1,0,1} -> {0,-1,1}
# Currently: -1=missing, 0=Normal, 1=Abnormal
# New:        0=missing, -1=Normal, 1=Abnormal
dat$teacher_report <- ifelse(dat$teacher_report == -1, 0,
                      ifelse(dat$teacher_report == 0, -1, 1))
dat$y <- dat$teacher_report

# Recode other binary variables: {0,1} -> {-1,1}
dat$health <- 2 * dat$health - 1
dat$father <- 2 * dat$father - 1
dat$parent_report <- 2 * dat$parent_report - 1

# Recompute interactions with recoded variables
dat$fp <- dat$father * dat$parent_report
dat$fh <- dat$father * dat$health
dat$hp <- dat$health * dat$parent_report
dat$fhp <- dat$father * dat$health * dat$parent_report

n <- nrow(dat)

cat("=== Recoding test: {0,1} -> {-1,1} ===\n\n")
cat("Variable ranges after recoding:\n")
cat(sprintf("  teacher_report: {%s}\n", paste(sort(unique(dat$teacher_report)), collapse=", ")))
cat(sprintf("  health: {%s}\n", paste(sort(unique(dat$health)), collapse=", ")))
cat(sprintf("  father: {%s}\n", paste(sort(unique(dat$father)), collapse=", ")))
cat(sprintf("  parent_report: {%s}\n", paste(sort(unique(dat$parent_report)), collapse=", ")))
cat(sprintf("  r: {%s}\n\n", paste(sort(unique(dat$r)), collapse=", ")))

# Note: mu_cc uses original teacher_report values for the mean
# With {-1,1} coding, the population mean target changes
# E[Y] where Y in {-1,1} = 2*P(Abnormal) - 1
mu_cc_orig <- mean(dat_orig$teacher_report[dat_orig$r == 1])
mu_cc_recoded <- mean(dat$teacher_report[dat$r == 1])
cat(sprintf("Complete case mean (original {0,1}): %.4f\n", mu_cc_orig))
cat(sprintf("Complete case mean (recoded {-1,1}): %.4f\n", mu_cc_recoded))
cat(sprintf("Relationship: recoded = 2*original - 1 = %.4f\n\n", 2*mu_cc_orig - 1))

W <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

ps_specifications <- list(
  formula.list = list(
    r ~ teacher_report + health + father + teacher_report:health + teacher_report:father,
    r ~ teacher_report + health + parent_report + teacher_report:health + teacher_report:parent_report,
    r ~ teacher_report + parent_report + father + teacher_report:parent_report + teacher_report:father
  ),
  h_alpha.list = list(
    c("health", "father", "parent_report", "fp", "fh", "hp"),
    c("health", "father", "parent_report", "fp", "fh", "hp"),
    c("health", "father", "parent_report", "fp", "fh", "hp")
  ),
  outcome = "teacher_report",
  inv_link = function(eta) 1 / (1 + exp(eta))
)

h_nu <- function(data) {
  cbind(health = data$health, father = data$father, parent_report = data$parent_report,
        fp = data$fp, fh = data$fh, hp = data$hp, fhp = data$fhp)
}

# Fit and diagnose
ebmr <- EBMRAlgorithmFast4$new("teacher_report", ps_specifications, dat, W)

cat("=== Model diagnostics ===\n\n")
for (i in 1:3) {
  pf <- ebmr$ps_fit.list[[i]]
  ps <- pf$fitted.values
  cat(sprintf("Model %d:\n", i))
  cat(sprintf("  alpha = (%s)\n", paste(round(pf$coefficients, 4), collapse=", ")))
  cat(sprintf("  SE    = (%s)\n", paste(round(pf$se, 4), collapse=", ")))
  cat(sprintf("  PS range: [%.6f, %.6f], >0.99: %d, <0.01: %d\n",
      min(ps), max(ps), sum(ps > 0.99), sum(ps < 0.01)))
  mu_single <- mean(dat$r * dat$teacher_report / ps)
  cat(sprintf("  IPW mean (recoded): %.4f, (back to original): %.4f\n\n",
      mu_single, (mu_single + 1) / 2))
}

# All model combinations
all_model_sets <- list(c(1), c(2), c(3), c(1,2), c(1,3), c(2,3), c(1,2,3))
model_set_labels <- c("100", "010", "001", "110", "101", "011", "111")

cat("=== All model combinations ===\n\n")
cat(sprintf("  %-8s %12s %12s %12s %12s\n",
            "Label", "mu(recoded)", "mu(original)", "SE", "w.hat"))
cat("  ", paste(rep("-", 62), collapse = ""), "\n")

results_list <- list()
for (j in 1:7) {
  res <- ebmr$EBMR_IPW(h_nu = h_nu, model_indices = all_model_sets[[j]], true_ps = NULL)
  mu_orig_scale <- (res$mu_ipw + 1) / 2
  se_orig_scale <- res$se_ipw / 2
  w_str <- paste(round(res$w.hat, 4), collapse=",")
  cat(sprintf("  %-8s %12.4f %12.4f %12.4f     (%s)\n",
      model_set_labels[j], res$mu_ipw, mu_orig_scale, se_orig_scale, w_str))
  results_list[[j]] <- list(mu = res$mu_ipw, se = res$se_ipw,
                             mu_orig = mu_orig_scale, se_orig = se_orig_scale)
}
cat("  ", paste(rep("-", 62), collapse = ""), "\n")

# Bootstrap comparison
cat("\n=== Perturbation bootstrap (B=1000) ===\n\n")

B <- 1000
n_cores <- min(detectCores() - 1, 10)

cl <- makeCluster(n_cores)
registerDoSNOW(cl)
clusterExport(cl, c("dat", "n", "ps_specifications", "all_model_sets",
                     "model_set_labels", "W", "h_nu"), envir = environment())
clusterEvalQ(cl, devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE))

pb <- txtProgressBar(max = B, style = 3)
progress <- function(nn) setTxtProgressBar(pb, nn)
opts <- list(progress = progress)

perturb_mat <- foreach(
  b = 1:B, .combine = 'rbind', .options.snow = opts,
  .packages = c("stringr", "Matrix", "numDeriv")
) %dopar% {
  set.seed(12345 + b)
  wt <- rexp(n, rate = 1)
  r_idx <- dat$r == 1
  mu_cc_b <- sum(wt[r_idx] * dat$teacher_report[r_idx]) / sum(wt[r_idx])

  mu_vec <- setNames(rep(NA, 7), model_set_labels)
  for (j in seq_along(all_model_sets)) {
    ms <- all_model_sets[[j]]
    ps_b <- list(formula.list = ps_specifications$formula.list[ms],
                 h_alpha.list = ps_specifications$h_alpha.list[ms],
                 outcome = ps_specifications$outcome,
                 inv_link = ps_specifications$inv_link)
    mu_vec[j] <- tryCatch({
      eb <- EBMRAlgorithmFast4$new("teacher_report", ps_b, dat, W, wt = wt)
      rb <- eb$EBMR_IPW(h_nu = h_nu, true_ps = NULL, se.fit = FALSE, wt = wt)
      rb$mu_ipw
    }, error = function(e) NA)
  }
  c(CC = mu_cc_b, mu_vec)
}
close(pb)
stopCluster(cl)

# Results on original scale
cat("\n\n=== Results (on original {0,1} scale) ===\n\n")

compute_se_capped <- function(x) {
  x <- x[!is.na(x)]; n_tot <- length(x)
  if (n_tot < 10) return(list(se = NA, n_out = NA))
  Q1 <- quantile(x, 0.25); Q3 <- quantile(x, 0.75); IQR_val <- Q3 - Q1
  is_out <- x < (Q1 - 3 * IQR_val) | x > (Q3 + 3 * IQR_val)
  max_rm <- floor(0.01 * n_tot)
  if (sum(is_out) > max_rm && max_rm > 0) {
    dm <- abs(x - median(x)); oi <- which(is_out)
    keep <- oi[order(dm[oi], decreasing = TRUE)[(max_rm + 1):length(oi)]]
    is_out[keep] <- FALSE
  }
  list(se = sd(x[!is_out]), n_out = sum(is_out))
}

labels <- c("CC", model_set_labels)
# Convert back to original scale: mu_orig = (mu_recoded + 1)/2, se_orig = se_recoded/2
estimates_orig <- c(mu_cc_orig,
                    sapply(results_list, function(r) r$mu_orig))
analytical_orig <- c(sd(dat_orig$teacher_report[dat_orig$r == 1]) / sqrt(sum(dat_orig$r)),
                     sapply(results_list, function(r) r$se_orig))

cat(sprintf("  %-8s %10s %12s %12s %8s\n",
            "Label", "Estimate", "Analytical", "Ptb(trim)", "Out"))
cat("  ", paste(rep("-", 56), collapse = ""), "\n")

for (jj in 1:8) {
  # Convert perturbation results to original scale
  ptb_orig <- (perturb_mat[, jj] + 1) / 2
  bt <- compute_se_capped(ptb_orig)
  cat(sprintf("  %-8s %10.4f %12.4f %12.4f %8d\n",
              labels[jj], estimates_orig[jj], analytical_orig[jj], bt$se, bt$n_out))
}
cat("  ", paste(rep("-", 56), collapse = ""), "\n")
