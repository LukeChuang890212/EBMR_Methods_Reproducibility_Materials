## Test: Enrich h_alpha and h_nu with ratio terms X1/X2, X1/X3, X2/X3
## X1=health+1, X2=father+1, X3=parent_report+1 (values {1,2})
## Ratios take values {0.5, 1, 2} — nonlinear, not in the binary span
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

# Add ratio terms
dat$X1 <- dat$health + 1
dat$X2 <- dat$father + 1
dat$X3 <- dat$parent_report + 1
dat$X1_X2 <- dat$X1 / dat$X2
dat$X1_X3 <- dat$X1 / dat$X3
dat$X2_X3 <- dat$X2 / dat$X3
dat$X1_X2_X3 <- dat$X1 / (dat$X2 * dat$X3)

n <- nrow(dat)

cat("=== Ratio enrichment test ===\n\n")
cat("Ratio term distinct values:\n")
cat(sprintf("  X1/X2: {%s}\n", paste(sort(unique(dat$X1_X2)), collapse=", ")))
cat(sprintf("  X1/X3: {%s}\n", paste(sort(unique(dat$X1_X3)), collapse=", ")))
cat(sprintf("  X2/X3: {%s}\n", paste(sort(unique(dat$X2_X3)), collapse=", ")))
cat(sprintf("  X1/(X2*X3): {%s}\n\n", paste(sort(unique(dat$X1_X2_X3)), collapse=", ")))

W <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

ps_specifications <- list(
  formula.list = list(
    r ~ teacher_report + health + father + teacher_report:health + teacher_report:father,
    r ~ teacher_report + health + parent_report + teacher_report:health + teacher_report:parent_report,
    r ~ teacher_report + parent_report + father + teacher_report:parent_report + teacher_report:father
  ),
  # h_alpha: original 6 terms + 3 ratio terms
  h_alpha.list = list(
    function(dat) cbind(health=dat$health, father=dat$father, parent_report=dat$parent_report,
                        fp=dat$fp, fh=dat$fh, hp=dat$hp,
                        X1_X2=dat$X1_X2, X1_X3=dat$X1_X3, X2_X3=dat$X2_X3),
    function(dat) cbind(health=dat$health, father=dat$father, parent_report=dat$parent_report,
                        fp=dat$fp, fh=dat$fh, hp=dat$hp,
                        X1_X2=dat$X1_X2, X1_X3=dat$X1_X3, X2_X3=dat$X2_X3),
    function(dat) cbind(health=dat$health, father=dat$father, parent_report=dat$parent_report,
                        fp=dat$fp, fh=dat$fh, hp=dat$hp,
                        X1_X2=dat$X1_X2, X1_X3=dat$X1_X3, X2_X3=dat$X2_X3)
  ),
  outcome = "teacher_report",
  inv_link = function(eta) 1 / (1 + exp(eta))
)

# h_nu: original 7 terms + 4 ratio terms
h_nu <- function(data) {
  cbind(health = data$health, father = data$father, parent_report = data$parent_report,
        fp = data$fp, fh = data$fh, hp = data$hp, fhp = data$fhp,
        X1_X2 = data$X1_X2, X1_X3 = data$X1_X3, X2_X3 = data$X2_X3,
        X1_X2_X3 = data$X1_X2_X3)
}

# Fit
ebmr <- EBMRAlgorithmFast4$new("teacher_report", ps_specifications, dat, W)

cat("=== Model diagnostics ===\n\n")
for (i in 1:3) {
  pf <- ebmr$ps_fit.list[[i]]
  ps <- pf$fitted.values
  cat(sprintf("Model %d: k=%d, h_dim=%d, overid=%d\n", i,
      ncol(pf$design_matrix), ncol(pf$h_x), ncol(pf$h_x) - ncol(pf$design_matrix)))
  cat(sprintf("  alpha = (%s)\n", paste(round(pf$coefficients, 4), collapse=", ")))
  cat(sprintf("  PS range: [%.6f, %.6f], >0.99: %d, <0.01: %d\n",
      min(ps), max(ps), sum(ps > 0.99), sum(ps < 0.01)))
  cat(sprintf("  IPW mean: %.4f\n\n", mean(dat$r * dat$teacher_report / ps)))
}

# All combinations
all_model_sets <- list(c(1), c(2), c(3), c(1,2), c(1,3), c(2,3), c(1,2,3))
model_set_labels <- c("100", "010", "001", "110", "101", "011", "111")

cat("=== Point estimates ===\n\n")
cat(sprintf("  %-8s %10s %10s %20s\n", "Label", "Estimate", "SE", "w.hat"))
cat("  ", paste(rep("-", 52), collapse=""), "\n")

mu_cc <- mean(dat$teacher_report[dat$r == 1])
se_cc <- sd(dat$teacher_report[dat$r == 1]) / sqrt(sum(dat$r))
cat(sprintf("  %-8s %10.4f %10.4f\n", "CC", mu_cc, se_cc))

point_mu <- point_se <- numeric(7)
for (j in 1:7) {
  res <- ebmr$EBMR_IPW(h_nu = h_nu, model_indices = all_model_sets[[j]], true_ps = NULL)
  point_mu[j] <- res$mu_ipw; point_se[j] <- res$se_ipw
  w_str <- paste(round(res$w.hat, 4), collapse=",")
  cat(sprintf("  %-8s %10.4f %10.4f     (%s)\n",
      model_set_labels[j], res$mu_ipw, res$se_ipw, w_str))
}
cat("  ", paste(rep("-", 52), collapse=""), "\n")

# Perturbation bootstrap
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

# Results
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

cat("\n\n=== Results ===\n\n")
labels <- c("CC", model_set_labels)
estimates <- c(mu_cc, point_mu)
analytical <- c(se_cc, point_se)

cat(sprintf("  %-8s %10s %12s %12s %8s %8s\n",
            "Label", "Estimate", "Analytical", "Ptb(trim)", "Out", "NA"))
cat("  ", paste(rep("-", 62), collapse=""), "\n")
for (jj in 1:8) {
  bt <- compute_se_capped(perturb_mat[, jj])
  n_na <- sum(is.na(perturb_mat[, jj]))
  cat(sprintf("  %-8s %10.4f %12.4f %12.4f %8d %8d\n",
              labels[jj], estimates[jj], analytical[jj], bt$se, bt$n_out, n_na))
}
cat("  ", paste(rep("-", 62), collapse=""), "\n")
