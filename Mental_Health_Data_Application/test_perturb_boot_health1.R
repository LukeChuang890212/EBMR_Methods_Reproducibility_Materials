## Test: Perturbation bootstrap for health=1 analysis
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Mental_Health_Data_Application")
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)
library(Matrix); library(numDeriv)
library(parallel); library(foreach); library(doSNOW)
source("MHD_functions.R")

original_data <- read.csv("data_application.csv")
percent <- original_data$Percentage
n_full <- 2486
class_n <- round(n_full * percent / 100)
dat_full <- gen_data(original_data, class_n, n_full)
dat_full$y <- dat_full$teacher_report
dat <- dat_full[dat_full$health == 1, ]
n <- nrow(dat)

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

all_model_sets <- list(c(1), c(2), c(3), c(1,2), c(1,3), c(2,3), c(1,2,3))
model_set_labels <- c("100", "010", "001", "110", "101", "011", "111")

# Point estimates
ebmr <- EBMRAlgorithmFast4$new("teacher_report", ps_specifications, dat, W)
point_mu <- point_se <- numeric(7)
for (j in 1:7) {
  res <- ebmr$EBMR_IPW(h_nu = h_nu, model_indices = all_model_sets[[j]], true_ps = NULL)
  point_mu[j] <- res$mu_ipw
  point_se[j] <- res$se_ipw
}

cat("=== Perturbation bootstrap for health=1 ===\n\n")

B <- 1000
n_cores <- min(detectCores() - 1, 10)
cat(sprintf("B=%d, cores=%d\n\n", B, n_cores))

cl <- makeCluster(n_cores)
registerDoSNOW(cl)
clusterExport(cl, c("dat", "n", "ps_specifications", "all_model_sets",
                     "model_set_labels", "W", "h_nu"),
              envir = environment())
clusterEvalQ(cl, devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE))

pb <- txtProgressBar(max = B, style = 3)
progress <- function(nn) setTxtProgressBar(pb, nn)
opts <- list(progress = progress)

perturb_mat <- foreach(
  b = 1:B,
  .combine = 'rbind',
  .options.snow = opts,
  .packages = c("stringr", "Matrix", "numDeriv")
) %dopar% {
  set.seed(12345 + b)
  wt <- rexp(n, rate = 1)

  # Weighted complete case
  r_idx <- dat$r == 1
  mu_cc_b <- sum(wt[r_idx] * dat$teacher_report[r_idx]) / sum(wt[r_idx])

  mu_vec <- setNames(rep(NA, 7), model_set_labels)

  for (j in seq_along(all_model_sets)) {
    model_set <- all_model_sets[[j]]
    label <- model_set_labels[j]

    ps_spec_b <- list(
      formula.list = ps_specifications$formula.list[model_set],
      h_alpha.list = ps_specifications$h_alpha.list[model_set],
      outcome = ps_specifications$outcome,
      inv_link = ps_specifications$inv_link
    )

    mu_vec[j] <- tryCatch({
      ebmr_b <- EBMRAlgorithmFast4$new("teacher_report", ps_spec_b, dat, W, wt = wt)
      result_b <- ebmr_b$EBMR_IPW(h_nu = h_nu, true_ps = NULL, se.fit = FALSE, wt = wt)
      result_b$mu_ipw
    }, error = function(e) NA)
  }

  c(CC = mu_cc_b, mu_vec)
}

close(pb)
stopCluster(cl)

# Outlier removal functions
compute_trimmed_se_capped <- function(x) {
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

compute_trimmed_se_nocap <- function(x) {
  x <- x[!is.na(x)]; n_tot <- length(x)
  if (n_tot < 10) return(list(se = NA, n_out = NA))
  Q1 <- quantile(x, 0.25); Q3 <- quantile(x, 0.75); IQR_val <- Q3 - Q1
  is_out <- x < (Q1 - 3 * IQR_val) | x > (Q3 + 3 * IQR_val)
  list(se = sd(x[!is_out]), n_out = sum(is_out))
}

# Results
cat("\n\n=== Results ===\n\n")
labels <- c("CC", model_set_labels)
estimates <- c(mean(dat$teacher_report[dat$r == 1]), point_mu)
analytical <- c(sd(dat$teacher_report[dat$r == 1]) / sqrt(sum(dat$r)), point_se)

cat(sprintf("  %-8s %10s %12s %12s %8s %12s %8s\n",
            "Label", "Estimate", "Analytical", "Trim(1%%cap)", "Out", "Trim(nocap)", "Out"))
cat("  ", paste(rep("-", 78), collapse = ""), "\n")

for (j in 1:8) {
  bt1 <- compute_trimmed_se_capped(perturb_mat[, j])
  bt2 <- compute_trimmed_se_nocap(perturb_mat[, j])
  cat(sprintf("  %-8s %10.4f %12.4f %12.4f %8d %12.4f %8d\n",
              labels[j], estimates[j], analytical[j],
              bt1$se, bt1$n_out, bt2$se, bt2$n_out))
}
cat("  ", paste(rep("-", 78), collapse = ""), "\n")
