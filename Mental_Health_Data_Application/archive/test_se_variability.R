## Examine variability of analytical SE across bootstrap reps
## Compare analytical SE distribution with bootstrap SE
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
h_alpha <- c("health", "father", "parent_report", "fp", "fh", "hp", "fhp")
h_nu <- function(data) cbind(health=data$health, father=data$father, parent_report=data$parent_report,
                              fp=data$fp, fh=data$fh, hp=data$hp, fhp=data$fhp)

ps_spec <- list(
  formula.list = list(
    r ~ teacher_report + health + father + teacher_report:health + teacher_report:father + father:health,
    r ~ teacher_report + parent_report + health + teacher_report:parent_report + teacher_report:health + parent_report:health,
    r ~ teacher_report + parent_report + father + teacher_report:parent_report + teacher_report:father + parent_report:father
  ),
  h_alpha.list = list(h_alpha, h_alpha, h_alpha),
  outcome = "teacher_report", inv_link = function(eta) 1/(1+exp(eta)), optimizer = "L-BFGS-B"
)

# Point estimates
ebmr <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat, W)
alpha_init <- lapply(ebmr$ps_fit.list, function(pf) pf$coefficients)
res_pt <- ebmr$EBMR_IPW(h_nu = h_nu, true_ps = NULL)
nu_init <- res_pt$nu.hat

all_sets <- list(c(1), c(2), c(3), c(1,2), c(1,3), c(2,3), c(1,2,3))
labels <- c("100", "010", "001", "110", "101", "011", "111")

# Get point analytical SE for each estimator
point_se <- numeric(7)
for (j in 1:7) {
  res <- ebmr$EBMR_IPW(h_nu = h_nu, model_indices = all_sets[[j]], true_ps = NULL)
  point_se[j] <- res$se_ipw
}

B <- 1000
cat(sprintf("=== SE variability across bootstrap (B=%d) ===\n\n", B))

library(parallel); library(foreach); library(doSNOW)
n_cores <- min(detectCores() - 1, 10)
cl <- makeCluster(n_cores)
registerDoSNOW(cl)
clusterExport(cl, c("dat", "n", "ps_spec", "W", "h_nu", "all_sets", "labels",
                     "alpha_init", "nu_init"), envir = environment())
clusterEvalQ(cl, devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE))

pb <- txtProgressBar(max = B, style = 3)
progress <- function(nn) setTxtProgressBar(pb, nn)
opts <- list(progress = progress)

boot_raw <- foreach(b = 1:B, .combine = 'rbind', .options.snow = opts,
                     .packages = c("stringr", "Matrix", "numDeriv")) %dopar% {
  set.seed(12345 + b)
  idx <- sample(1:n, n, replace = TRUE)
  dat_b <- dat[idx, ]

  row <- rep(NA, 14)  # 7 mu + 7 se
  tryCatch({
    ps_spec_b <- list(
      formula.list = ps_spec$formula.list,
      h_alpha.list = ps_spec$h_alpha.list,
      alpha_init.list = alpha_init,
      outcome = ps_spec$outcome,
      inv_link = ps_spec$inv_link
    )
    ebmr_b <- EBMRAlgorithmFast4$new("teacher_report", ps_spec_b, dat_b, W)

    for (j in 1:7) {
      tryCatch({
        nu_init_b <- nu_init[all_sets[[j]]]
        res_b <- ebmr_b$EBMR_IPW(h_nu = h_nu, model_indices = all_sets[[j]],
                                   nu_init = nu_init_b, true_ps = NULL, se.fit = TRUE)
        row[j] <- res_b$mu_ipw
        row[j + 7] <- res_b$se_ipw
      }, error = function(e) {})
    }
  }, error = function(e) {})
  row
}

close(pb)
stopCluster(cl)

mu_mat <- boot_raw[, 1:7]
se_mat <- boot_raw[, 8:14]

cat("\n=== Results ===\n\n")
cat(sprintf("%-8s %10s %10s %10s %10s %10s %10s\n",
            "Label", "Point SE", "Mean SE_b", "Med SE_b", "SD SE_b", "ESD mu_b", "Ratio"))
cat(paste(rep("-", 72), collapse=""), "\n")

for (j in 1:7) {
  valid <- !is.na(mu_mat[, j]) & !is.na(se_mat[, j])
  mu_v <- mu_mat[valid, j]
  se_v <- se_mat[valid, j]
  esd <- sd(mu_v)
  cat(sprintf("%-8s %10.4f %10.4f %10.4f %10.4f %10.4f %10.3f\n",
      labels[j], point_se[j], mean(se_v), median(se_v), sd(se_v), esd,
      mean(se_v) / esd))
}

cat("\n=== SE distribution quantiles ===\n\n")
cat(sprintf("%-8s %10s %10s %10s %10s %10s %10s %10s\n",
            "Label", "Min", "Q5", "Q25", "Median", "Q75", "Q95", "Max"))
cat(paste(rep("-", 82), collapse=""), "\n")

for (j in 1:7) {
  valid <- !is.na(se_mat[, j])
  se_v <- se_mat[valid, j]
  cat(sprintf("%-8s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",
      labels[j], min(se_v), quantile(se_v, 0.05), quantile(se_v, 0.25),
      median(se_v), quantile(se_v, 0.75), quantile(se_v, 0.95), max(se_v)))
}
