setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
source("Simulation.r")
library(EBMRalgorithmFast4)

data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
n_val <- 2000

ps_spec <- get_ps_spec("9-alt1")
inv_link_func <- ps_spec[["inv_link"]]

# Model 3 (misspecified): r ~ y + u2 + z2
formula_m3 <- ps_spec[["formula.list"]][[3]]
# Model 1 (correct): r ~ y + u1 + u2
formula_m1 <- ps_spec[["formula.list"]][[1]]

test_reps <- c(1, 2, 3, 4, 5, 8, 28, 92, 133, 155)

cat("=== COMPARING NORMAL vs OUTLIER REPS ===\n\n")
cat(sprintf("%-5s %-8s %8s %8s %8s %8s %8s %8s %8s %8s %8s\n",
    "Rep", "Type", "miss%", "cor(y,r)", "cor(u2,r)", "cor(z2,r)",
    "sd(y)", "sd(u2)", "sd(z2)", "mean(y|r=1)", "mean(y|r=0)"))

for (rep_i in test_reps) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]
  is_outlier <- rep_i %in% c(8, 28, 83, 92, 133, 155)

  r <- dat[["r"]]
  y <- dat[["y"]]
  u1 <- dat[["u1"]]
  u2 <- dat[["u2"]]
  z1 <- dat[["z1"]]
  z2 <- dat[["z2"]]

  cat(sprintf("%-5d %-8s %7.1f%% %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",
      rep_i, ifelse(is_outlier, "OUTLIER", "normal"),
      100 * (1 - mean(r)),
      cor(y, r), cor(u2, r), cor(z2, r),
      sd(y), sd(u2), sd(z2),
      mean(y[r == 1]), mean(y[r == 0])))
}

cat("\n=== DETAILED COMPARISON: DATA CHARACTERISTICS ===\n\n")
for (rep_i in test_reps) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]
  is_outlier <- rep_i %in% c(8, 28, 83, 92, 133, 155)

  r <- dat[["r"]]
  y <- dat[["y"]]
  u2 <- dat[["u2"]]
  z2 <- dat[["z2"]]
  n <- nrow(dat)

  # Design matrix for model 3: r ~ y + u2 + z2
  X3 <- model.matrix(formula_m3, dat)
  # Design matrix for model 1: r ~ y + u1 + u2
  X1 <- model.matrix(formula_m1, dat)

  # h_x (moment conditions)
  h_x <- cbind(1, u1 = dat[["u1"]], u2 = u2, z1 = dat[["z1"]], z2 = z2)

  # Check: for the correctly specified model, E[(r/pi - 1)*h_x] should be ~0
  # For misspecified, it won't be exactly 0, but how far off?

  # Fit simple logistic regression for each model to get baseline PS
  # Model 1 (correct): r ~ y + u1 + u2
  glm1 <- glm(r ~ y + dat[["u1"]] + u2, family = binomial(link = "logit"))
  ps1 <- 1 / (1 + exp(predict(glm1)))  # logistic complement

  # Model 3 (misspecified): r ~ y + u2 + z2
  glm3 <- glm(r ~ y + u2 + z2, family = binomial(link = "logit"))
  ps3 <- 1 / (1 + exp(predict(glm3)))  # logistic complement

  # Moment conditions with MLE-based PS
  g1_mle <- (r / ps1 - 1) * h_x
  g3_mle <- (r / ps3 - 1) * h_x
  gbar1 <- colMeans(g1_mle)
  gbar3 <- colMeans(g3_mle)

  cat(sprintf("\n--- Rep %d [%s] ---\n", rep_i, ifelse(is_outlier, "OUTLIER", "normal")))
  cat(sprintf("  Model 1 MLE |g_bar|: %s  ||g_bar||=%.6f\n",
      paste(sprintf("%.4f", gbar1), collapse=", "), sqrt(sum(gbar1^2))))
  cat(sprintf("  Model 3 MLE |g_bar|: %s  ||g_bar||=%.6f\n",
      paste(sprintf("%.4f", gbar3), collapse=", "), sqrt(sum(gbar3^2))))

  # Condition number of X'X for each model
  cat(sprintf("  cond(X1'X1)=%.1f, cond(X3'X3)=%.1f\n",
      kappa(crossprod(X1)), kappa(crossprod(X3))))

  # How much do moment conditions improve when going from MLE to GMM?
  # The GMM can reduce g_bar by finding extreme alpha — check the "room" for improvement
  cat(sprintf("  Model 3: max|g_bar_mle|=%.6f, sum|g_bar_mle|=%.6f\n",
      max(abs(gbar3)), sum(abs(gbar3))))

  # Key: correlation between IV (z1, z2) and outcome y among respondents
  r1_idx <- r == 1
  cat(sprintf("  Among r=1: cor(y,z1)=%.4f, cor(y,z2)=%.4f, cor(u1,z2)=%.4f, cor(u2,z2)=%.4f\n",
      cor(y[r1_idx], dat[["z1"]][r1_idx]),
      cor(y[r1_idx], z2[r1_idx]),
      cor(dat[["u1"]][r1_idx], z2[r1_idx]),
      cor(u2[r1_idx], z2[r1_idx])))

  # Check separation: can y + u2 + z2 perfectly predict r?
  cat(sprintf("  GLM3 coefs: %s\n",
      paste(sprintf("%.4f", coef(glm3)), collapse=", ")))
  cat(sprintf("  GLM3 deviance=%.1f, null_dev=%.1f, ratio=%.4f\n",
      deviance(glm3), glm3[["null.deviance"]], deviance(glm3)/glm3[["null.deviance"]]))
}
