## Diagnose: Why is alpha unstable across bootstrap reps?
## Test different optimization strategies for stability
## Focus on the 111 ensemble
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
h_nu <- function(data) cbind(health=data$health, father=data$father, parent_report=data$parent_report,
                              fp=data$fp, fh=data$fh, hp=data$hp, fhp=data$fhp)

ps_spec_full <- list(
  formula.list = list(
    r ~ teacher_report + health + father + teacher_report:health + teacher_report:father,
    r ~ teacher_report + health + parent_report + teacher_report:health + teacher_report:parent_report,
    r ~ teacher_report + parent_report + father + teacher_report:parent_report + teacher_report:father
  ),
  h_alpha.list = list(h_alpha, h_alpha, h_alpha),
  outcome = "teacher_report",
  inv_link = function(eta) 1 / (1 + exp(eta))
)

B <- 200
cat("=== Testing optimizer strategies on 111 ensemble (B=200) ===\n\n")

# Strategy 1: Fast4 L-BFGS-B (default)
# Strategy 2: Fast4 constrained_nr
# Strategy 3: Fast2 Gauss-Newton (simulate by using Fast2 package)
# Strategy 4: Fast4 constrained_nr with smaller trust_radius
# Strategy 5: Two-step GMM (limit outer iterations)

strategies <- list(
  list(name = "Fast4 L-BFGS-B", optimizer = "L-BFGS-B"),
  list(name = "Fast4 constrained_nr", optimizer = "constrained_nr")
)

for (strat in strategies) {
  cat(sprintf("--- %s ---\n", strat$name))

  mus <- nu1s <- nu2s <- nu3s <- rep(NA, B)
  alphas_m1 <- matrix(NA, B, 6)
  alphas_m2 <- matrix(NA, B, 6)
  alphas_m3 <- matrix(NA, B, 6)

  for (b in 1:B) {
    set.seed(12345 + b)
    idx <- sample(1:n, n, replace = TRUE)
    dat_b <- dat[idx, ]

    tryCatch({
      ps_spec <- ps_spec_full
      ps_spec$optimizer <- strat$optimizer
      ebmr_b <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat_b, W)
      res_b <- ebmr_b$EBMR_IPW(h_nu = h_nu, true_ps = NULL, se.fit = FALSE)

      mus[b] <- res_b$mu_ipw
      nu1s[b] <- res_b$nu.hat[1]
      nu2s[b] <- res_b$nu.hat[2]
      nu3s[b] <- res_b$nu.hat[3]
      alphas_m1[b, ] <- ebmr_b$ps_fit.list[[1]]$coefficients
      alphas_m2[b, ] <- ebmr_b$ps_fit.list[[2]]$coefficients
      alphas_m3[b, ] <- ebmr_b$ps_fit.list[[3]]$coefficients
    }, error = function(e) {})
  }

  valid <- !is.na(mus)
  cat(sprintf("  Valid: %d/%d\n", sum(valid), B))
  cat(sprintf("  mu: mean=%.4f, sd(bootstrap SE)=%.4f\n", mean(mus, na.rm=TRUE), sd(mus, na.rm=TRUE)))

  # Weight distribution
  w1 <- nu1s^2 / (nu1s^2 + nu2s^2 + nu3s^2)
  cat(sprintf("  w1 (M1 weight): mean=%.4f, sd=%.4f, >0.9: %d/%d\n",
      mean(w1, na.rm=TRUE), sd(w1, na.rm=TRUE),
      sum(w1 > 0.9, na.rm=TRUE), sum(valid)))

  # Alpha stability per model
  for (m in 1:3) {
    alphas <- switch(m, alphas_m1, alphas_m2, alphas_m3)
    alpha_sds <- apply(alphas, 2, sd, na.rm=TRUE)
    cat(sprintf("  M%d alpha sd: %s\n", m, paste(round(alpha_sds, 2), collapse=", ")))
  }

  # Correlation between alpha instability and mu
  alpha3_m1 <- alphas_m1[valid, 3]  # teacher_report:health interaction
  cat(sprintf("  cor(alpha_m1[3], mu): %.3f\n", cor(alpha3_m1, mus[valid], use="complete.obs")))
  cat(sprintf("  cor(w1, mu): %.3f\n\n", cor(w1[valid], mus[valid], use="complete.obs")))
}

# Now test Fast2 for comparison
cat("--- Fast2 Gauss-Newton (reference) ---\n")
tryCatch({
  library(EBMRalgorithmFast2)

  ps_spec_f2 <- list(
    formula.list = ps_spec_full$formula.list,
    h_x_names.list = ps_spec_full$h_alpha.list,
    outcome = "teacher_report",
    inv_link = ps_spec_full$inv_link
  )

  mus_f2 <- nu1s_f2 <- nu2s_f2 <- nu3s_f2 <- rep(NA, B)
  alphas_m1_f2 <- matrix(NA, B, 6)

  for (b in 1:B) {
    set.seed(12345 + b)
    idx <- sample(1:n, n, replace = TRUE)
    dat_b <- dat[idx, ]

    tryCatch({
      ebmr_b <- EBMRalgorithmFast2$new("teacher_report", ps_spec_f2, dat_b, W)
      res_b <- ebmr_b$EBMR_IPW(h_nu = h_nu, se.fit = FALSE)

      mus_f2[b] <- res_b$mu_ipw
      nu1s_f2[b] <- res_b$nu.hat[1]
      nu2s_f2[b] <- res_b$nu.hat[2]
      nu3s_f2[b] <- res_b$nu.hat[3]
      alphas_m1_f2[b, ] <- ebmr_b$ps_fit.list[[1]]$coefficients
    }, error = function(e) {})
  }

  valid_f2 <- !is.na(mus_f2)
  cat(sprintf("  Valid: %d/%d\n", sum(valid_f2), B))
  cat(sprintf("  mu: mean=%.4f, sd(bootstrap SE)=%.4f\n", mean(mus_f2, na.rm=TRUE), sd(mus_f2, na.rm=TRUE)))

  w1_f2 <- nu1s_f2^2 / (nu1s_f2^2 + nu2s_f2^2 + nu3s_f2^2)
  cat(sprintf("  w1 (M1 weight): mean=%.4f, sd=%.4f, >0.9: %d/%d\n",
      mean(w1_f2, na.rm=TRUE), sd(w1_f2, na.rm=TRUE),
      sum(w1_f2 > 0.9, na.rm=TRUE), sum(valid_f2)))

  alpha_sds_f2 <- apply(alphas_m1_f2, 2, sd, na.rm=TRUE)
  cat(sprintf("  M1 alpha sd: %s\n", paste(round(alpha_sds_f2, 2), collapse=", ")))
  cat(sprintf("  cor(w1, mu): %.3f\n", cor(w1_f2[valid_f2], mus_f2[valid_f2], use="complete.obs")))

}, error = function(e) {
  cat(sprintf("  Fast2 not available: %s\n", e$message))
})
