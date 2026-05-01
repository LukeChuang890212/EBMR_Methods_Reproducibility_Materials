## Test: Initialize J=3 nu from J=2 solutions
## For each bootstrap rep, first solve 110 and 101, then use their solutions
## as starting points for 111
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Mental_Health_Data_Application")
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)
library(Matrix); library(numDeriv); library(nloptr)
source("MHD_functions.R")

original_data <- read.csv("data_application.csv")
percent <- original_data$Percentage
n <- 2486
class_n <- round(n * percent / 100)
dat <- gen_data(original_data, class_n, n)
dat$y <- dat$teacher_report

W_fn <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_alpha <- c("health", "father", "parent_report", "fp", "fh", "hp", "fhp")
h_nu <- function(data) cbind(health=data$health, father=data$father, parent_report=data$parent_report,
                              fp=data$fp, fh=data$fh, hp=data$hp, fhp=data$fhp)

ps_spec <- list(
  formula.list = list(
    r ~ teacher_report + health + father + teacher_report:health + teacher_report:father,
    r ~ teacher_report + health + parent_report + teacher_report:health + teacher_report:parent_report,
    r ~ teacher_report + parent_report + father + teacher_report:parent_report + teacher_report:father
  ),
  h_alpha.list = list(h_alpha, h_alpha, h_alpha),
  outcome = "teacher_report", inv_link = function(eta) 1/(1+exp(eta)), optimizer = "L-BFGS-B"
)

# Analytical SE from package
ebmr_pt <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat, W_fn)
res_pt <- ebmr_pt$EBMR_IPW(h_nu = h_nu, true_ps = NULL)
anal_se <- res_pt$se_ipw
cat(sprintf("Analytical SE: %.4f\n\n", anal_se))

B <- 300
cat(sprintf("=== J=2 init for J=3 nu (B=%d) ===\n\n", B))

mus_j2init <- mus_default <- rep(NA, B)

for (b in 1:B) {
  if (b %% 100 == 0) cat(sprintf("  rep %d...\n", b))
  set.seed(12345 + b)
  idx <- sample(1:n, n, replace = TRUE)
  dat_b <- dat[idx, ]

  tryCatch({
    ebmr_b <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat_b, W_fn)

    # Default 111
    res_d <- ebmr_b$EBMR_IPW(h_nu = h_nu, true_ps = NULL, se.fit = FALSE)
    mus_default[b] <- res_d$mu_ipw

    # Get J=2 solutions
    r110 <- ebmr_b$EBMR_IPW(h_nu = h_nu, model_indices = c(1,2), true_ps = NULL, se.fit = FALSE)
    r101 <- ebmr_b$EBMR_IPW(h_nu = h_nu, model_indices = c(1,3), true_ps = NULL, se.fit = FALSE)
    r011 <- ebmr_b$EBMR_IPW(h_nu = h_nu, model_indices = c(2,3), true_ps = NULL, se.fit = FALSE)

    # Build J=3 inits from J=2 solutions
    inits <- list(
      c(r110$nu.hat[1], r110$nu.hat[2], 0.001),  # from 110
      c(r101$nu.hat[1], 0.001, r101$nu.hat[2]),   # from 101
      c(0.001, r011$nu.hat[1], r011$nu.hat[2])    # from 011
    )

    # Try each init, compute objective, pick best
    ps_mat <- do.call(cbind, lapply(ebmr_b$ps_fit.list, function(pf) pf$fitted.values))
    h_x <- cbind(1, h_nu(dat_b))
    n_b <- nrow(dat_b)

    best_obj <- Inf; best_mu <- NA
    for (ini in inits) {
      tryCatch({
        res_s <- ebmr_b$EBMR_IPW(h_nu = h_nu, nu_init = ini, true_ps = NULL, se.fit = FALSE)
        ps_nu <- as.vector(ps_mat %*% res_s$nu.hat)
        g_mat <- as.vector(dat_b$r / ps_nu - 1) * h_x
        G <- colMeans(g_mat)
        obj <- sum(G^2)  # W=I objective for comparison
        if (obj < best_obj) { best_obj <- obj; best_mu <- res_s$mu_ipw }
      }, error = function(e) {})
    }
    mus_j2init[b] <- best_mu
  }, error = function(e) {})
}

cat(sprintf("\nDefault:  boot_se=%.4f, ratio=%.3f\n",
    sd(mus_default, na.rm=TRUE), anal_se / sd(mus_default, na.rm=TRUE)))
cat(sprintf("J2 init:  boot_se=%.4f, ratio=%.3f\n",
    sd(mus_j2init, na.rm=TRUE), anal_se / sd(mus_j2init, na.rm=TRUE)))
cat(sprintf("Agree: %d/%d (%.1f%%)\n",
    sum(abs(mus_default - mus_j2init) < 0.001, na.rm=TRUE),
    sum(!is.na(mus_default) & !is.na(mus_j2init)),
    100*mean(abs(mus_default - mus_j2init) < 0.001, na.rm=TRUE)))
