## Fix: Initialize J=3 nu from J=2 solutions instead of uniform (1/3,1/3,1/3)
## Test on 300 bootstrap reps
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
    r ~ teacher_report + health + father + teacher_report:health + teacher_report:father,
    r ~ teacher_report + health + parent_report + teacher_report:health + teacher_report:parent_report,
    r ~ teacher_report + parent_report + father + teacher_report:parent_report + teacher_report:father
  ),
  h_alpha.list = list(h_alpha, h_alpha, h_alpha),
  outcome = "teacher_report", inv_link = function(eta) 1/(1+exp(eta)), optimizer = "L-BFGS-B"
)

# Point estimate
ebmr <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat, W)
res_default <- ebmr$EBMR_IPW(h_nu = h_nu, true_ps = NULL)
cat(sprintf("Default 111: mu=%.4f, se=%.4f, w=(%s)\n",
    res_default$mu_ipw, res_default$se_ipw, paste(round(res_default$w.hat, 4), collapse=",")))

# Multi-start point estimate
starts <- list(
  c(1/3, 1/3, 1/3),      # uniform
  c(1, 0.001, 0.001),     # M1
  c(0.001, 1, 0.001),     # M2
  c(0.001, 0.001, 1)      # M3
)

# Also get J=2 solutions to use as starts
r110 <- ebmr$EBMR_IPW(h_nu = h_nu, model_indices = c(1,2), true_ps = NULL, se.fit = FALSE)
r101 <- ebmr$EBMR_IPW(h_nu = h_nu, model_indices = c(1,3), true_ps = NULL, se.fit = FALSE)
starts[[5]] <- c(r110$nu.hat[1], r110$nu.hat[2], 0.001)  # from 110
starts[[6]] <- c(r101$nu.hat[1], 0.001, r101$nu.hat[2])  # from 101

compute_nu_obj <- function(ebmr_obj, nu_hat) {
  ps_mat <- do.call(cbind, lapply(ebmr_obj$ps_fit.list, function(pf) pf$fitted.values))
  h_x <- cbind(1, h_nu(ebmr_obj$data))
  nn <- nrow(ebmr_obj$data)
  r_vec <- ebmr_obj$data$r
  ps_nu <- as.vector(ps_mat %*% nu_hat)
  g_mat <- as.vector(r_vec / ps_nu - 1) * h_x
  G_vec <- colMeans(g_mat)
  W_hat <- tryCatch(solve(crossprod(g_mat) / nn), error = function(e) diag(ncol(h_x)))
  as.numeric(t(G_vec) %*% W_hat %*% G_vec)
}

cat("\nPoint estimate with multi-start:\n")
best_obj <- Inf; best_res <- NULL
for (si in seq_along(starts)) {
  tryCatch({
    res_s <- ebmr$EBMR_IPW(h_nu = h_nu, nu_init = starts[[si]], true_ps = NULL, se.fit = FALSE)
    obj <- compute_nu_obj(ebmr, res_s$nu.hat)
    cat(sprintf("  Start %d: w=(%s), mu=%.4f, obj=%.4e\n", si,
        paste(round(res_s$w.hat, 4), collapse=","), res_s$mu_ipw, obj))
    if (obj < best_obj) { best_obj <- obj; best_res <- res_s }
  }, error = function(e) {})
}
cat(sprintf("  BEST: mu=%.4f, obj=%.4e\n\n", best_res$mu_ipw, best_obj))

# Bootstrap: compare default vs multi-start
B <- 300
cat(sprintf("=== Bootstrap (B=%d): default vs multi-start ===\n\n", B))

mus_default <- mus_multi <- rep(NA, B)

for (b in 1:B) {
  if (b %% 100 == 0) cat(sprintf("  rep %d...\n", b))
  set.seed(12345 + b)
  idx <- sample(1:n, n, replace = TRUE)
  dat_b <- dat[idx, ]

  tryCatch({
    ebmr_b <- EBMRAlgorithmFast4$new("teacher_report", ps_spec, dat_b, W)

    # Default (uniform init)
    res_d <- ebmr_b$EBMR_IPW(h_nu = h_nu, true_ps = NULL, se.fit = FALSE)
    mus_default[b] <- res_d$mu_ipw

    # Multi-start: uniform + M1 + M2 + M3 + from J=2 solutions
    r110_b <- ebmr_b$EBMR_IPW(h_nu = h_nu, model_indices = c(1,2), true_ps = NULL, se.fit = FALSE)
    r101_b <- ebmr_b$EBMR_IPW(h_nu = h_nu, model_indices = c(1,3), true_ps = NULL, se.fit = FALSE)

    starts_b <- list(
      c(1/3, 1/3, 1/3),
      c(1, 0.001, 0.001),
      c(0.001, 1, 0.001),
      c(0.001, 0.001, 1),
      c(r110_b$nu.hat[1], r110_b$nu.hat[2], 0.001),
      c(r101_b$nu.hat[1], 0.001, r101_b$nu.hat[2])
    )

    best_obj_b <- Inf; best_mu_b <- NA
    for (si in seq_along(starts_b)) {
      tryCatch({
        res_s <- ebmr_b$EBMR_IPW(h_nu = h_nu, nu_init = starts_b[[si]], true_ps = NULL, se.fit = FALSE)
        obj <- compute_nu_obj(ebmr_b, res_s$nu.hat)
        if (obj < best_obj_b) { best_obj_b <- obj; best_mu_b <- res_s$mu_ipw }
      }, error = function(e) {})
    }
    mus_multi[b] <- best_mu_b
  }, error = function(e) {})
}

cat(sprintf("\nDefault:     boot_se=%.4f, ratio=%.3f\n",
    sd(mus_default, na.rm=TRUE), res_default$se_ipw / sd(mus_default, na.rm=TRUE)))
cat(sprintf("Multi-start: boot_se=%.4f, ratio=%.3f\n",
    sd(mus_multi, na.rm=TRUE), res_default$se_ipw / sd(mus_multi, na.rm=TRUE)))
cat(sprintf("\nAgree: %d/%d (%.1f%%)\n",
    sum(abs(mus_default - mus_multi) < 0.001, na.rm=TRUE),
    sum(!is.na(mus_default) & !is.na(mus_multi)),
    100*mean(abs(mus_default - mus_multi) < 0.001, na.rm=TRUE)))
