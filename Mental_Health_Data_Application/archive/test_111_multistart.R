## Test: multistart nu for 111 ensemble across bootstrap reps
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

ps_specifications <- list(
  formula.list = list(
    r ~ teacher_report + health + father + teacher_report:health + teacher_report:father + father:health,
    r ~ teacher_report + parent_report + health + teacher_report:parent_report + teacher_report:health + parent_report:health,
    r ~ teacher_report + parent_report + father + teacher_report:parent_report + teacher_report:father + parent_report:father
  ),
  h_alpha.list = list(
    c("health", "father", "parent_report", "fp", "fh", "hp", "fhp"),
    c("health", "father", "parent_report", "fp", "fh", "hp", "fhp"),
    c("health", "father", "parent_report", "fp", "fh", "hp", "fhp")
  ),
  outcome = "teacher_report",
  inv_link = function(eta) 1 / (1 + exp(eta)),
  optimizer = "L-BFGS-B"
)

h_nu <- function(data) cbind(health=data$health, father=data$father, parent_report=data$parent_report,
                              fp=data$fp, fh=data$fh, hp=data$hp, fhp=data$fhp)

J <- 3
nu_starts <- list(rep(1/J, J))
for (k in 1:J) { s <- rep(0.001, J); s[k] <- 1; nu_starts[[k+1]] <- s }

B <- 500
cat("=== 111 multistart nu (B=500) ===\n\n")

# Single start (default) vs multistart
mus_single <- mus_multi <- rep(NA, B)
w_single <- w_multi <- matrix(NA, B, 3)

for (b in 1:B) {
  if (b %% 100 == 0) cat(sprintf("  rep %d...\n", b))
  set.seed(12345 + b)
  idx <- sample(1:n, n, replace = TRUE)
  dat_b <- dat[idx, ]

  tryCatch({
    ebmr_b <- EBMRAlgorithmFast4$new("teacher_report", ps_specifications, dat_b, W)

    # Single start (default)
    res_s <- ebmr_b$EBMR_IPW(h_nu = h_nu, true_ps = NULL, se.fit = FALSE)
    mus_single[b] <- res_s$mu_ipw
    w_single[b, ] <- res_s$w.hat

    # Multistart: try all starts, pick lowest objective
    ps_mat <- do.call(cbind, lapply(ebmr_b$ps_fit.list, function(pf) pf$fitted.values))
    h_x2 <- h_nu(dat_b)
    h_x <- cbind(1, h_x2)
    n_b <- nrow(dat_b)
    r_b <- dat_b$r

    best_obj <- Inf
    best_res <- NULL

    for (si in seq_along(nu_starts)) {
      tryCatch({
        res_m <- ebmr_b$EBMR_IPW(h_nu = h_nu, nu_init = nu_starts[[si]],
                                   true_ps = NULL, se.fit = FALSE)
        # Compute objective
        ps_nu <- as.vector(ps_mat %*% res_m$nu.hat)
        g_mat <- as.vector(r_b / ps_nu - 1) * h_x
        G_vec <- colMeans(g_mat)
        W_hat <- tryCatch(solve(crossprod(g_mat) / n_b), error = function(e) diag(ncol(h_x)))
        obj <- as.numeric(t(G_vec) %*% W_hat %*% G_vec)

        if (obj < best_obj) {
          best_obj <- obj
          best_res <- res_m
        }
      }, error = function(e) {})
    }

    if (!is.null(best_res)) {
      mus_multi[b] <- best_res$mu_ipw
      w_multi[b, ] <- best_res$w.hat
    }
  }, error = function(e) {})
}

cat("\n=== Results ===\n\n")

valid_s <- !is.na(mus_single)
valid_m <- !is.na(mus_multi)

cat(sprintf("Single start: valid=%d, mu_sd=%.4f\n", sum(valid_s), sd(mus_single, na.rm=TRUE)))
cat(sprintf("Multi  start: valid=%d, mu_sd=%.4f\n\n", sum(valid_m), sd(mus_multi, na.rm=TRUE)))

# Weight comparison
cat("Single start weights:\n")
for (j in 1:3) cat(sprintf("  w%d: mean=%.4f, sd=%.4f, >0.9: %d\n",
    j, mean(w_single[valid_s,j]), sd(w_single[valid_s,j]), sum(w_single[valid_s,j] > 0.9)))

cat("\nMulti start weights:\n")
for (j in 1:3) cat(sprintf("  w%d: mean=%.4f, sd=%.4f, >0.9: %d\n",
    j, mean(w_multi[valid_m,j]), sd(w_multi[valid_m,j]), sum(w_multi[valid_m,j] > 0.9)))

# Check: how often do single and multi agree?
agree <- valid_s & valid_m & abs(mus_single - mus_multi) < 0.001
cat(sprintf("\nAgree (|mu diff| < 0.001): %d/%d (%.1f%%)\n",
    sum(agree), sum(valid_s & valid_m), 100*mean(agree[valid_s & valid_m])))

disagree <- valid_s & valid_m & abs(mus_single - mus_multi) >= 0.001
if (sum(disagree) > 0) {
  cat(sprintf("Disagree: %d reps\n", sum(disagree)))
  cat(sprintf("  Single mu_sd in disagree reps: %.4f\n", sd(mus_single[disagree])))
  cat(sprintf("  Multi  mu_sd in disagree reps: %.4f\n", sd(mus_multi[disagree])))
}
