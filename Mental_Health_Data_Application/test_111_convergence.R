## Test: Is 111 ensemble nu instability caused by loose convergence?
## h_nu = h_alpha = (health, father, parent_report, fp, fh, hp, fhp)
## Track nu objective and gradient at convergence across bootstrap reps
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

# h_nu = h_alpha (same)
h_nu <- function(data) cbind(health=data$health, father=data$father, parent_report=data$parent_report,
                              fp=data$fp, fh=data$fh, hp=data$hp, fhp=data$fhp)

B <- 300
cat("=== 111 nu convergence diagnostics (B=300) ===\n\n")

mus <- rep(NA, B)
w_mat <- matrix(NA, B, 3)
nu_obj <- nu_grad <- rep(NA, B)

for (b in 1:B) {
  if (b %% 100 == 0) cat(sprintf("  rep %d...\n", b))
  set.seed(12345 + b)
  idx <- sample(1:n, n, replace = TRUE)
  dat_b <- dat[idx, ]

  tryCatch({
    ebmr_b <- EBMRAlgorithmFast4$new("teacher_report", ps_specifications, dat_b, W)
    res_b <- ebmr_b$EBMR_IPW(h_nu = h_nu, true_ps = NULL, se.fit = FALSE)
    mus[b] <- res_b$mu_ipw
    w_mat[b, ] <- res_b$w.hat

    # Get nu GMM convergence info from the ensemble fit
    # Recompute nu objective and gradient manually
    ps_mat <- do.call(cbind, lapply(ebmr_b$ps_fit.list, function(pf) pf$fitted.values))
    h_x <- cbind(1, h_nu(dat_b))  # ensemble adds intercept
    n_b <- nrow(dat_b)
    r_b <- dat_b$r
    nu_hat <- res_b$nu.hat

    ps_nu <- as.vector(ps_mat %*% nu_hat)
    g_mat <- as.vector(r_b / ps_nu - 1) * h_x
    G_vec <- colMeans(g_mat)
    W_hat <- tryCatch(solve(crossprod(g_mat) / n_b), error = function(e) diag(ncol(h_x)))
    nu_obj[b] <- as.numeric(t(G_vec) %*% W_hat %*% G_vec)

    # Gradient of nu objective
    neg_r_ps2 <- -r_b / (ps_nu^2)
    Gamma_nu <- crossprod(h_x * neg_r_ps2, ps_mat) / n_b
    nu_grad[b] <- max(abs(2 * as.vector(crossprod(Gamma_nu, W_hat %*% G_vec))))
  }, error = function(e) {})
}

valid <- !is.na(mus)
cat(sprintf("\nValid: %d/%d\n", sum(valid), B))
cat(sprintf("mu: mean=%.4f, sd=%.4f\n\n", mean(mus, na.rm=TRUE), sd(mus, na.rm=TRUE)))

# Nu convergence quality
cat("=== Nu convergence quality ===\n")
cat(sprintf("  nu_obj:  mean=%.2e, med=%.2e, max=%.2e\n",
    mean(nu_obj, na.rm=TRUE), median(nu_obj, na.rm=TRUE), max(nu_obj, na.rm=TRUE)))
cat(sprintf("  nu_grad: mean=%.2e, med=%.2e, max=%.2e\n",
    mean(nu_grad, na.rm=TRUE), median(nu_grad, na.rm=TRUE), max(nu_grad, na.rm=TRUE)))
cat(sprintf("  nu_grad > 1e-4: %d\n", sum(nu_grad > 1e-4, na.rm=TRUE)))
cat(sprintf("  nu_grad > 1e-3: %d\n", sum(nu_grad > 1e-3, na.rm=TRUE)))
cat(sprintf("  nu_grad > 1e-2: %d\n\n", sum(nu_grad > 1e-2, na.rm=TRUE)))

# Does convergence quality correlate with weight switching?
m1_dom <- valid & w_mat[,1] > 0.9
non_m1 <- valid & w_mat[,1] <= 0.9
cat("=== M1-dominated vs non-M1 reps ===\n")
cat(sprintf("  M1-dominated (%d reps): nu_obj mean=%.2e, nu_grad mean=%.2e\n",
    sum(m1_dom), mean(nu_obj[m1_dom], na.rm=TRUE), mean(nu_grad[m1_dom], na.rm=TRUE)))
cat(sprintf("  Non-M1 (%d reps):       nu_obj mean=%.2e, nu_grad mean=%.2e\n",
    sum(non_m1), mean(nu_obj[non_m1], na.rm=TRUE), mean(nu_grad[non_m1], na.rm=TRUE)))

cat(sprintf("\n  cor(nu_obj, mu): %.3f\n", cor(nu_obj[valid], mus[valid], use="complete.obs")))
cat(sprintf("  cor(nu_grad, mu): %.3f\n", cor(nu_grad[valid], mus[valid], use="complete.obs")))
cat(sprintf("  cor(nu_grad, w1): %.3f\n", cor(nu_grad[valid], w_mat[valid,1], use="complete.obs")))
