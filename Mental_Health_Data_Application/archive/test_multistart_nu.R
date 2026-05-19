## Test: Multi-start nu estimation for health=1 analysis
## Try J+1 starting values: (1/J,...,1/J) and (1,0,...,0), (0,1,...,0), (0,0,...,1)
## Pick the one with lowest GMM objective
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
    r ~ teacher_report + father,
    r ~ teacher_report + parent_report,
    r ~ father + parent_report
  ),
  h_alpha.list = list(
    c("father", "parent_report"),
    c("father", "parent_report"),
    c("father", "parent_report")
  ),
  outcome = "teacher_report",
  inv_link = function(eta) 1 / (1 + exp(eta))
)

h_nu <- function(data) {
  cbind(father = data$father, parent_report = data$parent_report,
        fp = data$father * data$parent_report)
}

all_model_sets <- list(c(1), c(2), c(3), c(1,2), c(1,3), c(2,3), c(1,2,3))
model_set_labels <- c("100", "010", "001", "110", "101", "011", "111")

# Fit PS models once
ebmr <- EBMRAlgorithmFast4$new("teacher_report", ps_specifications, dat, W)

cat("=== Multi-start nu estimation for health=1 ===\n\n")

# For each multi-model combination, try multiple nu starts
for (j in 4:7) {  # only ensembles (110, 101, 011, 111)
  model_set <- all_model_sets[[j]]
  label <- model_set_labels[j]
  J <- length(model_set)

  # Generate J+1 starting values
  starts <- list(rep(1/J, J))  # uniform
  for (k in 1:J) {
    s <- rep(0.001, J)
    s[k] <- 1
    starts[[k + 1]] <- s
  }

  cat(sprintf("--- %s (J=%d, %d starts) ---\n", label, J, length(starts)))

  best_obj <- Inf
  best_result <- NULL

  for (s_idx in seq_along(starts)) {
    nu_init <- starts[[s_idx]]
    tryCatch({
      res <- ebmr$EBMR_IPW(h_nu = h_nu, model_indices = model_set,
                             nu_init = nu_init, true_ps = NULL)

      # Compute GMM objective for nu
      ps_mat <- do.call(cbind, lapply(ebmr$ps_fit.list[model_set], function(pf) pf$fitted.values))
      h_x2 <- h_nu(dat)
      h_x <- cbind(1, h_x2)
      ps_nu <- as.vector(ps_mat %*% res$nu.hat)
      g_mat <- as.vector(dat$r / ps_nu - 1) * h_x
      G_vec <- colMeans(g_mat)
      W_hat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(ncol(h_x)))
      obj <- as.numeric(t(G_vec) %*% W_hat %*% G_vec)

      cat(sprintf("  Start %d: nu=(%s) -> nu_hat=(%s), w=(%s), obj=%.4e, mu=%.4f, se=%.4f\n",
          s_idx, paste(round(nu_init, 3), collapse=","),
          paste(round(res$nu.hat, 4), collapse=","),
          paste(round(res$w.hat, 4), collapse=","),
          obj, res$mu_ipw, res$se_ipw))

      if (obj < best_obj) {
        best_obj <- obj
        best_result <- res
      }
    }, error = function(e) {
      cat(sprintf("  Start %d: ERROR: %s\n", s_idx, e$message))
    })
  }

  cat(sprintf("  BEST: mu=%.4f, se=%.4f, w=(%s), obj=%.4e\n\n",
      best_result$mu_ipw, best_result$se_ipw,
      paste(round(best_result$w.hat, 4), collapse=","), best_obj))
}

# Now run perturbation bootstrap with multi-start
cat("\n=== Perturbation bootstrap with multi-start nu ===\n\n")

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

  r_idx <- dat$r == 1
  mu_cc_b <- sum(wt[r_idx] * dat$teacher_report[r_idx]) / sum(wt[r_idx])

  mu_vec <- setNames(rep(NA, 7), model_set_labels)

  for (j in seq_along(all_model_sets)) {
    model_set <- all_model_sets[[j]]
    label <- model_set_labels[j]
    J_local <- length(model_set)

    ps_spec_b <- list(
      formula.list = ps_specifications$formula.list[model_set],
      h_alpha.list = ps_specifications$h_alpha.list[model_set],
      outcome = ps_specifications$outcome,
      inv_link = ps_specifications$inv_link
    )

    tryCatch({
      ebmr_b <- EBMRAlgorithmFast4$new("teacher_report", ps_spec_b, dat, W, wt = wt)

      if (J_local == 1) {
        # Single model, no multi-start needed
        res_b <- ebmr_b$EBMR_IPW(h_nu = h_nu, true_ps = NULL, se.fit = FALSE, wt = wt)
        mu_vec[j] <- res_b$mu_ipw
      } else {
        # Multi-start: try J+1 starting values, pick lowest objective
        starts <- list(rep(1/J_local, J_local))
        for (k in 1:J_local) {
          s <- rep(0.001, J_local)
          s[k] <- 1
          starts[[k + 1]] <- s
        }

        best_obj <- Inf
        best_mu <- NA

        for (s_idx in seq_along(starts)) {
          tryCatch({
            res_s <- ebmr_b$EBMR_IPW(h_nu = h_nu, nu_init = starts[[s_idx]],
                                       true_ps = NULL, se.fit = FALSE, wt = wt)
            # Compute objective
            ps_mat <- do.call(cbind, lapply(ebmr_b$ps_fit.list, function(pf) pf$fitted.values))
            h_x2 <- h_nu(dat)
            h_x <- cbind(1, h_x2)
            ps_nu <- as.vector(ps_mat %*% res_s$nu.hat)
            g_mat <- as.vector(wt * (dat$r / ps_nu - 1)) * h_x
            G_vec <- colMeans(g_mat)
            W_hat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(ncol(h_x)))
            obj <- as.numeric(t(G_vec) %*% W_hat %*% G_vec)

            if (obj < best_obj) {
              best_obj <- obj
              best_mu <- res_s$mu_ipw
            }
          }, error = function(e) {})
        }
        mu_vec[j] <- best_mu
      }
    }, error = function(e) {})
  }

  c(CC = mu_cc_b, mu_vec)
}

close(pb)
stopCluster(cl)

# Point estimates with multi-start
ebmr_pt <- EBMRAlgorithmFast4$new("teacher_report", ps_specifications, dat, W)
point_mu <- point_se <- numeric(7)
for (j in 1:7) {
  model_set <- all_model_sets[[j]]
  J_local <- length(model_set)
  if (J_local == 1) {
    res <- ebmr_pt$EBMR_IPW(h_nu = h_nu, model_indices = model_set, true_ps = NULL)
    point_mu[j] <- res$mu_ipw; point_se[j] <- res$se_ipw
  } else {
    starts <- list(rep(1/J_local, J_local))
    for (k in 1:J_local) { s <- rep(0.001, J_local); s[k] <- 1; starts[[k+1]] <- s }
    best_obj <- Inf; best_res <- NULL
    for (s_idx in seq_along(starts)) {
      tryCatch({
        res <- ebmr_pt$EBMR_IPW(h_nu = h_nu, model_indices = model_set,
                                  nu_init = starts[[s_idx]], true_ps = NULL)
        ps_mat <- do.call(cbind, lapply(ebmr_pt$ps_fit.list[model_set], function(pf) pf$fitted.values))
        h_x <- cbind(1, h_nu(dat))
        ps_nu <- as.vector(ps_mat %*% res$nu.hat)
        g_mat <- as.vector(dat$r / ps_nu - 1) * h_x
        G_vec <- colMeans(g_mat)
        W_hat <- tryCatch(solve(crossprod(g_mat) / n), error = function(e) diag(ncol(h_x)))
        obj <- as.numeric(t(G_vec) %*% W_hat %*% G_vec)
        if (obj < best_obj) { best_obj <- obj; best_res <- res }
      }, error = function(e) {})
    }
    point_mu[j] <- best_res$mu_ipw; point_se[j] <- best_res$se_ipw
  }
}

mu_cc <- mean(dat$teacher_report[dat$r == 1])
se_cc <- sd(dat$teacher_report[dat$r == 1]) / sqrt(sum(dat$r))

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

compute_se_nocap <- function(x) {
  x <- x[!is.na(x)]; n_tot <- length(x)
  if (n_tot < 10) return(list(se = NA, n_out = NA))
  Q1 <- quantile(x, 0.25); Q3 <- quantile(x, 0.75); IQR_val <- Q3 - Q1
  is_out <- x < (Q1 - 3 * IQR_val) | x > (Q3 + 3 * IQR_val)
  list(se = sd(x[!is_out]), n_out = sum(is_out))
}

cat("\n\n=== Results (multi-start nu) ===\n\n")
labels <- c("CC", model_set_labels)
estimates <- c(mu_cc, point_mu)
analytical <- c(se_cc, point_se)

cat(sprintf("  %-8s %10s %12s %12s %8s %12s %8s\n",
            "Label", "Estimate", "Analytical", "Ptb(1%%cap)", "Out", "Ptb(nocap)", "Out"))
cat("  ", paste(rep("-", 78), collapse = ""), "\n")

for (jj in 1:8) {
  bt1 <- compute_se_capped(perturb_mat[, jj])
  bt2 <- compute_se_nocap(perturb_mat[, jj])
  n_na <- sum(is.na(perturb_mat[, jj]))
  cat(sprintf("  %-8s %10.4f %12.4f %12.4f %8d %12.4f %8d  (NA:%d)\n",
              labels[jj], estimates[jj], analytical[jj],
              bt1$se, bt1$n_out, bt2$se, bt2$n_out, n_na))
}
cat("  ", paste(rep("-", 78), collapse = ""), "\n")
