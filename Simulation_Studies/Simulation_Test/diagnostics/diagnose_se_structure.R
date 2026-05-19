setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Data_Generation.r")
source("config/scenarios.R")
library(EBMRalgorithmFast4)
source("Basic_setup.r")

# Identify outlier reps with normal alpha
f <- "Simulation_Results/EBMR_IPW_setting3-miss50-scenario9-3_2_n2000_replicate1000_test59.RDS"
res <- readRDS(f)
mu <- res[1,]; se <- res[3,]
valid <- !is.na(mu) & !is.na(se)
mu_v <- mu[valid]; se_v <- se[valid]

q1 <- quantile(mu_v, 0.25); q3 <- quantile(mu_v, 0.75); iqr <- q3 - q1
q1s <- quantile(se_v, 0.25); q3s <- quantile(se_v, 0.75); iqrs <- q3s - q1s
is_out <- (mu_v < q1 - 3*iqr | mu_v > q3 + 3*iqr) | (se_v < q1s - 3*iqrs | se_v > q3s + 3*iqrs)

alpha_mat <- rbind(res["alpha.hat1", valid], res["alpha.hat2", valid],
                   res["alpha.hat3", valid], res["alpha.hat4", valid])
alpha_norm <- apply(alpha_mat, 2, function(a) sqrt(sum(a^2)))

normal_alpha_out <- which(valid)[which(is_out & alpha_norm <= 5)]
normal_alpha_out <- normal_alpha_out[order(se[normal_alpha_out], decreasing = TRUE)]
normal_reps <- c(1, 2, 3, 4, 5)

check_reps <- c(normal_alpha_out[1:5], normal_reps)

data_file <- misspecified_model_all_data_file.list$setting3$miss50[[1]]
all_data <- readRDS(data_file)
n_val <- 2000

ps_spec <- get_ps_spec("9-alt1")
subset_ps_spec <- list(
  formula.list = ps_spec$formula.list[2],
  h_alpha.list = ps_spec$h_alpha.list[2],
  inv_link = ps_spec$inv_link,
  outcome = ps_spec$outcome
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

for (rep_i in check_reps) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]
  label <- if (rep_i %in% normal_alpha_out) "OUT" else "NORM"

  tryCatch({
    ebmr <- EBMRAlgorithmFast4$new("y", subset_ps_spec, dat, W_func)

    # Use EBMR_IPW to get the actual SE, then also manually decompose
    res_ipw <- ebmr$EBMR_IPW(
      h_nu = function(dat) cbind(u1 = dat$u1, u2 = dat$u2, z1 = dat$z1, z2 = dat$z2),
      se.fit = TRUE, type = "HT"
    )

    ps_fit <- ebmr$ps_fit.list[[1]]
    ps_fitted <- ps_fit$fitted.values
    r <- as.numeric(dat$r)
    y <- dat$y
    n <- length(y)
    alpha_hat <- ps_fit$coefficients

    cat(sprintf("=== Rep %d (%s) ===\n", rep_i, label))
    cat(sprintf("  ||alpha||=%.2f, mu=%.4f, se=%.4f\n",
                sqrt(sum(alpha_hat^2)), res_ipw$mu_ipw, res_ipw$se_ipw))
    cat(sprintf("  PS: min=%.6f, max=%.4f, n(ps<0.01)=%d, n(ps<0.05)=%d\n",
                min(ps_fitted), max(ps_fitted), sum(ps_fitted < 0.01), sum(ps_fitted < 0.05)))
    cat(sprintf("  IPW: max(r/ps)=%.1f, sum(r/ps)=%.1f\n",
                max(r/ps_fitted), sum(r/ps_fitted)))

    # Base term variance: var(r/ps*y)/n — this is the "naive" SE^2 ignoring alpha estimation
    base <- r / ps_fitted * y
    se_naive <- sqrt(var(base)/n)
    cat(sprintf("  se_naive (no alpha adj) = %.4f\n", se_naive))
    cat(sprintf("  se_actual / se_naive = %.2f\n", res_ipw$se_ipw / se_naive))

    # Distribution of base term
    cat(sprintf("  base: sd=%.4f, skew=%.2f, kurt=%.2f\n",
                sd(base), mean((base-mean(base))^3)/sd(base)^3,
                mean((base-mean(base))^4)/sd(base)^4))

    # Top contributors
    top5 <- order(abs(base - mean(base)), decreasing=TRUE)[1:5]
    cat(sprintf("  top5 |base-mean|: %.1f (ps=%.5f,r=%d,y=%.2f), %.1f, %.1f, %.1f, %.1f\n",
                abs(base[top5[1]]-mean(base)), ps_fitted[top5[1]], r[top5[1]], y[top5[1]],
                abs(base[top5[2]]-mean(base)),
                abs(base[top5[3]]-mean(base)),
                abs(base[top5[4]]-mean(base)),
                abs(base[top5[5]]-mean(base))))
    cat("\n")

  }, error = function(e) {
    cat(sprintf("=== Rep %d (%s) === ERROR: %s\n\n", rep_i, label, conditionMessage(e)))
  })
}

cat("Done!\n")
