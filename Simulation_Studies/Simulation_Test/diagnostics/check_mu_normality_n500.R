## Check normality of mu_ipw at n=500 for both miss50 and miss30.
## Hypothesis: non-normality of mu_ipw (heavier tails / skew) is the CP killer at small n.
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({ source("Basic_setup.r"); source("Data_Generation.r") })
set.seed(12345); mu_true <- mean(setting2.A1(1e7)$y)
cat(sprintf("mu_true = %.4f\n\n", mu_true))

skew <- function(x) { x <- x-mean(x); mean(x^3) / mean(x^2)^1.5 }
kurt <- function(x) { x <- x-mean(x); mean(x^4) / mean(x^2)^2 - 3 }

check_file <- function(label, f) {
  cat(strrep("=", 70), "\n", sep="")
  cat(sprintf("%s\n", label))
  cat(strrep("=", 70), "\n", sep="")
  x <- readRDS(f)
  mu <- x["mu_ipw", ]; se <- x["se_ipw", ]
  valid <- !is.na(mu) & !is.na(se) & is.finite(mu) & is.finite(se)
  mu <- mu[valid]; se <- se[valid]

  cat(sprintf("n=%d  mean(mu)=%.4f  median(mu)=%.4f  sd(mu)=%.4f\n",
              length(mu), mean(mu), median(mu), sd(mu)))

  s <- skew(mu); k <- kurt(mu); n <- length(mu)
  se_s <- sqrt(6/n); se_k <- sqrt(24/n)
  cat(sprintf("Skewness=%+.3f (z=%.2f)  Excess kurtosis=%+.3f (z=%.2f)\n",
              s, s/se_s, k, k/se_k))
  cat(sprintf("  (|z|>1.96 -> reject normality at 5%%)\n"))

  sw <- shapiro.test(mu)
  jb_stat <- n * (s^2/6 + k^2/24)
  jb_p <- pchisq(jb_stat, df=2, lower.tail=FALSE)
  cat(sprintf("Shapiro-Wilk: W=%.4f p=%.4g   Jarque-Bera: JB=%.2f p=%.4g\n",
              sw$statistic, sw$p.value, jb_stat, jb_p))

  cat("\nEmpirical quantiles of (mu - mu_true)/se vs N(0,1):\n")
  zsc <- (mu - mu_true) / se
  qs <- c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)
  for (q in qs) {
    cat(sprintf("  q=%.3f  empirical=%+7.3f   N(0,1)=%+7.3f   t(2)=%+7.3f   t(5)=%+7.3f\n",
                q, quantile(zsc, q), qnorm(q), qt(q, 2), qt(q, 5)))
  }

  cat("\nCoverage with different critical values:\n")
  for (info in list(list(c=1.96, n="N(0,1) z=1.96 (95%)"),
                    list(c=qt(0.975, 5),  n="t(df=5) (heavy tail)"),
                    list(c=qt(0.975, 2),  n="t(df=2)"),
                    list(c=qt(0.95, 1.96^2/qnorm(0.975)^2 * Inf), n="adaptive (skip)"))) {
    if (is.null(info$c) || !is.finite(info$c)) next
    ci_lo <- mu - info$c*se; ci_hi <- mu + info$c*se
    cp <- mean(ci_lo <= mu_true & mu_true <= ci_hi)
    cat(sprintf("  crit=%.4f  CP=%.3f  (%s)\n", info$c, cp, info$n))
  }
  # Empirical-quantile CI (oracle upper bound)
  emp_lo <- quantile(zsc, 0.025); emp_hi <- quantile(zsc, 0.975)
  ci_lo <- mu - emp_hi*se; ci_hi <- mu - emp_lo*se
  cp_emp <- mean(ci_lo <= mu_true & mu_true <= ci_hi)
  cat(sprintf("  CI from empirical z-quantiles [%.2f, %.2f] -> CP=%.3f (oracle, upper bound)\n",
              emp_lo, emp_hi, cp_emp))

  # Save Q-Q plot
  png_path <- sprintf("Simulation_Test/qq_%s.png", gsub(" ", "_", label))
  png(png_path, width=800, height=400)
  par(mfrow=c(1,2))
  hist(mu, breaks=40, main=paste("mu_ipw:", label), xlab="mu_ipw")
  abline(v=mu_true, col="red", lwd=2)
  qqnorm(mu, main="Normal Q-Q"); qqline(mu, col="blue")
  dev.off()
  cat(sprintf("\nSaved: %s\n\n", png_path))
}

check_file("n500_miss50",
           "Simulation_Results/EBMR_IPW_setting2-miss50-scenario8-1_1_n500_replicate1000_test65.RDS")
check_file("n500_miss30",
           "Simulation_Results/EBMR_IPW_setting2-miss30-scenario8-1_1_n500_replicate1000_test65.RDS")
check_file("n2000_miss30 (reference)",
           "Simulation_Results/EBMR_IPW_setting2-miss30-scenario8-1_1_n2000_replicate1000_test65.RDS")
