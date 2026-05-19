## Check whether non-normality of mu_ipw drives CP shortfall in
## EBMR_IPW_setting2-miss30-scenario8-1_1_n2000_replicate1000_test65.
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({ source("Basic_setup.r"); source("Data_Generation.r") })

f <- "Simulation_Results/EBMR_IPW_setting2-miss30-scenario8-1_1_n2000_replicate1000_test65.RDS"
x <- readRDS(f)
set.seed(12345); mu_true <- mean(setting2.A2(1e7)$y)
mu <- x["mu_ipw", ]; se <- x["se_ipw", ]
valid <- !is.na(mu) & !is.na(se) & is.finite(mu) & is.finite(se)
mu <- mu[valid]; se <- se[valid]

cat(sprintf("n=%d  mu_true=%.4f  mean(mu)=%.4f  median(mu)=%.4f  sd(mu)=%.4f\n",
            length(mu), mu_true, mean(mu), median(mu), sd(mu)))

# Skewness & excess kurtosis (Fisher: normal -> skew=0, kurt=0)
skew <- function(x) { x <- x-mean(x); mean(x^3) / mean(x^2)^1.5 }
kurt <- function(x) { x <- x-mean(x); mean(x^4) / mean(x^2)^2 - 3 }
se_skew <- sqrt(6/length(mu))   # asymptotic SE under normality
se_kurt <- sqrt(24/length(mu))
cat(sprintf("Skewness=%+.3f (SE %.3f under H0)  Excess kurtosis=%+.3f (SE %.3f)\n",
            skew(mu), se_skew, kurt(mu), se_kurt))
cat(sprintf("  z(skew)=%.2f  z(kurt)=%.2f  (|z|>1.96 -> reject normality at 5%%)\n",
            skew(mu)/se_skew, kurt(mu)/se_kurt))

# Shapiro-Wilk + Jarque-Bera
sw <- shapiro.test(mu)
jb_stat <- length(mu) * (skew(mu)^2/6 + kurt(mu)^2/24)
jb_p <- pchisq(jb_stat, df=2, lower.tail=FALSE)
cat(sprintf("\nShapiro-Wilk:  W=%.4f  p=%.4g\n", sw$statistic, sw$p.value))
cat(sprintf("Jarque-Bera:   JB=%.2f  p=%.4g\n", jb_stat, jb_p))

# Quantile comparison: empirical quantiles of (mu - mu_true)/se vs standard normal
zsc <- (mu - mu_true) / se
qs <- c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)
cat(sprintf("\n%-12s | %s\n", "quantile",
            "empirical (mu-mu_true)/se vs N(0,1)"))
for (q in qs) {
  cat(sprintf("  q=%5.3f    | empirical=%+.3f   N(0,1)=%+.3f\n",
              q, quantile(zsc, q), qnorm(q)))
}

# Symmetric CI coverage with t-quantile at different df
cat("\nCoverage with different CI critical values:\n")
for (crit in c(1.96, 2.00, qt(0.975, 1999), qt(0.975, 999))) {
  ci_lo <- mu - crit*se; ci_hi <- mu + crit*se
  cp <- mean(ci_lo <= mu_true & mu_true <= ci_hi)
  cat(sprintf("  crit=%.4f -> CP=%.3f\n", crit, cp))
}

# What if we use empirical quantiles of (mu-mu_true)/se to build CI?
emp_lo <- quantile(zsc, 0.025); emp_hi <- quantile(zsc, 0.975)
ci_lo <- mu - emp_hi*se; ci_hi <- mu - emp_lo*se
cp_emp <- mean(ci_lo <= mu_true & mu_true <= ci_hi)
cat(sprintf("\nCI from empirical z-quantiles [%.2f, %.2f] -> CP=%.3f (oracle-tuned, upper bound)\n",
            emp_lo, emp_hi, cp_emp))

# Save histogram
png("Simulation_Test/mu_ipw_dist_8-1_1_miss30.png", width=800, height=400)
par(mfrow=c(1,2))
hist(mu, breaks=40, main="mu_ipw distribution", xlab="mu_ipw")
abline(v=mu_true, col="red", lwd=2)
qqnorm(mu, main="Normal Q-Q plot of mu_ipw"); qqline(mu, col="blue")
dev.off()
cat("\nSaved histogram + Q-Q plot to Simulation_Test/mu_ipw_dist_8-1_1_miss30.png\n")
