setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

f <- "Simulation_Results/EBMR_IPW_setting3-miss50-scenario9-3_2_n2000_replicate1000_test59.RDS"
res <- readRDS(f)
mu <- res[1, ]; se <- res[3, ]
valid <- !is.na(mu) & !is.na(se)
mu_v <- mu[valid]; se_v <- se[valid]

png("plot_mu_se_distribution.png", width = 1200, height = 800)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# 1. Histogram of mu_ipw
hist(mu_v, breaks = 80, main = "mu_ipw distribution", xlab = "mu_ipw", col = "steelblue")
abline(v = median(mu_v), col = "red", lwd = 2)
# IQR bounds
q1 <- quantile(mu_v, 0.25); q3 <- quantile(mu_v, 0.75); iqr <- q3 - q1
abline(v = q1 - 3*iqr, col = "orange", lwd = 2, lty = 2)
abline(v = q3 + 3*iqr, col = "orange", lwd = 2, lty = 2)
legend("topright", c("median", "3*IQR bounds"), col = c("red", "orange"), lwd = 2, lty = c(1, 2))

# 2. Histogram of se_ipw
hist(se_v, breaks = 80, main = "se_ipw distribution", xlab = "se_ipw", col = "steelblue")
abline(v = median(se_v), col = "red", lwd = 2)

# 3. Scatter: mu_ipw vs se_ipw
plot(mu_v, se_v, pch = 16, cex = 0.6, col = rgb(0, 0, 0.7, 0.5),
     main = "mu_ipw vs se_ipw", xlab = "mu_ipw", ylab = "se_ipw")
# Mark 3*IQR mu outliers
is_mu_out <- mu_v < (q1 - 3*iqr) | mu_v > (q3 + 3*iqr)
points(mu_v[is_mu_out], se_v[is_mu_out], pch = 16, col = "red", cex = 1)
legend("topright", c("normal", "mu outlier (3*IQR)"),
       col = c(rgb(0, 0, 0.7, 0.5), "red"), pch = 16)

# 4. Zoomed scatter (exclude extreme tails for readability)
xlim <- quantile(mu_v, c(0.01, 0.99))
ylim <- quantile(se_v, c(0, 0.99))
plot(mu_v, se_v, pch = 16, cex = 0.6, col = rgb(0, 0, 0.7, 0.5),
     main = "mu_ipw vs se_ipw (zoomed)", xlab = "mu_ipw", ylab = "se_ipw",
     xlim = xlim, ylim = ylim)
points(mu_v[is_mu_out], se_v[is_mu_out], pch = 16, col = "red", cex = 1)

dev.off()
cat("Saved: plot_mu_se_distribution.png\n")
