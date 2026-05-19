## Re-confirm s1.B1 n=2000, M3 only, constrained_nr cond=1e3, 1000 reps.
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({ source("Basic_setup.r"); source("Data_Generation.r"); source("Simulation.r") })
devtools::load_all("../EIPS", quiet = TRUE)
library(parallel); library(foreach); library(doSNOW)

W_sm <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_full <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2)
h_nu_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1_u2=dat$u1*dat$u2)
inv_link_fn <- function(eta) 1/(1+exp(eta))

ps_spec <- list(
  formula.list = list(r ~ y + u2 + z2),
  h_alpha.list = list(h_full),
  inv_link = inv_link_fn, outcome = "y",
  optimizer = "constrained_nr", cond_threshold = 1e3
)

clean_iqr <- function(x, mult = 2, max_pct = 0.01) {
  Q1 <- quantile(x, 0.25); Q3 <- quantile(x, 0.75); iqr <- Q3 - Q1
  is_out <- x < (Q1 - mult*iqr) | x > (Q3 + mult*iqr)
  max_rm <- floor(max_pct * length(x))
  if (sum(is_out) > max_rm && max_rm > 0) {
    d <- abs(x - median(x)); oi <- which(is_out)
    keep <- oi[order(d[oi], decreasing = TRUE)[(max_rm+1):length(oi)]]
    is_out[keep] <- FALSE
  }
  !is_out
}

n_cores <- min(detectCores() - 1, 10)
n_reps <- 1000
nn <- 2000
mu_true <- get_mu_true("setting1")
cat(sprintf("\nConfig: s1.B1 n=%d, M3 only, cnr/cond=1e3, %d reps  (mu_true=%.4f)\n",
            nn, n_reps, mu_true))

t0 <- proc.time()
cl <- makeCluster(n_cores); registerDoSNOW(cl)
clusterExport(cl, c("nn","ps_spec","W_sm","h_full","h_nu_fn","inv_link_fn"),
              envir = environment())
clusterEvalQ(cl, {
  setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
  devtools::load_all("../EIPS", quiet = TRUE)
  source("Data_Generation.r")
})
pb <- txtProgressBar(max = n_reps, style = 3)
opts <- list(progress = function(n) setTxtProgressBar(pb, n))

res <- foreach(i = 1:n_reps, .combine = 'rbind', .options.snow = opts,
               .packages = c("stringr","Matrix")) %dopar% {
  out <- tryCatch({
    set.seed(i)
    dat <- setting1.B1(nn)
    ebmr <- EIPS$new("y", ps_spec, dat, W_sm)
    ipw <- ebmr$EBMR_IPW(h_nu_fn, model_indices = 1, se.fit = TRUE)
    a <- unname(ebmr$ps_fit.list[[1]]$coefficients)
    c(mu_ipw = ipw$mu_ipw, se_ipw = ipw$se_ipw,
      a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4],
      max_abs_alpha = max(abs(a)))
  }, error = function(e) c(mu_ipw = NA, se_ipw = NA,
                           a1 = NA, a2 = NA, a3 = NA, a4 = NA, max_abs_alpha = NA))
  out
}
close(pb); stopCluster(cl)
elapsed <- (proc.time() - t0)["elapsed"]

mu_v <- res[,"mu_ipw"]; se_v <- res[,"se_ipw"]
valid <- !is.na(mu_v) & !is.na(se_v) & is.finite(mu_v) & is.finite(se_v)
mu_v <- mu_v[valid]; se_v <- se_v[valid]
keep <- clean_iqr(mu_v); mu_c <- mu_v[keep]; se_c <- se_v[keep]
ci_lo <- mu_c - 1.96*se_c; ci_hi <- mu_c + 1.96*se_c
cp <- mean(ci_lo <= mu_true & mu_true <= ci_hi)
alpha_mat <- res[, c("a1","a2","a3","a4")]
n_zero <- sum(rowSums(abs(alpha_mat), na.rm = TRUE) == 0, na.rm = TRUE)
n_sep  <- sum(res[, "max_abs_alpha"] > 50, na.rm = TRUE)
sd_a <- apply(alpha_mat, 2, sd, na.rm = TRUE)
cat(sprintf("\nNA=%d valid=%d Out=%d stuck@0=%d sep(>50)=%d\n",
            sum(!valid), sum(valid), sum(!keep), n_zero, n_sep))
cat(sprintf("Bias=%+.4f ESD=%.4f ESE=%.4f ESE/ESD=%.3f CP=%.3f\n",
            mean(mu_c) - mu_true, sd(mu_c), mean(se_c),
            mean(se_c)/sd(mu_c), cp))
cat(sprintf("sd(alpha) = (%.3f, %.3f, %.3f, %.3f)   max|alpha|=%.1f\n",
            sd_a[1], sd_a[2], sd_a[3], sd_a[4], max(res[,"max_abs_alpha"], na.rm=TRUE)))
cat(sprintf("median(alpha) = (%.3f, %.3f, %.3f, %.3f)\n",
            median(alpha_mat[,1], na.rm=TRUE), median(alpha_mat[,2], na.rm=TRUE),
            median(alpha_mat[,3], na.rm=TRUE), median(alpha_mat[,4], na.rm=TRUE)))
cat(sprintf("\nelapsed: %.0fs\n", elapsed))
saveRDS(list(res = res, ps_spec = ps_spec, n_reps = n_reps, nn = nn,
             mu_true = mu_true, summary = list(
               Bias = mean(mu_c) - mu_true, ESD = sd(mu_c),
               ESE = mean(se_c), ratio = mean(se_c)/sd(mu_c), CP = cp)),
        "Simulation_Test/confirm_s1_n2000_cnr1e3_results.RDS")
cat("\nDONE\n")
