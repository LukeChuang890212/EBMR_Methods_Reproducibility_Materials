## Check whether the recently regenerated 9-3 _23 file used the nu cnr/1e2 override.
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({ source("Basic_setup.r"); source("Simulation.r") })

for (sc in c("9-2", "9-3")) {
  f <- sprintf("Simulation_Results/EBMR_IPW_setting2-miss50-scenario%s_23_n2000_replicate1000_test65.RDS", sc)
  cat(sprintf("\n========== %s ==========\n", basename(f)))
  if (!file.exists(f)) { cat("MISSING\n"); next }
  x <- readRDS(f)
  cat(sprintf("dim: %d x %d (rows x reps)\n", nrow(x), ncol(x)))
  cat(sprintf("rownames: %s\n", paste(rownames(x), collapse=", ")))
  # mu_ipw is row 1; se_ipw is row 3; w.hat rows are at the end (last 2 for J=2)
  w_rows <- grep("^w\\.hat", rownames(x))
  cat(sprintf("w.hat rows: %s\n", paste(rownames(x)[w_rows], collapse=", ")))
  if (length(w_rows) >= 1) {
    w1_v <- x[w_rows[1], ]
    valid <- !is.na(w1_v) & is.finite(w1_v)
    w1_v <- w1_v[valid]
    cat(sprintf("w_M2 distribution: <0.1 n=%d, 0.1-0.9 n=%d, >=0.9 n=%d (mean=%.3f, sd=%.3f, n_valid=%d)\n",
                sum(w1_v < 0.1), sum(w1_v >= 0.1 & w1_v < 0.9), sum(w1_v >= 0.9),
                mean(w1_v), sd(w1_v), length(w1_v)))
  }
  # Clean & report ratio
  cl <- clean_sim_result(x, multiplier = 2, max_pct = 0.01, verbose = FALSE)
  sc_clean <- cl$result
  mu_v <- sc_clean[1, ]; se_v <- sc_clean[3, ]
  mu_true <- 0.8047
  cat(sprintf("[cleaned] n=%d  Bias=%+.4f  ESD=%.4f  ESE(mean)=%.4f  ratio=%.3f\n",
              ncol(sc_clean), mean(mu_v)-mu_true, sd(mu_v), mean(se_v), mean(se_v)/sd(mu_v)))
  cat("Expected if override took effect: w_M2 >= 0.9 in ~99% of reps, ratio ~ 0.95\n")
  cat("Expected if override DID NOT fire: bimodal w_M2 (some near 0, many near 1), ratio ~ 0.18\n")
}
