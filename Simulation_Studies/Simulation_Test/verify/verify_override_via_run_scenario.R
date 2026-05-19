## End-to-end verification: call run_scenario() exactly as Simulation_demo.R does
## for scenario 9-3 setting2 test65, then inspect the resulting _23 file's M3 alpha
## distribution to confirm the cnr/1e4 override actually took effect.
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EBMRalgorithmFast4", quiet = TRUE)

cat("\n##### Step 1: call run_scenario as Simulation_demo.R does #####\n")
run_scenario("9-3", setting = "setting2", version = "test65", type = "HT")

cat("\n##### Step 2: inspect the produced _23 file's M3 alpha distribution #####\n")
f <- "Simulation_Results/EBMR_IPW_setting2-miss50-scenario9-3_23_n2000_replicate1000_test65.RDS"
if (!file.exists(f)) {
  cat("FILE NOT FOUND:", f, "\n"); quit(status = 1)
}
x <- readRDS(f)
cat(sprintf("dim: %s    n_reps: %d\n", paste(dim(x), collapse=" x "), ncol(x)))
cat("rownames:\n"); print(rownames(x))

# For _23 the output is 24 rows: 4 (mu/se) + 4 (M2 alpha) + 4 (M3 alpha) + 4 (M2 se_alpha) + 4 (M3 se_alpha) + 2 nu + 2 w
# Identify M2 vs M3 alpha rows (M2 is slot 1, M3 is slot 2 in the subset)
alpha_rows <- grep("^alpha.hat", rownames(x))
cat(sprintf("\nalpha rows: %d\n", length(alpha_rows)))
# First 4 = M2 (slot 1), next 4 = M3 (slot 2)
a_M2 <- x[alpha_rows[1:4], , drop = FALSE]
a_M3 <- x[alpha_rows[5:8], , drop = FALSE]
max_a_M2 <- apply(a_M2, 2, function(v) max(abs(v), na.rm = TRUE))
max_a_M3 <- apply(a_M3, 2, function(v) max(abs(v), na.rm = TRUE))
cat(sprintf("\nM2 (formula u1_z2): max|alpha| median=%.2f q90=%.2f q99=%.2f max=%.2f\n",
            median(max_a_M2), quantile(max_a_M2, 0.9), quantile(max_a_M2, 0.99), max(max_a_M2)))
cat(sprintf("M3 (formula u2_z2): max|alpha| median=%.2f q90=%.2f q99=%.2f max=%.2f\n",
            median(max_a_M3), quantile(max_a_M3, 0.9), quantile(max_a_M3, 0.99), max(max_a_M3)))
cat("\nExpected if override APPLIED to M3: max|alpha| <= ~50 across all reps (cnr/1e4 caps it).\n")
cat("Expected if override NOT applied:   M3 max|alpha| up to ~250 in some reps (default L-BFGS-B).\n")

# ESE/ESD via clean_sim_result
mu_true <- get_mu_true("setting2")
cleaned <- clean_sim_result(x, multiplier = 2, max_pct = 0.01, verbose = FALSE)
sc <- cleaned$result
mu_v <- sc["mu_ipw", ]; se_v <- sc["se_ipw", ]
cat(sprintf("\nESE/ESD summary  (n_clean=%d Out=%d):\n", cleaned$n_successful, cleaned$n_outliers))
cat(sprintf("  Bias=%+.4f  ESD=%.4f  ESE(mean)=%.4f  ESE(median)=%.4f  ratio=%.3f\n",
            mean(mu_v)-mu_true, sd(mu_v), mean(se_v), median(se_v), mean(se_v)/sd(mu_v)))
cat("\nExpected with override: ratio close to 1 (was 0.18 in the pre-override saved test65).\n")
