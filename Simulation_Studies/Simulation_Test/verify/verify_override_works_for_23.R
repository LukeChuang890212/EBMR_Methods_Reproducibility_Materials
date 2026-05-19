## Verify that when run_scenario builds the _23 (M2+M3) ps_spec for
## scenario=9-3, setting=setting2, miss=miss50, n=2000, the M3 slot gets
## constrained_nr / cond=1e4 from optimizer_overrides.txt — and run 20 reps
## to confirm M3's alpha really stays small (max|alpha| ~50 not ~250).
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EPS", quiet = TRUE)

# Mirror exactly what run_scenario does for _23:
scenario_id <- "9-3"; setting <- "setting2"; miss_rate <- "miss50"; current_n <- 2000
ps_spec <- get_ps_spec(SCENARIOS[[scenario_id]]$ps_spec_id)  # 9-alt1
model_set <- c(2, 3)        # _23 = M2 + M3 in absolute indices
J_sub <- length(model_set)
subset_ps_spec <- list(
  formula.list = ps_spec$formula.list[model_set],
  h_alpha.list = ps_spec$h_alpha.list[model_set],
  inv_link = ps_spec$inv_link, outcome = ps_spec$outcome,
  optimizer = ps_spec$optimizer, cond_threshold = ps_spec$cond_threshold
)

# Replicate the override resolution (lines 879-927 of scenarios.R)
override_file <- "config/optimizer_overrides.txt"
override_lines <- readLines(override_file)
override_lines <- override_lines[!grepl("^#|^\\s*$", override_lines)]
opt_list <- vector("list", J_sub); cond_list <- vector("list", J_sub)
has_override <- FALSE
for (ol in override_lines) {
  fields <- trimws(strsplit(ol, "\\|")[[1]])
  if (length(fields) >= 7) {
    o_scenario <- fields[1]; o_setting <- fields[2]
    o_miss <- fields[3]; o_n <- as.numeric(fields[4])
    o_model_idx <- as.numeric(trimws(strsplit(fields[5], ",")[[1]]))
    o_optimizer <- fields[6]; o_cond <- as.numeric(fields[7])
    if (o_scenario == scenario_id && o_setting == setting &&
        o_miss == miss_rate && o_n == current_n) {
      for (k in seq_along(model_set)) {
        if (model_set[k] %in% o_model_idx) {
          opt_list[[k]] <- o_optimizer
          cond_list[[k]] <- o_cond
          has_override <- TRUE
        }
      }
    }
  }
}
for (k in seq_along(model_set)) {
  if (is.null(opt_list[[k]])) opt_list[[k]] <- "L-BFGS-B"
  if (is.null(cond_list[[k]])) cond_list[[k]] <- 1e8
}
subset_ps_spec$optimizer <- opt_list
subset_ps_spec$cond_threshold <- cond_list

cat("===== Resolved ps_spec for scenario 9-3 setting2 miss50 n=2000 _23 =====\n")
cat("model_set =", model_set, "  (M2 -> slot 1, M3 -> slot 2)\n")
cat("formulas:\n"); for (j in 1:J_sub) cat("  slot", j, ": "); print(subset_ps_spec$formula.list)
cat("optimizer:\n"); print(subset_ps_spec$optimizer)
cat("cond_threshold:\n"); print(subset_ps_spec$cond_threshold)
cat("\nExpected: slot 1 (M2) = L-BFGS-B/1e8, slot 2 (M3) = constrained_nr/1e4\n\n")

# Now run 20 reps and inspect M3's alpha distribution
all_data <- readRDS("Simulation_Data/setting2.B1_n2000_replicate1000.RDS")
nn <- 2000; n_check <- 20
W_sm <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))

max_a_M2 <- numeric(n_check); max_a_M3 <- numeric(n_check)
for (i in 1:n_check) {
  dat <- all_data[((i-1)*nn+1):(i*nn), ]
  ebmr <- EPS$new("y", subset_ps_spec, dat, W_sm)
  a_M2 <- unname(ebmr$ps_fit.list[[1]]$coefficients)
  a_M3 <- unname(ebmr$ps_fit.list[[2]]$coefficients)
  max_a_M2[i] <- max(abs(a_M2))
  max_a_M3[i] <- max(abs(a_M3))
}
cat(sprintf("Across %d reps:\n", n_check))
cat(sprintf("  M2 max|alpha|: median=%.2f  q90=%.2f  max=%.2f\n",
            median(max_a_M2), quantile(max_a_M2, 0.9), max(max_a_M2)))
cat(sprintf("  M3 max|alpha|: median=%.2f  q90=%.2f  max=%.2f\n",
            median(max_a_M3), quantile(max_a_M3, 0.9), max(max_a_M3)))
cat("\nExpected: M3 max|alpha| stays in low single digits (cond=1e4 rejects bigger steps).\n")
cat("If the override is NOT applied, M3 max|alpha| would be much larger (>50, up to ~250).\n")
