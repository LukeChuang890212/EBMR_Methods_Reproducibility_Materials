## End-to-end diagnostic: does the nu:23 override fire when Simulation_demo.r
## calls run_scenario("9-3", setting="setting2", version="test65", type="HT") ?
##
## Strategy: re-source everything fresh (mirroring how Simulation_demo.r boots),
## then run a 20-rep mini-version that exercises the exact code path. The script
## reports at each layer whether the expected behavior is observed.
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

cat("===== Layer 1: installed package has new args? =====\n")
library(EBMRalgorithmFast4)
fmls <- names(formals(EBMRalgorithmFast4:::EBMR_IPW))
ok1 <- all(c("nu_optimizer", "nu_cond_threshold") %in% fmls)
cat(sprintf("EBMR_IPW formals include nu_optimizer/nu_cond_threshold: %s\n", ok1))
cat(sprintf("  installed at: %s\n", find.package("EBMRalgorithmFast4")))
if (!ok1) stop("Package is stale -- run devtools::install('EBMRalgorithmFast4') first.")

cat("\n===== Layer 2: Simulation.r passes nu params through to EBMR_IPW? =====\n")
source("Simulation.r")
sim_fmls <- names(formals(simulate))
ok2 <- all(c("nu_optimizer", "nu_cond_threshold") %in% sim_fmls)
cat(sprintf("simulate() formals include nu params: %s\n", ok2))
sim_body <- deparse(body(simulate))
ok2b <- any(grepl("nu_optimizer\\s*=\\s*nu_optimizer", sim_body))
cat(sprintf("simulate() body forwards nu_optimizer to EBMR_IPW: %s\n", ok2b))
if (!ok2 || !ok2b) stop("Simulation.r is stale -- re-source it (or restart R).")

cat("\n===== Layer 3: scenarios.R resolver recognizes nu:23 ? =====\n")
suppressMessages({ source("Basic_setup.r"); source("Data_Generation.r")
                   source("config/scenarios.R") })
sc_src <- readLines("config/scenarios.R")
ok3 <- any(grepl("startsWith\\(tolower\\(trimws\\(o_model_field\\)\\), \"nu\"\\)", sc_src))
cat(sprintf("scenarios.R contains nu resolver branch: %s\n", ok3))
ok3b <- any(grepl("nu_optimizer\\s*=\\s*nu_optimizer_override", sc_src))
cat(sprintf("scenarios.R forwards nu_optimizer_override to simulate(): %s\n", ok3b))
if (!ok3 || !ok3b) stop("scenarios.R is stale -- re-source it.")

cat("\n===== Layer 4: optimizer_overrides.txt has the entry? =====\n")
ov <- readLines("config/optimizer_overrides.txt")
ov <- ov[!grepl("^#|^\\s*$", ov)]
nu23_lines <- ov[grepl("nu:23", ov)]
cat(sprintf("Lines matching 'nu:23': %d\n", length(nu23_lines)))
for (l in nu23_lines) cat("  ", l, "\n")
ok4 <- length(nu23_lines) >= 1
if (!ok4) stop("No nu:23 entries in optimizer_overrides.txt.")

cat("\n===== Layer 5: 20-rep run via simulate() with the same args run_scenario builds =====\n")
# Mirror run_scenario's _23 path for 9-3 setting2 miss50 n=2000
scenario_id <- "9-3"; setting <- "setting2"; miss_rate <- "miss50"; current_n <- 2000
ps_spec <- get_ps_spec(SCENARIOS[[scenario_id]]$ps_spec_id)
model_set <- c(2, 3); J_sub <- length(model_set)
subset_ps_spec <- list(
  formula.list = ps_spec$formula.list[model_set],
  h_alpha.list = ps_spec$h_alpha.list[model_set],
  inv_link = ps_spec$inv_link, outcome = ps_spec$outcome,
  optimizer = ps_spec$optimizer, cond_threshold = ps_spec$cond_threshold
)
# Resolve overrides as scenarios.R does (alpha + nu)
nu_opt_o <- "L-BFGS-B"; nu_cond_o <- 1e8; has_nu <- FALSE
opt_list <- vector("list", J_sub); cond_list <- vector("list", J_sub); has_alpha <- FALSE
for (ol in ov) {
  f <- trimws(strsplit(ol, "\\|")[[1]])
  if (length(f) < 7) next
  if (!(f[1]==scenario_id && f[2]==setting && f[3]==miss_rate && as.numeric(f[4])==current_n)) next
  if (startsWith(tolower(trimws(f[5])), "nu")) {
    sfx <- sub("^nu:?", "", tolower(trimws(f[5])))
    tgt <- if (grepl(",", sfx)) sort(as.numeric(trimws(strsplit(sfx, ",")[[1]])))
           else sort(as.numeric(strsplit(sfx, "")[[1]]))
    if (J_sub > 1 && identical(sort(model_set), tgt)) {
      nu_opt_o <- f[6]; nu_cond_o <- as.numeric(f[7]); has_nu <- TRUE
    }
  } else {
    o_mi <- as.numeric(trimws(strsplit(f[5], ",")[[1]]))
    for (k in seq_along(model_set)) {
      if (model_set[k] %in% o_mi) {
        opt_list[[k]] <- f[6]; cond_list[[k]] <- as.numeric(f[7]); has_alpha <- TRUE
      }
    }
  }
}
if (has_alpha) {
  for (k in seq_along(model_set)) {
    if (is.null(opt_list[[k]])) opt_list[[k]] <- "L-BFGS-B"
    if (is.null(cond_list[[k]])) cond_list[[k]] <- 1e8
  }
  subset_ps_spec$optimizer <- opt_list
  subset_ps_spec$cond_threshold <- cond_list
}
cat(sprintf("Resolved subset_ps_spec optimizer: %s\n", paste(sapply(subset_ps_spec$optimizer, as.character), collapse=", ")))
cat(sprintf("Resolved subset_ps_spec cond:      %s\n", paste(sapply(subset_ps_spec$cond_threshold, as.character), collapse=", ")))
cat(sprintf("Resolved nu_optimizer:             %s\n", nu_opt_o))
cat(sprintf("Resolved nu_cond_threshold:        %g\n", nu_cond_o))
if (!has_nu) cat("  *** WARNING: nu override NOT picked up; check optimizer_overrides.txt entry. ***\n")

# Now actually call simulate() with 20 reps, no saved file
all_data <- readRDS("Simulation_Data/setting2.B1_n2000_replicate1000.RDS")
alpha.true <- misspecified_model_alpha.true.list$setting2$miss50[[1]]
ps_model.true <- function(dat, alpha.true) {
  X <- cbind(rep(1, nrow(dat)), dat$y, dat$u1, dat$u2)
  1 / (1 + exp(X %*% alpha.true))
}
cat("\nRunning simulate() with 20 reps (no save)...\n")
res <- simulate(
  all_data = all_data, ps_model.true = ps_model.true, alpha.true = alpha.true,
  ps_specifications = subset_ps_spec, n = current_n, replicate_num = 20,
  save_file = NULL, type = "HT", setting = setting,
  nu_optimizer = nu_opt_o, nu_cond_threshold = nu_cond_o
)
mu_v <- res["mu_ipw", ]
w1_v <- res["w.hat1", ]
if (is.null(w1_v)) w1_v <- res[grep("^w.hat", rownames(res))[1], ]
cat(sprintf("\n20-rep summary: mean(mu_ipw)=%.4f  sd=%.4f  range=[%.4f, %.4f]\n",
            mean(mu_v, na.rm=T), sd(mu_v, na.rm=T), min(mu_v, na.rm=T), max(mu_v, na.rm=T)))
cat(sprintf("w_M2 distribution: <0.1 n=%d, 0.1-0.9 n=%d, >=0.9 n=%d (mean=%.3f)\n",
            sum(w1_v < 0.1, na.rm=T), sum(w1_v >= 0.1 & w1_v < 0.9, na.rm=T),
            sum(w1_v >= 0.9, na.rm=T), mean(w1_v, na.rm=T)))
cat(sprintf("\nWith cnr/1e2 nu fix, expect: w_M2 >= 0.9 in ~19/20 reps.\n"))
cat(sprintf("Without fix (default):       w_M2 bimodal -- some near 0, some near 1.\n"))
cat("\nDONE\n")
