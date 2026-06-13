## Multi-seed multi-size verification: identical() across all seeds/sizes/funcs.
## Saves a summary to compare_data_generation_multiseed.log.

setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

log_file <- "Simulation_Test/verify/compare_data_generation_multiseed.log"
sink(log_file, split = TRUE)

env_orig <- new.env(parent = asNamespace("stats"))
env_test <- new.env(parent = asNamespace("stats"))
sys.source("Data_Generation.r",      envir = env_orig)
sys.source("Data_Generation_test.r", envir = env_test)

shared <- sort(intersect(ls(env_orig), ls(env_test)))
shared <- shared[
  vapply(shared,
         function(nm) is.function(env_orig[[nm]]) && is.function(env_test[[nm]]),
         logical(1))
]

seeds <- c(1L, 12345L, 999999L, 42L, 7L)
sizes <- c(100L, 500L, 2000L)

call_setting <- function(fn, n, seed) {
  set.seed(seed)
  if ("response.rate" %in% names(formals(fn))) {
    fn(n, response.rate = 0.5)
  } else {
    fn(n)
  }
}

cat(sprintf("Testing %d functions x %d sizes x %d seeds = %d total checks\n",
            length(shared), length(sizes), length(seeds),
            length(shared) * length(sizes) * length(seeds)))

n_total <- 0L; n_ok <- 0L; n_diff <- 0L; n_err <- 0L
diffs <- character()

for (seed in seeds) {
  for (n in sizes) {
    for (fn_name in shared) {
      a <- tryCatch(call_setting(env_orig[[fn_name]], n, seed),
                    error = function(e) e)
      b <- tryCatch(call_setting(env_test[[fn_name]], n, seed),
                    error = function(e) e)
      n_total <- n_total + 1L
      if (inherits(a, "error") || inherits(b, "error")) {
        n_err <- n_err + 1L
        diffs <- c(diffs,
                   sprintf("ERR  %-25s seed=%d n=%d", fn_name, seed, n))
      } else if (identical(a, b)) {
        n_ok <- n_ok + 1L
      } else {
        n_diff <- n_diff + 1L
        diffs <- c(diffs,
                   sprintf("DIFF %-25s seed=%d n=%d", fn_name, seed, n))
      }
    }
  }
}

cat(sprintf("\nResults: OK=%d  DIFF=%d  ERR=%d\n", n_ok, n_diff, n_err))
if (length(diffs)) {
  cat("\nFailures:\n")
  for (line in diffs) cat("  ", line, "\n", sep = "")
} else {
  cat("\nAll checks passed -- byte-equivalence confirmed across all seeds.\n")
}

sink()
