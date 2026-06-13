## compare_data_generation.R
##
## Verify that Data_Generation_test.r produces byte-identical data frames to
## Data_Generation.r for every shared setting function.
##
## Strategy:
##   1. Source Data_Generation.r into a clean environment -> capture functions.
##   2. Source Data_Generation_test.r into a clean environment -> capture too.
##   3. For each function name common to both environments, run with the same
##      RNG seed and compare results with identical() (byte-level equality).
##   4. Print one summary line per setting; print a brief diagnostic for any
##      mismatch.
##
## Run from anywhere; the script sets its own working directory.

setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

# Capture original functions in an isolated environment so they don't leak
# into the test versions (or vice versa) via the global env. Parent is the
# loaded stats namespace so rbinom/rnorm/etc. resolve.
env_orig  <- new.env(parent = asNamespace("stats"))
env_test  <- new.env(parent = asNamespace("stats"))
sys.source("Data_Generation.r",      envir = env_orig)
sys.source("Data_Generation_test.r", envir = env_test)

# Shared function names (those defined in BOTH environments).
names_orig  <- ls(env_orig)
names_test  <- ls(env_test)
shared      <- sort(intersect(names_orig, names_test))
only_orig   <- setdiff(names_orig, names_test)
only_test   <- setdiff(names_test, names_orig)

if (length(only_orig) > 0L) {
  cat("Functions defined ONLY in Data_Generation.r (not yet in test):\n  ",
      paste(only_orig, collapse = ", "), "\n\n", sep = "")
}
if (length(only_test) > 0L) {
  cat("Functions defined ONLY in Data_Generation_test.r (extra helpers etc.):\n  ",
      paste(only_test, collapse = ", "), "\n\n", sep = "")
}

# Compare for a few sample sizes; n=2000 is the workhorse; n=500 exercises the
# small-n branches of B-variant alpha0 = ifelse(n >= 1000, ...).
sample_sizes <- c(500L, 2000L)
seed_used    <- 12345L

# Helper: invoke a setting function with the seed pre-set, picking a sensible
# argument list depending on whether the function declares 'response.rate'.
call_setting <- function(fn, n) {
  set.seed(seed_used)
  if ("response.rate" %in% names(formals(fn))) {
    fn(n, response.rate = 0.5)
  } else {
    fn(n)
  }
}

# For each (fn, n), capture orig vs test and compare with identical().
# Mismatches get a short diagnostic line; column-wise max abs diff for numerics.
result_rows <- list()
for (fn_name in shared) {
  fn_o <- env_orig[[fn_name]]
  fn_t <- env_test[[fn_name]]
  if (!is.function(fn_o) || !is.function(fn_t)) next

  for (n in sample_sizes) {
    out_o <- tryCatch(call_setting(fn_o, n), error = function(e) e)
    out_t <- tryCatch(call_setting(fn_t, n), error = function(e) e)

    err_o <- inherits(out_o, "error")
    err_t <- inherits(out_t, "error")

    if (err_o || err_t) {
      status <- "ERR  "
      detail <- sprintf("orig_err=%s test_err=%s",
                        if (err_o) conditionMessage(out_o) else "-",
                        if (err_t) conditionMessage(out_t) else "-")
    } else if (identical(out_o, out_t)) {
      status <- "OK   "
      detail <- sprintf("nrow=%d ncol=%d", nrow(out_o), ncol(out_o))
    } else {
      status <- "DIFF "
      # Build a column-by-column diff summary
      cols_o <- names(out_o); cols_t <- names(out_t)
      if (!identical(cols_o, cols_t)) {
        detail <- sprintf("colnames differ: orig=[%s] test=[%s]",
                          paste(cols_o, collapse = ","),
                          paste(cols_t, collapse = ","))
      } else if (nrow(out_o) != nrow(out_t)) {
        detail <- sprintf("nrow differs: orig=%d test=%d", nrow(out_o), nrow(out_t))
      } else {
        diffs <- sapply(cols_o, function(c) {
          a <- out_o[[c]]; b <- out_t[[c]]
          if (is.numeric(a) && is.numeric(b)) {
            max(abs(a - b), na.rm = TRUE)
          } else {
            sum(a != b, na.rm = TRUE)
          }
        })
        worst <- which.max(diffs)
        detail <- sprintf("max|orig-test| per col: %s  (worst: %s = %.6g)",
                          paste(sprintf("%s=%.3g", cols_o, diffs), collapse = ", "),
                          cols_o[worst], diffs[worst])
      }
    }

    result_rows[[length(result_rows) + 1L]] <-
      data.frame(fn = fn_name, n = n, status = status, detail = detail,
                 stringsAsFactors = FALSE)
  }
}

results <- do.call(rbind, result_rows)
n_ok   <- sum(results$status == "OK   ")
n_diff <- sum(results$status == "DIFF ")
n_err  <- sum(results$status == "ERR  ")

cat(strrep("=", 80), "\n", sep = "")
cat(sprintf("compare_data_generation.R   seed=%d   %d functions x %d sizes\n",
            seed_used, length(shared), length(sample_sizes)))
cat(sprintf("  OK:   %3d\n", n_ok))
cat(sprintf("  DIFF: %3d\n", n_diff))
cat(sprintf("  ERR:  %3d\n", n_err))
cat(strrep("=", 80), "\n", sep = "")

# Print full table -- one row per (fn, n).
for (i in seq_len(nrow(results))) {
  cat(sprintf("  %-25s n=%-5d %s %s\n",
              results$fn[i], results$n[i], results$status[i], results$detail[i]))
}

# Quick top-line return value for programmatic use.
invisible(list(results = results, summary = c(OK = n_ok, DIFF = n_diff, ERR = n_err)))
