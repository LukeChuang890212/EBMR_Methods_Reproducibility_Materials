## Verify the nu:<subset> resolver in scenarios.R: nu:23 should match _23 only.
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

# Replicate the resolver block standalone, exercising it for different subsets.
test_resolver <- function(scenario_id, setting, miss_rate, current_n, model_set,
                         override_file = "config/optimizer_overrides.txt") {
  J_sub <- length(model_set)
  nu_optimizer_override <- "L-BFGS-B"
  nu_cond_override <- 1e8
  has_nu_override <- FALSE

  override_lines <- readLines(override_file)
  override_lines <- override_lines[!grepl("^#|^\\s*$", override_lines)]
  for (ol in override_lines) {
    fields <- trimws(strsplit(ol, "\\|")[[1]])
    if (length(fields) < 7) next
    o_scenario <- fields[1]; o_setting <- fields[2]
    o_miss <- fields[3]; o_n <- as.numeric(fields[4])
    o_model_field <- fields[5]
    o_optimizer <- fields[6]; o_cond <- as.numeric(fields[7])
    if (!(o_scenario == scenario_id && o_setting == setting &&
          o_miss == miss_rate && o_n == current_n)) next
    if (startsWith(tolower(trimws(o_model_field)), "nu")) {
      suffix <- sub("^nu:?", "", tolower(trimws(o_model_field)))
      if (nzchar(suffix)) {
        target_set <- if (grepl(",", suffix)) {
          sort(as.numeric(trimws(strsplit(suffix, ",")[[1]])))
        } else {
          sort(as.numeric(strsplit(suffix, "")[[1]]))
        }
        if (J_sub > 1 && identical(sort(model_set), target_set)) {
          nu_optimizer_override <- o_optimizer
          nu_cond_override <- o_cond
          has_nu_override <- TRUE
        }
      }
    }
  }
  list(nu_opt = nu_optimizer_override, nu_cond = nu_cond_override, hit = has_nu_override)
}

cases <- list(
  list(s="9-3", set="setting2", m="miss50", n=2000, ms=c(2,3),     expect_hit=TRUE),
  list(s="9-3", set="setting2", m="miss50", n=2000, ms=c(1,2),     expect_hit=FALSE),
  list(s="9-3", set="setting2", m="miss50", n=2000, ms=c(1,3),     expect_hit=FALSE),
  list(s="9-3", set="setting2", m="miss50", n=2000, ms=c(1,2,3),   expect_hit=FALSE),
  list(s="9-3", set="setting2", m="miss50", n=2000, ms=c(2),       expect_hit=FALSE), # J_sub=1
  list(s="9-2", set="setting2", m="miss50", n=2000, ms=c(2,3),     expect_hit=TRUE),
  list(s="9-4", set="setting2", m="miss50", n=2000, ms=c(2,3),     expect_hit=TRUE),
  list(s="8-2", set="setting2", m="miss50", n=2000, ms=c(2,3),     expect_hit=FALSE), # different scenario
  list(s="9-3", set="setting2", m="miss50", n=500,  ms=c(2,3),     expect_hit=FALSE), # different n
  list(s="9-3", set="setting2", m="miss30", n=2000, ms=c(2,3),     expect_hit=FALSE)  # different miss
)

cat(sprintf("%-25s | %-4s | %4s | %-7s | %s | %s\n",
            "scenario,set,miss,n,subset", "J_sub", "want", "got", "nu_opt", "nu_cond"))
cat(strrep("-", 95), "\n", sep="")
ok <- TRUE
for (c1 in cases) {
  r <- test_resolver(c1$s, c1$set, c1$m, c1$n, c1$ms)
  pass <- identical(r$hit, c1$expect_hit)
  ok <- ok && pass
  cat(sprintf("%-25s | %4d | %4s | %-7s | %-14s | %g %s\n",
              sprintf("%s,%s,%s,%d,_%s", c1$s, c1$set, c1$m, c1$n, paste(c1$ms, collapse="")),
              length(c1$ms),
              ifelse(c1$expect_hit, "HIT", "MISS"),
              ifelse(r$hit,        "HIT", "MISS"),
              r$nu_opt, r$nu_cond,
              ifelse(pass, "[OK]", "[FAIL]")))
}
cat(ifelse(ok, "\nALL PASS\n", "\nSOME FAIL\n"))
