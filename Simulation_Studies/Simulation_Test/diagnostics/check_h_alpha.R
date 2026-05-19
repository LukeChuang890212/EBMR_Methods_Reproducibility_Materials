setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
source("Basic_setup.r")
source("config/scenarios.R")

ps_spec <- get_ps_spec("9-alt1")
cat("Classes of h_alpha.list elements:\n")
for (i in seq_along(ps_spec[["h_alpha.list"]])) {
  cat(sprintf("  h_alpha.list[[%d]]: %s\n", i, class(ps_spec[["h_alpha.list"]][[i]])))
}

cat("\nUsing indices 1,3:\n")
ps_spec_13 <- list(
  formula.list = ps_spec[["formula.list"]][c(1, 3)],
  h_alpha.list = ps_spec[["h_alpha.list"]][c(1, 3)],
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]]
)
for (i in seq_along(ps_spec_13[["h_alpha.list"]])) {
  cat(sprintf("  ps_spec_13 h_alpha.list[[%d]]: %s\n", i, class(ps_spec_13[["h_alpha.list"]][[i]])))
  if (is.function(ps_spec_13[["h_alpha.list"]][[i]])) {
    cat("    -> is a function\n")
  } else if (is.character(ps_spec_13[["h_alpha.list"]][[i]])) {
    cat(sprintf("    -> character: %s\n", paste(ps_spec_13[["h_alpha.list"]][[i]], collapse=", ")))
  } else {
    cat(sprintf("    -> other type: %s\n", typeof(ps_spec_13[["h_alpha.list"]][[i]])))
    str(ps_spec_13[["h_alpha.list"]][[i]])
  }
}
