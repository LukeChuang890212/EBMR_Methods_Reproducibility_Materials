## Verify that the modified caption code in Simulation.r produces the expected
## LaTeX text for scenarios 7-1, 8-1, 9-3, and falls back correctly for others.
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({ source("Basic_setup.r"); source("Data_Generation.r"); source("Simulation.r") })

# Extract just the caption-building snippet (mirrors lines 619-656) for testing
test_caption <- function(setting, scenario, alpha_true.list) {
  mu.true <- if (setting == "setting1") 2.000 else if (setting == "setting2") 0.804 else 0.7
  mu_str <- sprintf("%.3f", round(mu.true, 3))
  alpha50 <- paste(sprintf("%.2f", round(alpha_true.list[[setting]][[1]][[1]], 2)), collapse = ", ")
  alpha30 <- paste(sprintf("%.2f", round(alpha_true.list[[setting]][[2]][[1]], 2)), collapse = ", ")

  outcome_type <- if (setting %in% c("setting1", "setting3")) "continuous" else "binary"
  model_desc <- if (scenario == "7-1") {
    paste0("a correctly specified model with ",
           "$\\pi_2(\\bm{X}, Y; \\bm{\\alpha}_2) = \\mathrm{expit}\\{(1, U_1, Z_1, Y)\\bm{\\alpha}_2\\}$ and ",
           "$\\pi_3(\\bm{X}, Y; \\bm{\\alpha}_3) = \\mathrm{expit}\\{(1, U_2, Z_1, Y)\\bm{\\alpha}_3\\}$.")
  } else if (scenario == "8-1") {
    paste0("a correctly specified model with ",
           "$\\pi_2(\\bm{X}, Y; \\bm{\\alpha}_2) = \\mathrm{expit}\\{(1, U_1, Z_1, Y)\\bm{\\alpha}_2\\}$ and ",
           "$\\pi_3(\\bm{X}, Y; \\bm{\\alpha}_3) = \\mathrm{expit}\\{(1, U_2, Z_2, Y)\\bm{\\alpha}_3\\}$.")
  } else if (scenario == "9-3") {
    paste0("a locally misspecified model with ",
           "$\\bm{h}_{\\bm{\\alpha}_1}(\\bm{X}) = (1, U_1, U_2, Z_1, Z_2, U_2^2)$.")
  } else {
    "a correctly specified model."
  }
  paste0(
    "Estimation of $\\mu_0 = ", mu_str, "$ is conducted under ", outcome_type, " outcome $Y$, ",
    "where the set of candidate models includes ", model_desc, " ",
    "The parameter vectors $\\bm{\\alpha}_0 = (", alpha50, ")^{\\top}$ and ",
    "$(", alpha30, ")^{\\top}$ correspond to missingness rates of 50\\% and 30\\%, respectively."
  )
}

cases <- list(
  list(scenario = "7-1", setting = "setting1", at = correct_model_alpha.true.list),
  list(scenario = "8-1", setting = "setting1", at = correct_model_alpha.true.list),
  list(scenario = "9-3", setting = "setting2", at = misspecified_model_alpha.true.list),
  list(scenario = "9-2", setting = "setting2", at = misspecified_model_alpha.true.list),
  list(scenario = "9-1", setting = "setting1", at = correct_model_alpha.true.list)
)

for (c in cases) {
  cat(strrep("=", 80), "\n", sep="")
  cat(sprintf("scenario=%s  setting=%s\n", c$scenario, c$setting))
  cat(strrep("-", 80), "\n", sep="")
  cat(test_caption(c$setting, c$scenario, c$at), "\n")
}
