setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

source("Basic_setup.r")
source("Data_Generation.r")
source("config/scenarios.R")
source("Simulation.r")
library(EBMRalgorithmFast4)

data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]][[1]]
all_data <- readRDS(data_file)
n_val <- 2000

ps_spec <- get_ps_spec("9-alt1")
subset_ps_spec <- list(
  formula.list = ps_spec[["formula.list"]][c(1, 3)],
  h_alpha.list = ps_spec[["h_alpha.list"]][c(1, 3)],
  inv_link = ps_spec[["inv_link"]],
  outcome = ps_spec[["outcome"]]
)

W_func <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
inv_link_func <- ps_spec[["inv_link"]]

# Run on more reps to get the full picture
n_reps <- 200
cat("=== MODEL COMPARISON: obj_M1 vs obj_M2 and ensemble weight ===\n\n")
cat(sprintf("%-5s %10s %10s %8s %8s %8s %8s %8s %8s\n",
    "Rep", "obj_M1", "obj_M2", "ratio", "w1", "||a1||", "||a2||", "mu_ipw", "type"))

results <- data.frame()
for (rep_i in 1:n_reps) {
  dat <- all_data[((rep_i - 1) * n_val + 1):(rep_i * n_val), ]

  tryCatch({
    ebmr <- EBMRAlgorithmFast4[["new"]]("y", subset_ps_spec, dat, W_func)
    res <- ebmr[["EBMR_IPW"]](
      h_nu = function(dat) cbind(u1 = dat[["u1"]], u2 = dat[["u2"]], z1 = dat[["z1"]], z2 = dat[["z2"]]),
      se.fit = FALSE, type = "HT"
    )

    fit1 <- ebmr[["ps_fit.list"]][[1]]
    fit2 <- ebmr[["ps_fit.list"]][[2]]
    obj1 <- fit1[["opt"]][["objective"]]
    obj2 <- fit2[["opt"]][["objective"]]
    a1_norm <- sqrt(sum(fit1[["coefficients"]]^2))
    a2_norm <- sqrt(sum(fit2[["coefficients"]]^2))
    w1 <- res[["w.hat"]][1]

    wrong <- w1 < 0.5
    cat(sprintf("%-5d %10.6f %10.6f %8.2f %8.4f %8.3f %8.3f %8.4f %s\n",
        rep_i, obj1, obj2, obj1/obj2, w1, a1_norm, a2_norm,
        res[["mu_ipw"]], ifelse(wrong, "WRONG", "")))

    results <- rbind(results, data.frame(
      rep = rep_i, obj1 = obj1, obj2 = obj2, ratio = obj1/obj2,
      w1 = w1, a1_norm = a1_norm, a2_norm = a2_norm,
      mu_ipw = res[["mu_ipw"]], wrong = wrong
    ))
  }, error = function(e) {
    cat(sprintf("%-5d ERROR: %s\n", rep_i, conditionMessage(e)))
  })
}

cat(sprintf("\n=== SUMMARY (n=%d) ===\n", nrow(results)))
cat(sprintf("Wrong weight reps: %d / %d\n", sum(results[["wrong"]]), nrow(results)))

cat("\n--- Normal reps ---\n")
normal <- results[!results[["wrong"]], ]
cat(sprintf("  obj1: mean=%.6f, median=%.6f\n", mean(normal[["obj1"]]), median(normal[["obj1"]])))
cat(sprintf("  obj2: mean=%.6f, median=%.6f\n", mean(normal[["obj2"]]), median(normal[["obj2"]])))
cat(sprintf("  ratio (obj1/obj2): mean=%.2f, median=%.2f\n", mean(normal[["ratio"]]), median(normal[["ratio"]])))
cat(sprintf("  ||a2||: mean=%.2f, median=%.2f, max=%.2f\n", mean(normal[["a2_norm"]]), median(normal[["a2_norm"]]), max(normal[["a2_norm"]])))

cat("\n--- Wrong weight reps ---\n")
wrong_r <- results[results[["wrong"]], ]
if (nrow(wrong_r) > 0) {
  cat(sprintf("  obj1: mean=%.6f, median=%.6f\n", mean(wrong_r[["obj1"]]), median(wrong_r[["obj1"]])))
  cat(sprintf("  obj2: mean=%.6f, median=%.6f\n", mean(wrong_r[["obj2"]]), median(wrong_r[["obj2"]])))
  cat(sprintf("  ratio (obj1/obj2): mean=%.2f, median=%.2f\n", mean(wrong_r[["ratio"]]), median(wrong_r[["ratio"]])))
  cat(sprintf("  ||a2||: mean=%.2f, median=%.2f, max=%.2f\n", mean(wrong_r[["a2_norm"]]), median(wrong_r[["a2_norm"]]), max(wrong_r[["a2_norm"]])))
}

# Key question: is obj2 lower in wrong reps, or is obj1 higher?
cat("\n--- Discriminant analysis ---\n")
cat(sprintf("  Normal: mean obj1=%.6f, mean obj2=%.6f\n", mean(normal[["obj1"]]), mean(normal[["obj2"]])))
if (nrow(wrong_r) > 0)
  cat(sprintf("  Wrong:  mean obj1=%.6f, mean obj2=%.6f\n", mean(wrong_r[["obj1"]]), mean(wrong_r[["obj2"]])))
