setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({source("Basic_setup.r"); source("Data_Generation.r"); source("config/scenarios.R"); source("Simulation.r"); library(EBMRalgorithmFast4)})

ps_spec <- get_ps_spec("9-alt1")
ps_spec_1 <- list(formula.list=ps_spec[["formula.list"]][1], h_alpha.list=ps_spec[["h_alpha.list"]][1],
                  inv_link=ps_spec[["inv_link"]], outcome=ps_spec[["outcome"]])
W_func  <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_nu_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1_u2=dat$u1*dat$u2)
all_data <- readRDS("Simulation_Data/Setting3.B1_n2000_replicate1000.RDS")
n_val <- 2000

Rprof("profile2.out", interval=0.005)
for (i in 1:100) {
  dat  <- all_data[((i-1)*n_val+1):(i*n_val), ]
  ebmr <- EBMRAlgorithmFast4$new("y", ps_spec_1, dat, W_func)
  ebmr$EBMR_IPW(h_nu=h_nu_fn, type="HT", se.fit=TRUE)
}
Rprof(NULL)

p <- summaryRprof("profile2.out")
cat("=== Top 20 by self time ===\n"); print(head(p$by.self, 20))
cat("\n=== Top 20 by total time ===\n"); print(head(p$by.total, 20))
