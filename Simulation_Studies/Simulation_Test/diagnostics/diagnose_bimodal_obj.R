## Investigate: are both basins true CUE minima, or is one a spurious fixed point?
## Compare CUE objective Q(nu) = G(nu)'W(nu)G(nu) at both solutions
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")
suppressMessages({
  source("Basic_setup.r"); source("Data_Generation.r")
  source("config/scenarios.R"); source("Simulation.r")
})
devtools::load_all("../EPS", quiet = TRUE)

ps_spec_base <- get_ps_spec("9")
h_alpha_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2)
W_sm <- function(g.matrix) solve(t(g.matrix) %*% g.matrix / nrow(g.matrix))
h_nu_fn <- function(dat) cbind(u1=dat$u1, u2=dat$u2, z1=dat$z1, z2=dat$z2, u1u2=dat$u1*dat$u2)
compute_pi <- function(eta) 1/(1+exp(eta))

data_file <- misspecified_model_all_data_file.list[["setting4"]][["miss50"]]
all_data <- NULL
for (fi in seq_along(data_file)) {
  if (file.exists(data_file[[fi]])) {
    test <- readRDS(data_file[[fi]])
    if (nrow(test) / 1000 >= 2) { all_data <- test; break }
  }
}

nn <- 2000
ps_spec_m23 <- list(
  formula.list = ps_spec_base$formula.list[2:3],
  h_alpha.list = list(h_alpha_fn, h_alpha_fn),
  inv_link = ps_spec_base$inv_link,
  outcome = ps_spec_base$outcome
)

cat("=== CUE objective at both nu basins ===\n\n")
cat(sprintf("%-5s | %-35s | %-35s | %s\n",
    "Rep", "M2-dom (init=0.99,0.01)", "M3-dom (init=0.01,0.99)", "Lower"))
cat(strrep("-", 85), "\n")

for (i in 1:20) {
  dat <- all_data[((i-1)*nn + 1):(i*nn), ]
  ebmr <- EPS$new("y", ps_spec_m23, dat, W_sm)
  ps_mat <- do.call(cbind, lapply(ebmr$ps_fit.list, function(pf) pf$fitted.values))
  r_vec <- as.vector(dat$r)
  h_x <- cbind(1, h_nu_fn(dat))

  # CUE objective for a given nu
  cue_obj <- function(nu) {
    w <- nu^2 / sum(nu^2)
    ps_nu <- as.vector(ps_mat %*% w)
    g_mat <- as.vector(r_vec / ps_nu - 1) * h_x
    G <- colMeans(g_mat)
    W <- tryCatch(solve(t(g_mat) %*% g_mat / nn), error = function(e) diag(ncol(h_x)))
    as.numeric(t(G) %*% W %*% G)
  }

  # Run from both inits
  ipw_m2 <- ebmr$EBMR_IPW(h_nu_fn, model_indices = 1:2, nu_init = c(0.99, 0.01), se.fit = FALSE)
  ipw_m3 <- ebmr$EBMR_IPW(h_nu_fn, model_indices = 1:2, nu_init = c(0.01, 0.99), se.fit = FALSE)

  obj_m2 <- cue_obj(ipw_m2$nu.hat)
  obj_m3 <- cue_obj(ipw_m3$nu.hat)

  lower <- if (abs(obj_m2 - obj_m3) < 1e-6) "TIE"
           else if (obj_m2 < obj_m3) "M2"
           else "M3"

  cat(sprintf("%5d | w=(%.3f,%.3f) Q=%.4e | w=(%.3f,%.3f) Q=%.4e | %s\n",
      i,
      ipw_m2$w.hat[1], ipw_m2$w.hat[2], obj_m2,
      ipw_m3$w.hat[1], ipw_m3$w.hat[2], obj_m3,
      lower))
}
