#------------------------------------------------------------------------------#
# Data_Generation_test.r — refactored data generators (test version)
#
# Goal: produce IDENTICAL data frames to Data_Generation.r for any given RNG seed.
# The redundant per-setting boilerplate is factored into a small number of helpers,
# and each setting function becomes a short, declarative composition.
#
# RNG-faithfulness contract: every rbinom/rnorm/sample call happens in the
# same order as the original, so set.seed(X); setting_name(n) produces the
# same data frame whether sourced from this file or the original.
#
# Convention shared with Data_Generation.r:
#   A suffix = correctly-specified PS    (A1=50% miss, A2=30% miss)
#   B suffix = misspecified PS via local exp-tilt   (B1=50% miss, B2=30% miss)
#   Cho_*    = Cho 2025 reference settings (uses x1,x2,x3 covariates)
#
# To verify equivalence:
#   source("Data_Generation.r");           orig <- as.list(.GlobalEnv)
#   source("Data_Generation_test.r");      refac <- as.list(.GlobalEnv)
#   then compare orig[[name]](n) vs refac[[name]](n) under the same set.seed().
# A ready-made comparison script lives at
#   Simulation_Test/verify/compare_data_generation.R
#------------------------------------------------------------------------------#

#==============================================================================
# Section 1. Low-level helpers
#==============================================================================

#' Generate (z1, z2, u1, u2) covariates with z1 ~ Bernoulli, z2 ~ N, u1 ~ Bernoulli,
#' u2 ~ N. Returns a named list to be unpacked by the caller. RNG calls happen
#' in the canonical order z1 -> z2 -> u1 -> u2.
gen_zu <- function(n, pz1, m_z2, sd_z2, pu1, m_u2, sd_u2) {
  z1 <- rbinom(n, size = 1, prob = pz1)
  z2 <- rnorm(n, mean = m_z2, sd = sd_z2)
  u1 <- rbinom(n, size = 1, prob = pu1)
  u2 <- rnorm(n, mean = m_u2, sd = sd_u2)
  list(z1 = z1, z2 = z2, u1 = u1, u2 = u2)
}

#' Same as gen_zu but with z1 ~ Normal (used in settings 11 and 12).
gen_zNu <- function(n, m_z1, sd_z1, m_z2, sd_z2, pu1, m_u2, sd_u2) {
  z1 <- rnorm(n, mean = m_z1, sd = sd_z1)
  z2 <- rnorm(n, mean = m_z2, sd = sd_z2)
  u1 <- rbinom(n, size = 1, prob = pu1)
  u2 <- rnorm(n, mean = m_u2, sd = sd_u2)
  list(z1 = z1, z2 = z2, u1 = u1, u2 = u2)
}

#' Cho-style covariates (x1, x2, x3). Default sd matches Cho_M1.A1/M2.A1 family
#' (sqrt(1/3)); pass sd_x3 to override the third covariate's sd.
gen_x123 <- function(n, mean = 1, sd = sqrt(1/3), sd_x3 = sd) {
  x1 <- rnorm(n, mean = mean, sd = sd)
  x2 <- rnorm(n, mean = mean, sd = sd)
  x3 <- rnorm(n, mean = mean, sd = sd_x3)
  list(x1 = x1, x2 = x2, x3 = x3)
}

#' Continuous outcome: y ~ N(mu_x, sd).
gen_y_normal <- function(mu_x, sd = 1) {
  rnorm(length(mu_x), mean = mu_x, sd = sd)
}

#' Binary outcome via standard logistic on linear predictor mu_x.
gen_y_binary <- function(mu_x) {
  prob <- exp(mu_x) / (1 + exp(mu_x))
  rbinom(length(mu_x), size = 1, prob = prob)
}

#' Complement-logistic propensity used in settings 1-19:
#'   P(r=1) = 1 / (1 + exp(alpha0 + ay*y + au1*u1 + au2*u2))
ps_complement <- function(alpha0, ay, au1, au2, y, u1, u2) {
  1 / (1 + exp(alpha0 + ay * y + au1 * u1 + au2 * u2))
}

#' Standard-logistic propensity (Cho settings): P(r=1) = expit(eta).
ps_standard <- function(eta) {
  exp(eta) / (1 + exp(eta))
}

#' Apply local-misspecification tilt to a propensity:
#'   ps' = ps * exp(mult * n^{-1/2} * tilt_value),  capped at 1 -> 0.95.
apply_tilt <- function(ps, n, tilt_value, mult = 1, cap = 0.95) {
  out <- ps * exp(mult * n^(-1/2) * tilt_value)
  out[out > 1] <- cap
  out
}

#' Sample r ~ Bernoulli(propensity).
sample_r <- function(propensity) {
  rbinom(length(propensity), size = 1, prob = propensity)
}

#' Build the standard (z1, z2, u1, u2, y, r) data frame.
df_zu <- function(cov, y, r) {
  data.frame(z1 = cov$z1, z2 = cov$z2, u1 = cov$u1, u2 = cov$u2, y = y, r = r)
}

#==============================================================================
# Section 2. Settings 1-4: cont. or binary Y, simple m(x) = 1 + u1+u2+z1+z2
# Setting1 (cont. Y, sd_z2=1), setting2 (binary Y, sd_z2=1),
# setting3 (cont. Y, sd_z2=2), setting4 (binary Y, sd_z2=2 with different alpha).
#==============================================================================

# --- setting1 (continuous y, ay=0.2 / au1=-0.8 / au2=+0.8) -------------------
setting1.A1 <- function(n) {
  cov <- gen_zu(n, pz1 = 0.4, m_z2 = 0, sd_z2 = 1, pu1 = 0.6, m_u2 = 0, sd_u2 = 1)
  mu_x <- 1 + cov$u1 + cov$u2 + cov$z1 + cov$z2
  y <- gen_y_normal(mu_x)
  ps <- ps_complement(0.0803, 0.2, -0.8, 0.8, y, cov$u1, cov$u2)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting1.A2 <- function(n) {
  cov <- gen_zu(n, 0.4, 0, 1, 0.6, 0, 1)
  mu_x <- 1 + cov$u1 + cov$u2 + cov$z1 + cov$z2
  y <- gen_y_normal(mu_x)
  ps <- ps_complement(-0.9646, 0.2, -0.8, 0.8, y, cov$u1, cov$u2)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting1.B1 <- function(n) {
  alpha0 <- ifelse(n >= 1000, 0.1765, 0.2730)
  cov <- gen_zu(n, 0.4, 0, 1, 0.6, 0, 1)
  mu_x <- 1 + cov$u1 + cov$u2 + cov$z1 + cov$z2
  y <- gen_y_normal(mu_x)
  ps <- ps_complement(alpha0, 0.2, -0.8, 0.8, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, y + cov$u1 + cov$u2, mult = 1)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting1.B2 <- function(n) {
  alpha0 <- ifelse(n >= 1000, -0.7843, -0.6167)
  cov <- gen_zu(n, 0.4, 0, 1, 0.6, 0, 1)
  mu_x <- 1 + cov$u1 + cov$u2 + cov$z1 + cov$z2
  y <- gen_y_normal(mu_x)
  ps <- ps_complement(alpha0, 0.2, -0.8, 0.8, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, y + cov$u1 + cov$u2, mult = 1)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

# --- setting2 (binary y; same m_x; PS coefs same as setting1; tilt mult = 10) -
setting2.A1 <- function(n) {
  cov <- gen_zu(n, 0.4, 0, 1, 0.6, 0, 1)
  mu_x <- 1 + cov$u1 + cov$u2 + cov$z1 + cov$z2
  y <- gen_y_binary(mu_x)
  ps <- ps_complement(0.3193, 0.2, -0.8, 0.8, y, cov$u1, cov$u2)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting2.A2 <- function(n) {
  cov <- gen_zu(n, 0.4, 0, 1, 0.6, 0, 1)
  mu_x <- 1 + cov$u1 + cov$u2 + cov$z1 + cov$z2
  y <- gen_y_binary(mu_x)
  ps <- ps_complement(-0.6734, 0.2, -0.8, 0.8, y, cov$u1, cov$u2)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting2.B1 <- function(n) {
  alpha0 <- ifelse(n >= 1000, 0.8794, 1.4141)
  cov <- gen_zu(n, 0.4, 0, 1, 0.6, 0, 1)
  mu_x <- 1 + cov$u1 + cov$u2 + cov$z1 + cov$z2
  y <- gen_y_binary(mu_x)
  ps <- ps_complement(alpha0, 0.2, -0.8, 0.8, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, y + cov$u1 + cov$u2, mult = 10)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting2.B2 <- function(n) {
  alpha0 <- ifelse(n >= 1000, 0.2012, 0.7290)
  cov <- gen_zu(n, 0.4, 0, 1, 0.6, 0, 1)
  mu_x <- 1 + cov$u1 + cov$u2 + cov$z1 + cov$z2
  y <- gen_y_binary(mu_x)
  ps <- ps_complement(alpha0, 0.2, -0.8, 0.8, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, y + cov$u1 + cov$u2, mult = 10)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

# --- setting3 (continuous y; sd_z2=2; ay=0.2/au1=-0.4/au2=+1.0) ---------------
setting3.A1 <- function(n) {
  cov <- gen_zu(n, 0.4, 0, 2, 0.6, 0, 1)
  mu_x <- 1 + cov$u1 + cov$u2 + cov$z1 + cov$z2
  y <- gen_y_normal(mu_x)
  ps <- ps_complement(-0.1581, 0.2, -0.4, 1.0, y, cov$u1, cov$u2)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting3.A2 <- function(n) {
  cov <- gen_zu(n, 0.4, 0, 2, 0.6, 0, 1)
  mu_x <- 1 + cov$u1 + cov$u2 + cov$z1 + cov$z2
  y <- gen_y_normal(mu_x)
  ps <- ps_complement(-1.2697, 0.2, -0.4, 1.0, y, cov$u1, cov$u2)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting3.B1 <- function(n) {
  alpha0 <- ifelse(n >= 1000, -0.0774, 0.0070)
  cov <- gen_zu(n, 0.4, 0, 2, 0.6, 0, 1)
  mu_x <- 1 + cov$u1 + cov$u2 + cov$z1 + cov$z2
  y <- gen_y_normal(mu_x)
  ps <- ps_complement(alpha0, 0.2, -0.4, 1.0, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, y + cov$u1 + cov$u2, mult = 1)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting3.B2 <- function(n) {
  alpha0 <- ifelse(n >= 1000, -1.0998, -0.9435)
  cov <- gen_zu(n, 0.4, 0, 2, 0.6, 0, 1)
  mu_x <- 1 + cov$u1 + cov$u2 + cov$z1 + cov$z2
  y <- gen_y_normal(mu_x)
  ps <- ps_complement(alpha0, 0.2, -0.4, 1.0, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, y + cov$u1 + cov$u2, mult = 1)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

# --- setting4 (binary y; sd_z2=2; ay=0.15/au1=-0.30/au2=+0.80) ---------------
setting4.A1 <- function(n) {
  cov <- gen_zu(n, 0.4, 0, 2, 0.6, 0, 1)
  mu_x <- 1 + cov$u1 + cov$u2 + cov$z1 + cov$z2
  y <- gen_y_binary(mu_x)
  ps <- ps_complement(-0.5173, 0.15, -0.30, 0.80, y, cov$u1, cov$u2)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting4.A2 <- function(n) {
  cov <- gen_zu(n, 0.4, 0, 2, 0.6, 0, 1)
  mu_x <- 1 + cov$u1 + cov$u2 + cov$z1 + cov$z2
  y <- gen_y_binary(mu_x)
  ps <- ps_complement(-1.4592, 0.15, -0.30, 0.80, y, cov$u1, cov$u2)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting4.B1 <- function(n) {
  alpha0 <- ifelse(n >= 1000, 0.0153, 0.5347)
  cov <- gen_zu(n, 0.4, 0, 2, 0.6, 0, 1)
  mu_x <- 1 + cov$u1 + cov$u2 + cov$z1 + cov$z2
  y <- gen_y_binary(mu_x)
  ps <- ps_complement(alpha0, 0.15, -0.30, 0.80, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, y + cov$u1 + cov$u2, mult = 10)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting4.B2 <- function(n) {
  alpha0 <- ifelse(n >= 1000, -0.6297, -0.0120)
  cov <- gen_zu(n, 0.4, 0, 2, 0.6, 0, 1)
  mu_x <- 1 + cov$u1 + cov$u2 + cov$z1 + cov$z2
  y <- gen_y_binary(mu_x)
  ps <- ps_complement(alpha0, 0.15, -0.30, 0.80, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, y + cov$u1 + cov$u2, mult = 10)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

#==============================================================================
# Section 3. Settings 8-10: small-coef m(x), unconventional PS coefficients
# Shared covariates (z1~Bern(0.3), u1~Bern(0.7), u2~N(0,sd)) with sd varying.
# Note: original setting8/9/10 .A1 declare a 'response.rate' arg that goes unused.
# Tilt expression in B settings is -y+u1-u2 (different sign pattern).
#==============================================================================

# --- setting8 (binary y; sd_u2=3) ---
setting8.A1 <- function(n, response.rate) {
  cov <- gen_zu(n, 0.3, 0, 1, 0.7, 0, 3)
  mu_x <- 0.2 + 0.5*cov$z1 + 0.5*cov$z2 + 0.5*cov$u1 + 0.5*cov$u2
  y <- gen_y_binary(mu_x)
  # Hardcoded form: 1/(1+exp(-0.04 + 0.2*y + 0.2*u1)) -> ay=0.2, au1=0.2, au2=0.
  ps <- ps_complement(-0.04, 0.2, 0.2, 0, y, cov$u1, cov$u2)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting8.A2 <- function(n) {
  cov <- gen_zu(n, 0.3, 0, 1, 0.7, 0, 3)
  mu_x <- 0.2 + 0.5*cov$z1 + 0.5*cov$z2 + 0.5*cov$u1 + 0.5*cov$u2
  y <- gen_y_binary(mu_x)
  ps <- ps_complement(-0.9, 0.2, 0.2, 0, y, cov$u1, cov$u2)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting8.B1 <- function(n) {
  cov <- gen_zu(n, 0.3, 0, 1, 0.7, 0, 3)
  mu_x <- 0.2 + 0.5*cov$z1 + 0.5*cov$z2 + 0.5*cov$u1 + 0.5*cov$u2
  y <- gen_y_binary(mu_x)
  ps <- ps_complement(-0.04, 0.2, 0.2, 0, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, -y + cov$u1 - cov$u2, mult = 1)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting8.B2 <- function(n, response.rate) {
  cov <- gen_zu(n, 0.3, 0, 1, 0.7, 0, 3)
  mu_x <- 0.2 + 0.5*cov$z1 + 0.5*cov$z2 + 0.5*cov$u1 + 0.5*cov$u2
  y <- gen_y_binary(mu_x)
  ps <- ps_complement(-0.9, 0.2, 0.2, 0, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, -y + cov$u1 - cov$u2, mult = 1)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

# --- setting9 (continuous y; sd_u2=2) ---
setting9.A1 <- function(n) {
  cov <- gen_zu(n, 0.3, 0, 1, 0.7, 0, 2)
  mu_x <- 0.2 + 0.5*cov$z1 + 0.5*cov$z2 + 0.5*cov$u1 + 0.5*cov$u2
  y <- gen_y_normal(mu_x)
  ps <- ps_complement(0.6, 0.2, -1, -0.5, y, cov$u1, cov$u2)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting9.A2 <- function(n) {
  cov <- gen_zu(n, 0.3, 0, 1, 0.7, 0, 2)
  mu_x <- 0.2 + 0.5*cov$z1 + 0.5*cov$z2 + 0.5*cov$u1 + 0.5*cov$u2
  y <- gen_y_normal(mu_x)
  ps <- ps_complement(-0.5, 0.2, -1, -0.5, y, cov$u1, cov$u2)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting9.B1 <- function(n) {
  cov <- gen_zu(n, 0.3, 0, 1, 0.7, 0, 2)
  mu_x <- 0.2 + 0.5*cov$z1 + 0.5*cov$z2 + 0.5*cov$u1 + 0.5*cov$u2
  y <- gen_y_normal(mu_x)
  ps <- ps_complement(0.6, 0.2, -1, -0.5, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, -y + cov$u1 - cov$u2, mult = 1)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting9.B2 <- function(n) {
  cov <- gen_zu(n, 0.3, 0, 1, 0.7, 0, 2)
  mu_x <- 0.2 + 0.5*cov$z1 + 0.5*cov$z2 + 0.5*cov$u1 + 0.5*cov$u2
  y <- gen_y_normal(mu_x)
  ps <- ps_complement(-0.5, 0.2, -1, -0.5, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, -y + cov$u1 - cov$u2, mult = 1)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

# --- setting10 (binary y; sd_u2=2; bigger z-loadings in m(x)) ---
setting10.A1 <- function(n, response.rate) {
  cov <- gen_zu(n, 0.3, 0, 1, 0.7, 0, 2)
  mu_x <- 0.2 + cov$z1 + cov$z2 + 0.5*cov$u1 + 0.5*cov$u2
  y <- gen_y_binary(mu_x)
  ps <- ps_complement(0.4, 0.4, -1, -0.5, y, cov$u1, cov$u2)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting10.A2 <- function(n) {
  cov <- gen_zu(n, 0.3, 0, 1, 0.7, 0, 2)
  mu_x <- 0.2 + cov$z1 + cov$z2 + 0.5*cov$u1 + 0.5*cov$u2
  y <- gen_y_binary(mu_x)
  ps <- ps_complement(-0.6, 0.4, -1, -0.5, y, cov$u1, cov$u2)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting10.B1 <- function(n) {
  cov <- gen_zu(n, 0.3, 0, 1, 0.7, 0, 2)
  mu_x <- 0.2 + cov$z1 + cov$z2 + 0.5*cov$u1 + 0.5*cov$u2
  y <- gen_y_binary(mu_x)
  ps <- ps_complement(0.4, 0.4, -1, -0.5, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, -y + cov$u1 - cov$u2, mult = 1)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting10.B2 <- function(n, response.rate) {
  cov <- gen_zu(n, 0.3, 0, 1, 0.7, 0, 2)
  mu_x <- 0.2 + cov$z1 + cov$z2 + 0.5*cov$u1 + 0.5*cov$u2
  y <- gen_y_binary(mu_x)
  ps <- ps_complement(-0.6, 0.4, -1, -0.5, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, -y + cov$u1 - cov$u2, mult = 1)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

#==============================================================================
# Section 4. Settings 11-12: z1 is NORMAL (not Bernoulli); tilt uses -0.5*y
#==============================================================================

# --- setting11 (continuous y) ---
setting11.A1 <- function(n) {
  cov <- gen_zNu(n, 0, 1.5, 0, 1, 0.7, 0, 2)
  mu_x <- 0.2 + 0.3*cov$z1 + 0.5*cov$z2 + 0.5*cov$u1 + 0.5*cov$u2
  y <- gen_y_normal(mu_x)
  ps <- ps_complement(1.3, 0.2, -2, -0.5, y, cov$u1, cov$u2)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting11.A2 <- function(n) {
  cov <- gen_zNu(n, 0, 1.5, 0, 1, 0.7, 0, 2)
  mu_x <- 0.2 + 0.3*cov$z1 + 0.5*cov$z2 + 0.5*cov$u1 + 0.5*cov$u2
  y <- gen_y_normal(mu_x)
  ps <- ps_complement(0.2, 0.2, -2, -0.5, y, cov$u1, cov$u2)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting11.B1 <- function(n) {
  cov <- gen_zNu(n, 0, 1.5, 0, 1, 0.7, 0, 2)
  mu_x <- 0.2 + 0.3*cov$z1 + 0.5*cov$z2 + 0.5*cov$u1 + 0.5*cov$u2
  y <- gen_y_normal(mu_x)
  ps <- ps_complement(1.3, 0.2, -2, -0.5, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, -0.5*y + cov$u1 - cov$u2, mult = 1)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting11.B2 <- function(n) {
  cov <- gen_zNu(n, 0, 1.5, 0, 1, 0.7, 0, 2)
  mu_x <- 0.2 + 0.3*cov$z1 + 0.5*cov$z2 + 0.5*cov$u1 + 0.5*cov$u2
  y <- gen_y_normal(mu_x)
  ps <- ps_complement(0.2, 0.2, -2, -0.5, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, -0.5*y + cov$u1 - cov$u2, mult = 1)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

# --- setting12 (binary y; m_x uses z2 with coef=1) ---
setting12.A1 <- function(n, response.rate) {
  cov <- gen_zNu(n, 0, 1.5, 0, 1, 0.7, 0, 2)
  mu_x <- 0.2 + 0.5*cov$z1 + cov$z2 + 0.5*cov$u1 + 0.5*cov$u2
  y <- gen_y_binary(mu_x)
  ps <- ps_complement(1.2, 0.4, -2, -0.5, y, cov$u1, cov$u2)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting12.A2 <- function(n) {
  cov <- gen_zNu(n, 0, 1.5, 0, 1, 0.7, 0, 2)
  mu_x <- 0.2 + 0.5*cov$z1 + cov$z2 + 0.5*cov$u1 + 0.5*cov$u2
  y <- gen_y_binary(mu_x)
  ps <- ps_complement(0.05, 0.4, -2, -0.5, y, cov$u1, cov$u2)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting12.B1 <- function(n) {
  cov <- gen_zNu(n, 0, 1.5, 0, 1, 0.7, 0, 2)
  mu_x <- 0.2 + 0.5*cov$z1 + cov$z2 + 0.5*cov$u1 + 0.5*cov$u2
  y <- gen_y_binary(mu_x)
  ps <- ps_complement(1.2, 0.4, -2, -0.5, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, -0.5*y + cov$u1 - cov$u2, mult = 1)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting12.B2 <- function(n, response.rate) {
  cov <- gen_zNu(n, 0, 1.5, 0, 1, 0.7, 0, 2)
  mu_x <- 0.2 + 0.5*cov$z1 + cov$z2 + 0.5*cov$u1 + 0.5*cov$u2
  y <- gen_y_binary(mu_x)
  ps <- ps_complement(0.05, 0.4, -2, -0.5, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, -0.5*y + cov$u1 - cov$u2, mult = 1)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

#==============================================================================
# Section 5. Settings 13-14: m_x = 0.2+0.3*u1+0.3*u2+0.6*z1+0.6*z2; u2~N(1,2)
#==============================================================================

# --- setting13 (continuous y; au1=au2=-0.4) ---
setting13.A1 <- function(n) {
  cov <- gen_zu(n, 0.4, 0, 2, 0.6, 1, 2)
  mu_x <- 0.2 + 0.3*cov$u1 + 0.3*cov$u2 + 0.6*cov$z1 + 0.6*cov$z2
  y <- gen_y_normal(mu_x)
  ps <- ps_complement(0.1161, 0.2, -0.4, -0.4, y, cov$u1, cov$u2)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting13.A2 <- function(n) {
  cov <- gen_zu(n, 0.4, 0, 2, 0.6, 1, 2)
  mu_x <- 0.2 + 0.3*cov$u1 + 0.3*cov$u2 + 0.6*cov$z1 + 0.6*cov$z2
  y <- gen_y_normal(mu_x)
  ps <- ps_complement(-0.7780, 0.2, -0.4, -0.4, y, cov$u1, cov$u2)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting13.B1 <- function(n) {
  alpha0 <- ifelse(n >= 1000, 0.1764, 0.2396)
  cov <- gen_zu(n, 0.4, 0, 2, 0.6, 1, 2)
  mu_x <- 0.2 + 0.3*cov$u1 + 0.3*cov$u2 + 0.6*cov$z1 + 0.6*cov$z2
  y <- gen_y_normal(mu_x)
  ps <- ps_complement(alpha0, 0.2, -0.4, -0.4, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, y + cov$u1 + cov$u2, mult = 1)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting13.B2 <- function(n) {
  alpha0 <- ifelse(n >= 1000, -0.6807, -0.5863)
  cov <- gen_zu(n, 0.4, 0, 2, 0.6, 1, 2)
  mu_x <- 0.2 + 0.3*cov$u1 + 0.3*cov$u2 + 0.6*cov$z1 + 0.6*cov$z2
  y <- gen_y_normal(mu_x)
  ps <- ps_complement(alpha0, 0.2, -0.4, -0.4, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, y + cov$u1 + cov$u2, mult = 1)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

# --- setting14 (binary y; alpha0 in A1/A2 is hardcoded in PS; B1/B2 use n-dep) -
# Note: original A1/A2 define a local 'alpha0' via ifelse but the hardcoded value
# in response.prob is what gets used. We honour that: ignore the dead alpha0.
setting14.A1 <- function(n) {
  cov <- gen_zu(n, 0.4, 0, 2, 0.6, 1, 2)
  mu_x <- 0.2 + 0.3*cov$u1 + 0.3*cov$u2 + 0.6*cov$z1 + 0.6*cov$z2
  y <- gen_y_binary(mu_x)
  ps <- ps_complement(0.1154, 0.2, -0.4, -0.4, y, cov$u1, cov$u2)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting14.A2 <- function(n) {
  cov <- gen_zu(n, 0.4, 0, 2, 0.6, 1, 2)
  mu_x <- 0.2 + 0.3*cov$u1 + 0.3*cov$u2 + 0.6*cov$z1 + 0.6*cov$z2
  y <- gen_y_binary(mu_x)
  ps <- ps_complement(-0.7696, 0.2, -0.4, -0.4, y, cov$u1, cov$u2)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting14.B1 <- function(n) {
  alpha0 <- ifelse(n >= 1000, -0.6653, -0.5688)
  cov <- gen_zu(n, 0.4, 0, 2, 0.6, 1, 2)
  mu_x <- 0.2 + 0.3*cov$u1 + 0.3*cov$u2 + 0.6*cov$z1 + 0.6*cov$z2
  y <- gen_y_binary(mu_x)
  ps <- ps_complement(alpha0, 0.2, -0.4, -0.4, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, y + cov$u1 + cov$u2, mult = 1)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting14.B2 <- function(n) {
  alpha0 <- ifelse(n >= 1000, -0.6653, -0.5688)
  cov <- gen_zu(n, 0.4, 0, 2, 0.6, 1, 2)
  mu_x <- 0.2 + 0.3*cov$u1 + 0.3*cov$u2 + 0.6*cov$z1 + 0.6*cov$z2
  y <- gen_y_binary(mu_x)
  ps <- ps_complement(alpha0, 0.2, -0.4, -0.4, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, y + cov$u1 + cov$u2, mult = 1)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

#==============================================================================
# Section 6. Setting 15: binary y with nonlinear m(x) (u2^2, z2^2);
# m_x = 0.2+0.2*u1+0.2*u2+0.4*z1+0.4*z2+0.4*u2^2+0.4*z2^2; u2~N(3,2)
# B variants use PS = 1/(1+exp(a0+0.6*y+0.6*u1+0.6*u2)), tilt mult=10
#==============================================================================

setting15.A1 <- function(n) {
  cov <- gen_zu(n, 0.4, 0, 1, 0.6, 3, 2)
  mu_x <- 0.2 + 0.2*cov$u1 + 0.2*cov$u2 + 0.4*cov$z1 + 0.4*cov$z2 +
          0.4*cov$u2^2 + 0.4*cov$z2^2
  y <- gen_y_binary(mu_x)
  ps <- ps_complement(0.0917, 0.2, -0.4, -0.4, y, cov$u1, cov$u2)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting15.A2 <- function(n) {
  cov <- gen_zu(n, 0.4, 0, 1, 0.6, 3, 2)
  mu_x <- 0.2 + 0.2*cov$u1 + 0.2*cov$u2 + 0.4*cov$z1 + 0.4*cov$z2 +
          0.4*cov$u2^2 + 0.4*cov$z2^2
  y <- gen_y_binary(mu_x)
  ps <- ps_complement(-0.7946, 0.2, -0.4, -0.4, y, cov$u1, cov$u2)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting15.B1 <- function(n) {
  alpha0 <- ifelse(n >= 1000, -0.3480, 0.0973)
  cov <- gen_zu(n, 0.4, 0, 1, 0.6, 3, 2)
  mu_x <- 0.2 + 0.2*cov$u1 + 0.2*cov$u2 + 0.4*cov$z1 + 0.4*cov$z2 +
          0.4*cov$u2^2 + 0.4*cov$z2^2
  y <- gen_y_binary(mu_x)
  ps <- ps_complement(alpha0, 0.6, 0.6, 0.6, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, y + cov$u1 + cov$u2, mult = 10)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting15.B2 <- function(n) {
  alpha0 <- ifelse(n >= 1000, -1.0052, -0.4191)
  cov <- gen_zu(n, 0.4, 0, 1, 0.6, 3, 2)
  mu_x <- 0.2 + 0.2*cov$u1 + 0.2*cov$u2 + 0.4*cov$z1 + 0.4*cov$z2 +
          0.4*cov$u2^2 + 0.4*cov$z2^2
  y <- gen_y_binary(mu_x)
  ps <- ps_complement(alpha0, 0.6, 0.6, 0.6, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, y + cov$u1 + cov$u2, mult = 10)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

#==============================================================================
# Section 7. Settings 16-19 (B variants only): different tilt expressions and
# different multipliers. Shared m_x = 0.2+0.6*z1+0.6*z2+0.3*u1+0.3*u2;
# shared base PS = 1/(1+exp(a0+0.2*y-0.4*u1-0.4*u2)); u2 ~ N(3, 2).
# 16 (cont. y) and 17 (binary y) use mult=10, tilt=(y+u1+u2).
# 18 (cont. y) and 19 (binary y) use mult=1,  tilt=(u1+u2+z1+z2).
#==============================================================================

setting16.B1 <- function(n) {
  alpha0 <- ifelse(n >= 1000, 0.7340, 1.2000)
  cov <- gen_zu(n, 0.4, 0, 2, 0.6, 3, 2)
  mu_x <- 0.2 + 0.6*cov$z1 + 0.6*cov$z2 + 0.3*cov$u1 + 0.3*cov$u2
  y <- gen_y_normal(mu_x)
  ps <- ps_complement(alpha0, 0.2, -0.4, -0.4, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, y + cov$u1 + cov$u2, mult = 10)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting16.B2 <- function(n) {
  alpha0 <- ifelse(n >= 1000, -0.1460, 0.0560)
  cov <- gen_zu(n, 0.4, 0, 2, 0.6, 3, 2)
  mu_x <- 0.2 + 0.6*cov$z1 + 0.6*cov$z2 + 0.3*cov$u1 + 0.3*cov$u2
  y <- gen_y_normal(mu_x)
  ps <- ps_complement(alpha0, 0.2, -0.4, -0.4, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, y + cov$u1 + cov$u2, mult = 10)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting17.B1 <- function(n) {
  alpha0 <- ifelse(n >= 1000, 0.7370, 1.1980)
  cov <- gen_zu(n, 0.4, 0, 2, 0.6, 3, 2)
  mu_x <- 0.2 + 0.6*cov$z1 + 0.6*cov$z2 + 0.3*cov$u1 + 0.3*cov$u2
  y <- gen_y_binary(mu_x)
  ps <- ps_complement(alpha0, 0.2, -0.4, -0.4, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, y + cov$u1 + cov$u2, mult = 10)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting17.B2 <- function(n) {
  alpha0 <- ifelse(n >= 1000, -0.0730, 0.2730)
  cov <- gen_zu(n, 0.4, 0, 2, 0.6, 3, 2)
  mu_x <- 0.2 + 0.6*cov$z1 + 0.6*cov$z2 + 0.3*cov$u1 + 0.3*cov$u2
  y <- gen_y_binary(mu_x)
  ps <- ps_complement(alpha0, 0.2, -0.4, -0.4, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, y + cov$u1 + cov$u2, mult = 10)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting18.B1 <- function(n) {
  cov <- gen_zu(n, 0.4, 0, 2, 0.6, 3, 2)
  mu_x <- 0.2 + 0.6*cov$z1 + 0.6*cov$z2 + 0.3*cov$u1 + 0.3*cov$u2
  y <- gen_y_normal(mu_x)
  ps <- ps_complement(0.1900, 0.2, -0.4, -0.4, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, cov$u1 + cov$u2 + cov$z1 + cov$z2, mult = 1)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting18.B2 <- function(n) {
  cov <- gen_zu(n, 0.4, 0, 2, 0.6, 3, 2)
  mu_x <- 0.2 + 0.6*cov$z1 + 0.6*cov$z2 + 0.3*cov$u1 + 0.3*cov$u2
  y <- gen_y_normal(mu_x)
  ps <- ps_complement(-0.6600, 0.2, -0.4, -0.4, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, cov$u1 + cov$u2 + cov$z1 + cov$z2, mult = 1)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting19.B1 <- function(n) {
  cov <- gen_zu(n, 0.4, 0, 2, 0.6, 3, 2)
  mu_x <- 0.2 + 0.6*cov$z1 + 0.6*cov$z2 + 0.3*cov$u1 + 0.3*cov$u2
  y <- gen_y_binary(mu_x)
  ps <- ps_complement(0.2100, 0.2, -0.4, -0.4, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, cov$u1 + cov$u2 + cov$z1 + cov$z2, mult = 1)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

setting19.B2 <- function(n) {
  cov <- gen_zu(n, 0.4, 0, 2, 0.6, 3, 2)
  mu_x <- 0.2 + 0.6*cov$z1 + 0.6*cov$z2 + 0.3*cov$u1 + 0.3*cov$u2
  y <- gen_y_binary(mu_x)
  ps <- ps_complement(-0.6400, 0.2, -0.4, -0.4, y, cov$u1, cov$u2)
  ps <- apply_tilt(ps, n, cov$u1 + cov$u2 + cov$z1 + cov$z2, mult = 1)
  r <- sample_r(ps)
  df_zu(cov, y, r)
}

#==============================================================================
# Section 8. Cho settings (standard-logistic PS; covariates x1,x2,x3)
#==============================================================================

# --- Cho_M1 / Cho_M2 base (mean=1, sd=sqrt(1/3)) ---
Cho_M1.A1 <- function(n, response.rate) {
  cov <- gen_x123(n, mean = 1, sd = sqrt(1/3))
  mu_x <- 0.5 + cov$x1 + 0.5 * cov$x2
  y <- gen_y_normal(mu_x, sd = sqrt(1/3))
  ps <- ps_standard(-0.98 + 0.5*cov$x1 + 0.25*y)
  r <- sample_r(ps)
  data.frame(x1 = cov$x1, x2 = cov$x2, x12 = cov$x1^2, x22 = cov$x2^2, y = y, r = r)
}

Cho_M1.A2 <- function(n, response.rate) {
  cov <- gen_x123(n, mean = 1, sd = sqrt(1/3))
  # Original additionally draws an x4 that is never used. Honour that:
  x4 <- rnorm(n, mean = 1, sd = sqrt(1/3))
  mu_x <- 0.5 + cov$x1 + 0.5 * cov$x2
  y <- gen_y_normal(mu_x, sd = sqrt(1/3))
  ps <- ps_standard(-0.114 + 0.5*cov$x1 + 0.25*y)
  r <- sample_r(ps)
  data.frame(x1 = cov$x1, x2 = cov$x2, x3 = cov$x1^2, x4 = cov$x2^2, y = y, r = r)
}

Cho_M2.A1 <- function(n, response.rate) {
  x1 <- rnorm(n, mean = 1, sd = sqrt(1/3))
  x2 <- rnorm(n, mean = 1, sd = sqrt(1/3))
  mu_x <- 0.5 + x1 + 0.5 * x2
  y <- gen_y_normal(mu_x, sd = sqrt(1/3))
  ps <- ps_standard(0.02 + 0.5*x2 - 0.25*y)
  r <- sample_r(ps)
  data.frame(x1 = x1, x2 = x2, y = y, r = r)
}

Cho_M2.A2 <- function(n, response.rate) {
  x1 <- rnorm(n, mean = 1, sd = sqrt(1/3))
  x2 <- rnorm(n, mean = 1, sd = sqrt(1/3))
  mu_x <- 0.5 + x1 + 0.5 * x2
  y <- gen_y_normal(mu_x, sd = sqrt(1/3))
  ps <- ps_standard(0.865 + 0.5*x2 - 0.25*y)
  r <- sample_r(ps)
  data.frame(x1 = x1, x2 = x2, y = y, r = r)
}

# --- Cho_RM4 (mean of two PS branches, split on y <= 2) ---
Cho_RM4 <- function(n, response.rate) {
  x1 <- rnorm(n, mean = 1, sd = sqrt(1/3))
  x2 <- rnorm(n, mean = 1, sd = sqrt(1/3))
  mu_x <- 0.5 + x1 + 0.5 * x2
  y <- gen_y_normal(mu_x, sd = sqrt(1/3))
  eta <- ifelse(y <= 2, 1, 0) * (0.857 + 0.5*x1 - 0.25*y) +
         ifelse(y >  2, 1, 0) * (0.865 + 0.5*x2 - 0.25*y)
  ps <- ps_standard(eta)
  r <- sample_r(ps)
  data.frame(x1 = x1, x2 = x2, y = y, r = r)
}

# --- Cho_M{1,2}_gamma005 (mean=0, sd_x3=1; nonlinear correction 0.05*(x1^2-1/3)) -
Cho_M1_gamma005.A1 <- function(n, response.rate) {
  cov <- gen_x123(n, mean = 0, sd = sqrt(1/3), sd_x3 = 1)
  mu_x <- 0.5 + 0.5*cov$x1 + cov$x2 + 0.8*cov$x3 + 0.05*(cov$x1^2 - 1/3)
  y <- gen_y_normal(mu_x, sd = sqrt(1/3))
  ps <- ps_standard(0.1244 - 0.5*cov$x1 - 0.25*y)
  r <- sample_r(ps)
  data.frame(x1 = cov$x1, x2 = cov$x2, x3 = cov$x3, y = y, r = r)
}

Cho_M1_gamma005.A2 <- function(n, response.rate) {
  cov <- gen_x123(n, mean = 0, sd = sqrt(1/3), sd_x3 = 1)
  mu_x <- 0.5 + 0.5*cov$x1 + cov$x2 + 0.8*cov$x3 + 0.05*(cov$x1^2 - 1/3)
  y <- gen_y_normal(mu_x, sd = sqrt(1/3))
  ps <- ps_standard(1.0098 - 0.5*cov$x1 - 0.25*y)
  r <- sample_r(ps)
  data.frame(x1 = cov$x1, x2 = cov$x2, x3 = cov$x3, x4 = cov$x2^2, y = y, r = r)
}

Cho_M2_gamma005.A1 <- function(n, response.rate) {
  cov <- gen_x123(n, mean = 0, sd = sqrt(1/3), sd_x3 = 1)
  mu_x <- 0.5 + 0.5*cov$x1 + cov$x2 + 0.8*cov$x3 + 0.05*(cov$x1^2 - 1/3)
  y <- gen_y_normal(mu_x, sd = sqrt(1/3))
  ps <- ps_standard(-0.125 + 0.5*cov$x2 + 0.25*y)
  r <- sample_r(ps)
  data.frame(x1 = cov$x1, x2 = cov$x2, x3 = cov$x3, y = y, r = r)
}

Cho_M2_gamma005.A2 <- function(n, response.rate) {
  cov <- gen_x123(n, mean = 0, sd = sqrt(1/3), sd_x3 = 1)
  mu_x <- 0.5 + 0.5*cov$x1 + cov$x2 + 0.8*cov$x3 + 0.05*(cov$x1^2 - 1/3)
  y <- gen_y_normal(mu_x, sd = sqrt(1/3))
  ps <- ps_standard(0.7679 + 0.5*cov$x2 + 0.25*y)
  r <- sample_r(ps)
  data.frame(x1 = cov$x1, x2 = cov$x2, x3 = cov$x3, y = y, r = r)
}

# --- Cho_M{1,2}_gamma000: like gamma005 but no nonlinear x1^2-1/3 term ---
Cho_M1_gamma000.A1 <- function(n, response.rate) {
  cov <- gen_x123(n, mean = 0, sd = sqrt(1/3), sd_x3 = 1)
  mu_x <- 0.5 + 0.5*cov$x1 + cov$x2 + 0.8*cov$x3
  y <- gen_y_normal(mu_x, sd = sqrt(1/3))
  ps <- ps_standard(0.1246 - 0.5*cov$x1 - 0.25*y)
  r <- sample_r(ps)
  data.frame(x1 = cov$x1, x2 = cov$x2, x3 = cov$x3, y = y, r = r)
}

Cho_M1_gamma000.A2 <- function(n, response.rate) {
  # Original Cho_M1_gamma000.A2 returns 5 cols (x1,x2,x3,y,r) -- NOT x4 like the
  # gamma005 sibling. Honour that asymmetry.
  cov <- gen_x123(n, mean = 0, sd = sqrt(1/3), sd_x3 = 1)
  mu_x <- 0.5 + 0.5*cov$x1 + cov$x2 + 0.8*cov$x3
  y <- gen_y_normal(mu_x, sd = sqrt(1/3))
  ps <- ps_standard(1.0099 - 0.5*cov$x1 - 0.25*y)
  r <- sample_r(ps)
  data.frame(x1 = cov$x1, x2 = cov$x2, x3 = cov$x3, y = y, r = r)
}

Cho_M2_gamma000.A1 <- function(n, response.rate) {
  cov <- gen_x123(n, mean = 0, sd = sqrt(1/3), sd_x3 = 1)
  mu_x <- 0.5 + 0.5*cov$x1 + cov$x2 + 0.8*cov$x3
  y <- gen_y_normal(mu_x, sd = sqrt(1/3))
  ps <- ps_standard(-0.125 + 0.5*cov$x2 + 0.25*y)
  r <- sample_r(ps)
  data.frame(x1 = cov$x1, x2 = cov$x2, x3 = cov$x3, y = y, r = r)
}

Cho_M2_gamma000.A2 <- function(n, response.rate) {
  cov <- gen_x123(n, mean = 0, sd = sqrt(1/3), sd_x3 = 1)
  mu_x <- 0.5 + 0.5*cov$x1 + cov$x2 + 0.8*cov$x3
  y <- gen_y_normal(mu_x, sd = sqrt(1/3))
  ps <- ps_standard(0.7679 + 0.5*cov$x2 + 0.25*y)
  r <- sample_r(ps)
  data.frame(x1 = cov$x1, x2 = cov$x2, x3 = cov$x3, y = y, r = r)
}
