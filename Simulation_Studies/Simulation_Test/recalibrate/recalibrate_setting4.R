## Recalibrate alpha0 for setting4
## Propensity = 1/(1+exp(alpha0 + y - 0.5*u1 + 0.5*u2))
## m(x) = 1 + u1 + u2 + z1 + z2; z1~Bern(0.4), z2~N(0,1), u1~Bern(0.6), u2~N(0,1)
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

set.seed(12345)
N <- 1e7

z1 <- rbinom(N, size = 1, prob = 0.4)
z2 <- rnorm(N, mean = 0, sd = 1)
u1 <- rbinom(N, size = 1, prob = 0.6)
u2 <- rnorm(N, mean = 0, sd = 1)
eta <- 1 + u1 + u2 + z1 + z2
y <- rbinom(N, size = 1, prob = exp(eta)/(1+exp(eta)))

cat(sprintf("mu_true (mean y) = %.4f\n\n", mean(y)))

lin <- function(a0) a0 + y - 0.5*u1 + 0.5*u2

solve_alpha0_A <- function(target) {
  uniroot(function(a0) mean(1/(1+exp(lin(a0)))) - target, c(-10, 10))$root
}
solve_alpha0_B <- function(target, n_target) {
  f <- function(a0) {
    p <- 1/(1+exp(lin(a0))) * exp(10*n_target^(-1/2)*(y+u1+u2))
    p[p > 1] <- 0.95
    mean(p) - target
  }
  uniroot(f, c(-10, 10))$root
}

a_A1 <- solve_alpha0_A(0.5)
a_A2 <- solve_alpha0_A(0.7)
a_B1_2000 <- solve_alpha0_B(0.5, 2000)
a_B1_500  <- solve_alpha0_B(0.5, 500)
a_B2_2000 <- solve_alpha0_B(0.7, 2000)
a_B2_500  <- solve_alpha0_B(0.7, 500)

cat(sprintf("setting4.A1: alpha0 = %.4f\n", a_A1))
cat(sprintf("setting4.A2: alpha0 = %.4f\n", a_A2))
cat(sprintf("setting4.B1 n>=1000: alpha0 = %.4f\n", a_B1_2000))
cat(sprintf("setting4.B1 n<1000:  alpha0 = %.4f\n", a_B1_500))
cat(sprintf("setting4.B2 n>=1000: alpha0 = %.4f\n", a_B2_2000))
cat(sprintf("setting4.B2 n<1000:  alpha0 = %.4f\n", a_B2_500))

cat("\n=== Verification ===\n")
verify <- function(a0, tilt_n = NULL) {
  p <- 1/(1+exp(lin(a0)))
  if (!is.null(tilt_n)) { p <- p * exp(10*tilt_n^(-1/2)*(y+u1+u2)); p[p > 1] <- 0.95 }
  mean(p)
}
cat(sprintf("A1: %.4f (target 0.5)\n", verify(a_A1)))
cat(sprintf("A2: %.4f (target 0.7)\n", verify(a_A2)))
cat(sprintf("B1 n=2000: %.4f (target 0.5)\n", verify(a_B1_2000, 2000)))
cat(sprintf("B1 n=500:  %.4f (target 0.5)\n", verify(a_B1_500, 500)))
cat(sprintf("B2 n=2000: %.4f (target 0.7)\n", verify(a_B2_2000, 2000)))
cat(sprintf("B2 n=500:  %.4f (target 0.7)\n", verify(a_B2_500, 500)))
