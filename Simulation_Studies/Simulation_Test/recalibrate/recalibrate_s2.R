## Recalibrate alpha0 for setting2 (binary Y)
## Propensity = 1/(1+exp(alpha0 + 0.8*y - 0.5*u1 + 0.5*u2))
## m(x) = 1+u1+u2+z1+z2; z1~Bern(0.4), z2~N(0,1), u1~Bern(0.6), u2~N(0,1)
## B tilt: exp(10*n^(-1/2)*(y+u1+u2)), clip 0.95
setwd("c:/Users/stat-user/iCloudDrive/Desktop/EBMR/Simulation_Studies")

set.seed(12345)
N <- 1e7
z1 <- rbinom(N, 1, 0.4); z2 <- rnorm(N, 0, 1)
u1 <- rbinom(N, 1, 0.6); u2 <- rnorm(N, 0, 1)
mx <- 1 + u1 + u2 + z1 + z2
y <- rbinom(N, 1, exp(mx)/(1+exp(mx)))
cat(sprintf("setting2 mu_true (mean y) = %.4f\n\n", mean(y)))

lin <- function(a0) a0 + 0.5*y - 0.5*u1 + 0.5*u2
sA <- function(target) uniroot(function(a0) mean(1/(1+exp(lin(a0)))) - target, c(-15,15))$root
sB <- function(target, nt) {
  f <- function(a0) {
    p <- 1/(1+exp(lin(a0))) * exp(10*nt^(-1/2)*(y+u1+u2))
    p[p > 1] <- 0.95
    mean(p) - target
  }
  uniroot(f, c(-15,15))$root
}

A1 <- sA(0.5); A2 <- sA(0.7)
B1_2000 <- sB(0.5, 2000); B1_500 <- sB(0.5, 500)
B2_2000 <- sB(0.7, 2000); B2_500 <- sB(0.7, 500)

cat(sprintf("A1: %.4f   A2: %.4f\n", A1, A2))
cat(sprintf("B1 n>=1000: %.4f   B1 n<1000: %.4f\n", B1_2000, B1_500))
cat(sprintf("B2 n>=1000: %.4f   B2 n<1000: %.4f\n", B2_2000, B2_500))

cat("\n=== Verification ===\n")
v <- function(a0, nt=NULL) {
  p <- 1/(1+exp(lin(a0)))
  if (!is.null(nt)) { p <- p*exp(10*nt^(-1/2)*(y+u1+u2)); p[p>1] <- 0.95 }
  mean(p)
}
cat(sprintf("A1=%.4f A2=%.4f B1(2000)=%.4f B1(500)=%.4f B2(2000)=%.4f B2(500)=%.4f\n",
            v(A1), v(A2), v(B1_2000,2000), v(B1_500,500), v(B2_2000,2000), v(B2_500,500)))
