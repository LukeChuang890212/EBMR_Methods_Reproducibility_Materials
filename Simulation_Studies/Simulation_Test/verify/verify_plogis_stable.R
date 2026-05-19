## Verify plogis is stable without clamping
x <- c(-100, -50, -20, 0, 20, 50, 100, 700, -700)
cat("plogis(-x) [logistic_complement]:\n")
print(data.frame(x = x, plogis_minus_x = plogis(-x), clamped = plogis(-pmin(pmax(x, -20), 20))))
cat("\nDifference within (-20, 20):", max(abs(plogis(-x[abs(x) <= 20]) - plogis(-pmin(pmax(x[abs(x) <= 20], -20), 20)))), "\n")

# For very extreme x, clamping changes the result
cat("\nDifference for |x|>20 (this is WHY clamping was added):\n")
big <- c(50, 100, -100)
cat("plogis(-50) =", plogis(-50), "vs plogis(-20) =", plogis(-20), "\n")
cat("If eta=-50, true pi=", plogis(50), ", clamped pi=", plogis(20), "\n")
