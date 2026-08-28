## random number generation for a single Wald accumulator:
rwald(10, A=0.5, b=1, t0 = 0.5, v=c(1.2, 1))

# use somewhat plausible values for plotting:
A <- 0.2
b <- 0.5
t0 <- 0.3

# plot density:
curve(dwald(x, A=A, b=b, t0=t0, v=2), ylim = c(0, 4), xlim = c(0, 3),
      main="Density/PDF of the Wald accumulator", ylab="density",
      xlab="response time")
curve(dwald(x, A=0, b=b, t0=t0, v=2), add=TRUE, lty = 2)
legend("topright", legend=c("0.2", "0"),
       title = expression("Start point variability"~~italic(A)), lty = 1:2)

# plot cdf:
curve(pwald(x, A=A, b=b, t0=t0, v=2), xlim = c(0, 3), ylim = c(0, 1),
      ylab = "cumulative probability", xlab = "response time",
      main = "Distribution/CDF of the Wald accumulator")
curve(pwald(x, A=0, b=b, t0=t0, v=2), add=TRUE, lty = 2)
legend("bottomright", legend=c("0.2", "0"),
       title = expression("Start point variability"~~italic(A)), lty = 1:2)
