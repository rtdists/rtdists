
## random number generation using different distributions for v:
rlba_norm(10, A=0.5, b=1, v=c(1.2, 1), sv=c(1,1.5), t0 = 0.5)
rlba_gamma(10, A=0.5, b=1, v=c(1.2, 1), sv=c(1, 1.5), t0 = 0.5)
rlba_frechet(10, A=0.5, b=1, v=c(1.2, 1), sv=c(1, 1.5), t0 = 0.5)
rlba_lnorm(10, A=0.5, b=1, v=c(1.2, 1), sv=c(1, 1.5), t0 = 0.5)

# use somewhat plausible values for plotting:
A <- 0.2
b <- 0.5
v <- 1.0
t0 <- 0.3
sv <- 1.1  # note that frechet needs sv > 1.

# plot density:
curve(dlba_norm(x, A=A, b=b, v=v, t0=t0, sv=sv), ylim = c(0, 5),
      xlim=c(0,3), main="Density of LBA versions", ylab="density", xlab="response time")
curve(dlba_gamma(x, A=A, b=b, v=v, t0=t0, sv=sv), 
      add=TRUE, lty = 2)
curve(dlba_frechet(x, A=A, b=b, v=v, t0=t0, sv=sv), 
      add=TRUE, lty = 3)
curve(dlba_lnorm(x, A=A, b=b, v=v, t0=t0, sv=sv), 
      add=TRUE, lty = 4)
legend("topright", legend=c("Normal", "Gamma", "Frechet", "Log-Normal"), 
      title = expression("Distribution of"~~italic(v)), lty = 1:4)


# plot cdf:
curve(plba_norm(x, A=A, b=b, v=v, t0=t0, sv=sv), 
      xlim = c(0, 3),ylim = c(0,1), 
      ylab = "cumulative probability", xlab = "response time",
      main = "CDF of LBA versions")
curve(plba_gamma(x, A=A, b=b, v=v, t0=t0, sv=sv), 
      add=TRUE, lty = 2)
curve(plba_frechet(x, A=A, b=b, v=v, t0=t0, sv=sv), 
      add=TRUE, lty = 3)
curve(plba_lnorm(x, A=A, b=b, v=v, t0=t0, sv=sv), 
      add=TRUE, lty = 4)
legend("bottomright", legend=c("Normal", "Gamma", "Frechet", "Log-Normal"), 
       title = expression("Distribution of"~~italic(v)), lty = 1:4)
