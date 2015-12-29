
## random number generation using different distributions for v:
rlba_norm(10, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=c(0.2,0.3))
rlba_gamma(10, A=0.5, b=1, t0 = 0.5, shape_v=c(1.2, 1), scale_v=c(0.2,0.3))
rlba_frechet(10, A=0.5, b=1, t0 = 0.5, shape_v=c(1.2, 1), scale_v=c(0.2,0.3))
rlba_lnorm(10, A=0.5, b=1, t0 = 0.5, meanlog_v=c(1.2, 1), sdlog_v=c(0.2, 0.3))

# use somewhat plausible values for plotting:
A <- 0.2
b <- 0.5
t0 <- 0.3

# plot density:
curve(dlba_norm(x, A=A, b=b, t0=t0, mean_v = 1.0, sd_v = 0.5), ylim = c(0, 4),
      xlim=c(0,3), main="Density/PDF of LBA versions", ylab="density", xlab="response time")
curve(dlba_gamma(x, A=A, b=b, t0=t0, shape_v=1, scale_v=1), add=TRUE, lty = 2)
curve(dlba_frechet(x, A=A, b=b, t0=t0, shape_v=1,scale_v=1.0), add=TRUE, lty = 3)
curve(dlba_lnorm(x, A=A, b=b, t0=t0, meanlog_v = 0.5, sdlog_v = 0.5), add=TRUE, lty = 4)
legend("topright", legend=c("Normal", "Gamma", "Frechet", "Log-Normal"), 
      title = expression("Distribution of"~~italic(v)), lty = 1:4)


# plot cdf:
curve(plba_norm(x, A=A, b=b, t0=t0, mean_v=1.0, sd_v=1.0), 
      xlim = c(0, 3),ylim = c(0,1), 
      ylab = "cumulative probability", xlab = "response time",
      main = "Distribution/CDF of LBA versions")
curve(plba_gamma(x, A=A, b=b, t0=t0, shape_v=1,scale_v=1), add=TRUE, lty = 2)
curve(plba_frechet(x, A=A, b=b, t0=t0, shape=1, scale=1), add=TRUE, lty = 3)
curve(plba_lnorm(x, A=A, b=b, t0=t0, meanlog_v=0.5, sdlog_v = 0.5), add=TRUE, lty = 4)
legend("bottomright", legend=c("Normal", "Gamma", "Frechet", "Log-Normal"), 
       title = expression("Distribution of"~~italic(v)), lty = 1:4)
