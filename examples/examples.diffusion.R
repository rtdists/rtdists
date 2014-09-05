
## identical calls (but different random values)
rrd(10, a=1, v=2, t0=0.5)
rrd(10, a=1, v=2, t0=0.5, z=0.5, d=0, sz=0, sv=0, st0=0)


# plot density:
curve(drd(x, a=1, v=2, t0=0.5, boundary = "upper"), 
      xlim=c(0,3), main="Density of upper responses", ylab="density", xlab="response time")
curve(drd(x, a=1, v=2, t0=0.5, st0=0.2, boundary = "upper"), 
      add=TRUE, lty = 2)
legend("topright", legend=c("no", "yes"), title = "Starting Point Variability?", lty = 1:2)

# plot cdf:
curve(prd(x, a=1, v=2, t0=0.5, st0=0.2, boundary="u"), 
     xlim = c(0, 3),ylim = c(0,1), 
     ylab = "cumulative probability", xlab = "response time",
     main = "CDF of diffusion model with start point variability")
curve(prd(x, a=1, v=2, t0=0.5, st0=0.2, boundary="l"), 
     add=TRUE, lty = 2)
legend("topleft", legend=c("upper", "lower"), title="boundary", lty=1:2)

