
## identical calls (but different random values)
(rt1 <- rdiffusion(20, a=1, v=2, t0=0.5))
(rt2 <- rdiffusion(20, a=1, v=2, t0=0.5, z=0.5, d=0, sz=0, sv=0, st0=0))

# get density for random RTs:
ddiffusion(rt1$rt, rt1$response, a=1, v=2, t0=0.5)  # boundary is factor
ddiffusion(rt1$rt, as.numeric(rt1$response), a=1, v=2, t0=0.5) # boundary is numeric
ddiffusion(rt1$rt, as.character(rt1$response), a=1, v=2, t0=0.5) # boundary is character

ddiffusion(rt2$rt, rt2$response, a=1, v=2, t0=0.5)


# plot density:
curve(ddiffusion(x, a=1, v=2, t0=0.5, boundary = "upper"), 
      xlim=c(0,3), main="Density of upper responses", ylab="density", xlab="response time")
curve(ddiffusion(x, a=1, v=2, t0=0.5, st0=0.2, boundary = "upper"), 
      add=TRUE, lty = 2)
legend("topright", legend=c("no", "yes"), title = "Starting Point Variability?", lty = 1:2)

# plot cdf:
curve(pdiffusion(x, a=1, v=2, t0=0.5, st0=0.2, boundary="u"), 
     xlim = c(0, 3),ylim = c(0,1), 
     ylab = "cumulative probability", xlab = "response time",
     main = "CDF of diffusion model with start point variability")
curve(pdiffusion(x, a=1, v=2, t0=0.5, st0=0.2, boundary="l"), 
     add=TRUE, lty = 2)
legend("topleft", legend=c("upper", "lower"), title="boundary", lty=1:2)


### qLBA can only return values up to maximal predicted probability:
# maximum probability for a given set
pdiffusion(20, a=1, v=2, t0=0.5, st0=0.2, sz = 0.1, sv = 0.5, boundary="u")
# [1] 0.8705141
# pdiffusion(Inf, a=1, v=2, t0=0.5, st0=0.2, sz = 0.1, sv = 0.5, boundary="u") # equal but much slower

qdiffusion(0.87, a=1, v=2, t0=0.5, st0=0.2, sz = 0.1, sv = 0.5, boundary="u")
# [1] 1.769253

qdiffusion(0.88, a=1, v=2, t0=0.5, st0=0.2, sz = 0.1, sv = 0.5, boundary="u")
# NA
