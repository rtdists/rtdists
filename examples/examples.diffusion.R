
## identical calls (but different random values)
rt1 <- rdiffusion(500, a=1, v=2, t0=0.5)
head(rt1)
rt2 <- rdiffusion(500, a=1, v=2, t0=0.5, z=0.5, d=0, sz=0, sv=0, st0=0)
head(rt2)
  
# get density for random RTs (possible to specify arguments for pdiffusion in same way):
sum(log(ddiffusion(rt1$rt, rt1$response, a=1, v=2, t0=0.5)))  # response is factor
sum(log(ddiffusion(rt1$rt, as.numeric(rt1$response), a=1, v=2, t0=0.5))) # response is numeric
sum(log(ddiffusion(rt1$rt, as.character(rt1$response), a=1, v=2, t0=0.5))) # response is character
sum(log(ddiffusion(rt1, a=1, v=2, t0=0.5))) # response is data.frame


sum(log(ddiffusion(rt2$rt, rt2$response, a=1, v=2, t0=0.5)))

# can we recover the parameters?
ll_diffusion <- function(pars, rt, response) 
{
  densities <- tryCatch(
    ddiffusion(rt, response=response, 
               a=pars[1], v=pars[2], t0=pars[3], 
               sz=pars[4], 
               st0=pars[5], sv=pars[6]), 
    error = function(e) 0)
  if (any(densities == 0)) return(1e6)
  return(-sum(log(densities)))
}

\dontrun{
start <- c(runif(2, 0.5, 3), 0.1, runif(3, 0, 0.5))
names(start) <- c("a", "v", "t0", "sz", "st0", "sv")
recov <- nlminb(start, ll_diffusion, lower = 0, rt=rt1$rt, response=rt1$response)
round(recov$par, 3)
#    a     v    t0    sz   st0    sv 
#0.954 1.764 0.503 0.000 0.000 0.000 
}

# plot density:
curve(ddiffusion(x, a=1, v=2, t0=0.5, response = "upper"), 
      xlim=c(0,3), main="Density of upper responses", ylab="density", xlab="response time")
curve(ddiffusion(x, a=1, v=2, t0=0.5, st0=0.2, response = "upper"), 
      add=TRUE, lty = 2)
legend("topright", legend=c("no", "yes"), title = "Starting Point Variability?", lty = 1:2)

# plot cdf:
curve(pdiffusion(x, a=1, v=2, t0=0.5, st0=0.2, response="u"), 
     xlim = c(0, 3),ylim = c(0,1), 
     ylab = "cumulative probability", xlab = "response time",
     main = "CDF of diffusion model with start point variability")
curve(pdiffusion(x, a=1, v=2, t0=0.5, st0=0.2, response="l"), 
     add=TRUE, lty = 2)
legend("topleft", legend=c("upper", "lower"), title="response boundary", lty=1:2)

\dontrun{
### qLBA can only return values up to maximal predicted probability:
(max_p <- pdiffusion(Inf, a=1, v=2, t0=0.5, st0=0.2, sz = 0.1, sv = 0.5, response="u"))
# [1] 0.87
# (Note that with the current integration routine for pdiffusion use Inf and not smaller values.)

qdiffusion(0.87, a=1, v=2, t0=0.5, st0=0.2, sz = 0.1, sv = 0.5, response="u")
# [1] 1.945802

qdiffusion(0.88, a=1, v=2, t0=0.5, st0=0.2, sz = 0.1, sv = 0.5, response="u")
# NA with warning.

# to get predicted quantiles, scale required quantiles by maximally predicted response rate:
qs <- c(.1, .3, .5, .7, .9)
qdiffusion(qs*max_p, a=1, v=2, t0=0.5, st0=0.2, sz = 0.1, sv = 0.5, response="u")

# or set scale_p to TRUE which scales automatically by maximum p
# (but can be slow as it calculates max_p for each probability separately) 
qdiffusion(qs, a=1, v=2, t0=0.5, st0=0.2, sz = 0.1, sv = 0.5, response="u", scale_p = TRUE)


# qdiffusion also accepts a data.frame as first argument:
t3 <- data.frame(p = rep(c(0.05, 0.1, 0.87), 2), response = rep(c("upper", "lower"), each = 3))
#      p response
# 1 0.05    upper
# 2 0.10    upper
# 3 0.87    upper
# 4 0.05    lower
# 5 0.10    lower
# 6 0.87    lower
qdiffusion(t3, a=1, v=2, t0=0.5, st0=0.2, sz = 0.1, sv = 0.5)
}

## LBA and diffusion can be used interchangeably:
rt1 <- rLBA(500, A=0.5, b=1, t0 = 0.5, mean_v=c(2.4, 1.6), sd_v=c(1,1.2))
rt2 <- rdiffusion(500, a=1, v=2, t0=0.5)

# data can also be passed as data.frame (same is true for pLBA):
sum(log(dLBA(rt1, A=0.5, b=1, t0 = 0.5, mean_v=c(2.4, 1.6), sd_v=c(1,1.2))))
sum(log(dLBA(rt2, A=0.5, b=1, t0 = 0.5, mean_v=c(2.4, 1.6), sd_v=c(1,1.2))))

sum(log(ddiffusion(rt1, a=1, v=2, t0=0.5)))
sum(log(ddiffusion(rt2, a=1, v=2, t0=0.5)))


