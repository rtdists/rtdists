rm (list=ls())
# UNLOAD EXISTING RTDISTS, LOAD IN OLD VERSION FOR COMPARISON TESTING
devtools::unload(rtdists)
#install.packages ("rtdists", lib="tests\\temp_testing\\old_rtdists_0.6-6\\")
require (rtdists, lib.loc="tests\\temp_testing\\old_rtdists_0.6-6")

# Load in test bindings for new RCppFastDM PDF version
Rcpp::sourceCpp('src\\RFastDM.cpp', rebuild=TRUE)
source ("tests\\temp_testing\\RCppDensity.R")


# SCALAR TESTING - COMPARE CURVES
a   <- 2      
zr  <- 0.7
v   <- 1.3
t0  <- 0.2
d   <- 0
szr <- 0.3
sv  <- 0.4
st0 <- 0.5

bound <- "upper"
precision <- 3

curve( ddiffusion       (x, response=bound, a=a, v=v, t0=t0, z=zr, d=d, sz=szr, sv=sv, st0=st0, precision=precision), n=1001, add = FALSE, lty = 2, col="red",    xlim=c(0,5), ylim=c(0,1))
curve( rcpp_ddiffusion  (x, response=bound, a=a, v=v, t0=t0, z=zr,   d=d, sz=szr, sv=sv, st0=st0, precision=precision), n=1001, add = TRUE, lty = 2, col="blue",    xlim=c(0,5), ylim=c(0,1))


# VECTORISED TESTING - COMPARE CURVES
a   <- seq(1.0, 2.0, by=0.1)
zr  <- seq(0.3, 0.7, by=0.05)      
v   <- seq(1.0, 2.0, by=0.05)
t0  <- seq(0.1, 0.3, by=0.05)
d   <- seq(0.0, 0.1, by=0.05)
szr <- seq(0.3, 0.4, by=0.05)
sv  <- seq(0.3, 0.5, by=0.05)
st0 <- seq(0.1, 0.2, by=0.05)

curve( ddiffusion       (x, response=bound, a=a, v=v, t0=t0, z=zr, d=d, sz=szr, sv=sv, st0=st0, precision=precision), n=1001, add = FALSE, lty = 2, col="red",    xlim=c(0,5), ylim=c(0,1))
curve( rcpp_ddiffusion  (x, response=bound, a=a, v=v, t0=t0, z=zr,   d=d, sz=szr, sv=sv, st0=st0, precision=precision), n=1001, add = TRUE, lty = 2, col="blue",    xlim=c(0,5), ylim=c(0,1))


# SPEED CHECK
rts <- seq(0,5,length.out=1e6)

system.time({ 
  orig_ddiff <- ddiffusion (rts, response=bound, a=a, v=v, t0=t0, z=zr*a, d=d, sz=szr, sv=sv, st0=st0, precision=precision) 
})
system.time({ 
  rcpp_ddiff <- rcpp_ddiffusion (rts, response=bound, a=a, v=v, t0=t0, z=zr, d=d, sz=szr, sv=sv, st0=st0, precision=precision) 
})
