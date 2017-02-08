rm (list=ls())
# UNLOAD EXISTING RTDISTS, LOAD IN OLD VERSION FOR COMPARISON TESTING
devtools::unload(rtdists)
#install.packages ("rtdists", lib="tests\\temp_testing\\old_rtdists_0.6-6\\")
require (rtdists, lib.loc="tests\\temp_testing\\old_rtdists_0.6-6")

# Load in test bindings for new RCppFastDM Sampling version
Rcpp::sourceCpp('src\\RFastDM.cpp', rebuild=TRUE)
source ("tests\\temp_testing\\RCppSampling.R")

a   <- 2      
zr  <- 0.5
v   <- 1.5
t0  <- 0.2
d   <- 0
szr <- 0.2 # 0.5
sv  <- 0.4 # 0.5
st0 <- 0.6 # 0.5

precision <- 3

n <- 1e5

system.time ({
  orig <- rdiffusion      (n=n, a=a, v=v, t0=t0, z=zr, d=d, sz=szr, sv=sv, st0=st0, precision=precision)
})

system.time ({
  rcpp <- rcpp_rdiffusion (n=n, a=a, v=v, t0=t0, z=zr, d=d, sz=szr, sv=sv, st0=st0, precision=precision)
})

orig_upper <- orig[orig$response == "upper",]
rcpp_upper <- rcpp[rcpp$response == "upper",]

orig_lower <- orig[orig$response == "lower",]
rcpp_lower <- rcpp[rcpp$response == "lower",]

plot (density (orig_upper$rt), col="blue")
lines(density (rcpp_upper$rt), col="red")

lines(density (orig_lower$rt), col="green")
lines(density (rcpp_lower$rt), col="purple")

# VECTORISED TESTING - COMPARE CURVES
a   <- seq(1.0, 2.0, by=0.1)
zr  <- seq(0.3, 0.7, by=0.05)      
v   <- seq(1.0, 2.0, by=0.05)
t0  <- seq(0.1, 0.3, by=0.05)
d   <- seq(0.0, 0.1, by=0.05)
szr <- seq(0.3, 0.4, by=0.05)
sv  <- seq(0.3, 0.5, by=0.05)
st0 <- seq(0.1, 0.2, by=0.05)