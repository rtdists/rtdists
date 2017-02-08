# UNLOAD EXISTING RTDISTS, LOAD IN OLD VERSION FOR COMPARISON TESTING
devtools::unload(rtdists)
install.packages ("rtdists", lib="tests\\temp_testing\\old_rtdists_0.6-6\\")
require (rtdists, lib.loc="tests\\temp_testing\\old_rtdists_0.6-6")

# Load in test bindings for new RCppFastDM CDF versions
Rcpp::sourceCpp('src\\RFastDM.cpp', rebuild=TRUE)
source ("tests\\temp_testing\\RCppDistribution.R")


# SCALAR TESTING - COMPARE CURVES
# pdiffusion is integrated-PDF  (exactly equal to pwiener if vars=0; note use of scaled z)
# rcpp_fast_diffusion is 'fast' fast-dm implementation
# rcpp_pdiffusion is 'slow' fast-dm plot-cdf implementation
#
# ACCURACY NOTES:
#   UPPER BOUND:
#     1. all curves exactly overlay when zr = 0.5 (whether variances = 0 or not)
#     2. rcpp_diffusion exactly equals precise_rcpp_diffusion when variances = 0 and zr != 0.5
#     3. pdiffusion exactly equals precise_rcpp_diffusion (whether variances = 0 or not)
#     4. rdiffusion(EXISTING) ecdf exactly equals pdiffusion/precise_rcpp_diffusion (whether variances = 0 or not) 
#   LOWER BOUND:
#     1. all curves exactly overlay when zr = 0.5 (whether variances = 0 or not)
#     2. pdiffusion exactly equals precise_rcpp_diffusion (whether variances = 0 or not and whether zr = 0.5 or not)
#     3  rcpp_pdiffusion is way off for some parameters
#     4. ?rdiffusion ecdf is slightly lower when zr=0.5 and variances = 0 (==> likely sampling error)
#      

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

curve( pdiffusion              (x, response=bound, a=a, v=v, t0=t0, z=zr*a, d=d, sz=szr, sv=sv, st0=st0, precision=precision), n=1001, add = FALSE, lty = 2, col="red",    xlim=c(0,5), ylim=c(0,1))
curve( rcpp_fast_pdiffusion    (x, response=bound, a=a, v=v, t0=t0, z=zr,       d=d, sz=szr, sv=sv, st0=st0, precision=precision), n=1001, add = TRUE, lty = 2, col="green",    xlim=c(0,5), ylim=c(0,1))
curve( rcpp_pdiffusion         (x, response=bound, a=a, v=v, t0=t0, z=zr,       d=d, sz=szr, sv=sv, st0=st0, precision=precision), n=1001, add = TRUE, lty = 2, col="blue",    xlim=c(0,5), ylim=c(0,1))

# Note used of scaled z; EXISTING rdiffusion (though should be exactly equal)
g1 <- rdiffusion(n = 1e6, a=a, v=v, t0=t0, z=zr*a, d=d, sz=szr, sv=sv, st0=st0, precision=precision)
g1_prop <- mean(g1$response ==bound)
emp_cdf_rdiffusion <- ecdf(g1[g1$response == bound, "rt"])
curve(emp_cdf_rdiffusion(x)*g1_prop, 0, 3, ylim = c(0, 1), add = TRUE, col="orange", lty = 4)

# VECTORISED TESTING - COMPARE CURVES
a   <- seq(1.0, 2.0, by=0.1)  
zr  <- seq(0.3, 0.7, by=0.05)      
# v   <- seq(1.0, 2.0, by=0.05)
# t0  <- seq(0.1, 0.3, by=0.05)
# d   <- seq(0.0, 0.1, by=0.05)
# szr <- seq(0.3, 0.4, by=0.05)
# sv  <- seq(0.3, 0.5, by=0.05)
# st0 <- seq(0.1, 0.2, by=0.05)

bound <- "upper"
precision <- 3

curve( pdiffusion              (x, response=bound, a=a, v=v, t0=t0, z=zr*a, d=d, sz=szr, sv=sv, st0=st0, precision=precision), n=1001, add = FALSE, lty = 2, col="red",    xlim=c(0,5), ylim=c(0,1))
curve( rcpp_fast_pdiffusion    (x, response=bound, a=a, v=v, t0=t0, z=zr,       d=d, sz=szr, sv=sv, st0=st0, precision=precision), n=1001, add = TRUE, lty = 2, col="green",    xlim=c(0,5), ylim=c(0,1))
curve( rcpp_pdiffusion         (x, response=bound, a=a, v=v, t0=t0, z=zr,       d=d, sz=szr, sv=sv, st0=st0, precision=precision), n=1001, add = TRUE, lty = 2, col="blue",    xlim=c(0,5), ylim=c(0,1))

# Note used of scaled z; EXISTING rdiffusion (though should be exactly equal)
g1 <- rdiffusion(n = 1e6, a=a, v=v, t0=t0, z=zr*a, d=d, sz=szr, sv=sv, st0=st0, precision=precision)
g1_prop <- mean(g1$response ==bound)
emp_cdf_rdiffusion <- ecdf(g1[g1$response == bound, "rt"])
curve(emp_cdf_rdiffusion(x)*g1_prop, 0, 3, ylim = c(0, 1), add = TRUE, col="orange", lty = 4)


# CHECK SPEED
# With no variation (1e5): 
#   - (pdiffusion) integrated PDF:							                151.23
#   - (rcpp_pdiffusion) rcpp pdiffusion using F_get_F/F_get_N:	  1.23
#   - (precise_rcpp_diffusion) rcpp pdiffusion using F_get_val:	  1.22
# 
# With variance in z, v, t0 (1e5): 
# On 6700K
#   - integrated PDF:							            513.47		/ 348.9   / 852.3
#   - rcpp pdiffussion using F_get_F/F_get_N:	  1.39		/   1.40  /   1.39
#   - rcpp pdiffusion using F_get_val:			    1.34		/   1.31  /   1.23
#
# On Mac OSX ('2.2Ghz i7')
#   - integrated PDF:							            611.139   / 609.881
#   - rcpp pdiffussion using F_get_F/F_get_N:	  1.678   /   1.675 
#   - rcpp pdiffusion using F_get_val:			    1.739   /   1.790 
# 
# On 5930K
#   - integrated PDF:							            1101.4    / 1079.61
#   - rcpp pdiffussion using F_get_F/F_get_N:	  2.39    /    2.24  /   2.10  / 2.48 / 3.16
#   - rcpp pdiffusion using F_get_val:			    2.52    /    2.16  /   2.32  / 1.95 / 2.35

rts <- seq(0,5,length.out=1e5)

system.time({ 
  orig_pdiff <- pdiffusion (rts, response=bound, a=a, v=v, t0=t0, z=zr*a, d=d, sz=szr, sv=sv, st0=st0, precision=precision) 
})
system.time({ 
  rcpp_pdiff <- rcpp_fast_pdiffusion (rts, response=bound, a=a, v=v, t0=t0, z=zr, d=d, sz=szr, sv=sv, st0=st0, precision=precision) 
})
system.time({ 
  rcpp_pdiff_fast <- rcpp_pdiffusion (rts, response=bound, a=a, v=v, t0=t0, z=zr, d=d, sz=szr, sv=sv, st0=st0, precision=precision) 
})
