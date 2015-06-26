# Testing st0 implementation

library (devtools)

devtools::install_local ("rtdists-master (120615)/rtdists_0.4-1.tar.gz")
library (rtdists)

t1st0 <- rrd(1000, a=1, z=0.5, v=2, t0=1, d=0, sz = 0, sv = 0, st0 = 0)
t2st0 <- rrd(1000, a=1, z=0.5, v=2, t0=2, d=0, sz = 0, sv = 0, st0 = 0)

t1st1 <- rrd(1000, a=1, z=0.5, v=2, t0=1, d=0, sz = 0, sv = 0, st0 = 1)
t2st1 <- rrd(1000, a=1, z=0.5, v=2, t0=2, d=0, sz = 0, sv = 0, st0 = 1)

hist (t1st0$rt, breaks=100)
hist (t2st0$rt, breaks=100)
hist (t1st1$rt, breaks=100)
hist (t2st1$rt, breaks=100)


# t0 fix:
#   Henrik description for DDM:
#     t0: non-decision time or response time constant (in seconds). Average duration of all non-decisional processes
#     (encoding and response execution). Typical range: 0.1 < t0 < 0.5

library (devtools)

devtools::install_local ("rtdists-master (120615)/rtdists_0.3-1.tar.gz")
library (rtdists)


old_t0_001_st0_0 <- old_rrd(1000, a=1, z=0.5, v=2, t0= 1, d=0, sz = 0, sv = 0, st0 = 0)
old_t0_002_st0_0 <- old_rrd(1000, a=1, z=0.5, v=2, t0= 2, d=0, sz = 0, sv = 0, st0 = 0)
old_t0_001_st0_1 <- old_rrd(1000, a=1, z=0.5, v=2, t0= 1, d=0, sz = 0, sv = 0, st0 = 1)
old_t0_002_st0_1 <- old_rrd(1000, a=1, z=0.5, v=2, t0= 2, d=0, sz = 0, sv = 0, st0 = 1)

new_t0_001_st0_0 <- rrd(1000, a=1, z=0.5, v=2, t0= 1, d=0, sz = 0, sv = 0, st0 = 0)
new_t0_002_st0_0 <- rrd(1000, a=1, z=0.5, v=2, t0= 2, d=0, sz = 0, sv = 0, st0 = 0)
new_t0_001_st0_1 <- rrd(1000, a=1, z=0.5, v=2, t0= 1, d=0, sz = 0, sv = 0, st0 = 1)
new_t0_002_st0_1 <- rrd(1000, a=1, z=0.5, v=2, t0= 2, d=0, sz = 0, sv = 0, st0 = 1)


par(mfrow=c(2,2))
hist (old_t0_001_st0_0$rt, breaks=100, xlim=c(0,3), main="OLD DDM: t0 = 1, st0 = 0")
hist (old_t0_001_st0_1$rt, breaks=100, xlim=c(0,3), main="OLD DDM: t0 = 1, st0 = 1")

hist (new_t0_001_st0_0$rt, breaks=100, xlim=c(0,3), main="FIXED DDM: t0 = 1, st0 = 0")
hist (new_t0_001_st0_1$rt, breaks=100, xlim=c(0,3), main="FIXED DDM: t0 = 1, st0 = 1")


devtools::install_local ("rtdists-master (120615)/rtdists_0.4-2.tar.gz")
library (rtdists)
# 
# prd1_t0_1_st0_0 <- prd(1, a=1, z=0.5, v=2, t0= 1, d=0, sz = 0, sv = 0, st0 = 0)
# drd1_t0_1_st0_0 <- drd(1, a=1, z=0.5, v=2, t0= 1, d=0, sz = 0, sv = 0, st0 = 0)
# prd1_t0_1_st0_1 <- prd(1, a=1, z=0.5, v=2, t0= 1, d=0, sz = 0, sv = 0, st0 = 1)
# drd1_t0_1_st0_1 <- drd(1, a=1, z=0.5, v=2, t0= 1, d=0, sz = 0, sv = 0, st0 = 1)

prd0_t0_1_st0_0 <- prd(1.1, a=1, z=0.5, v=2, t0= 1, d=0, sz = 0, sv = 0, st0 = 0)
drd0_t0_1_st0_0 <- drd(1.1, a=1, z=0.5, v=2, t0= 1, d=0, sz = 0, sv = 0, st0 = 0)
prd0_t0_1_st0_1 <- prd(1.1, a=1, z=0.5, v=2, t0= 1, d=0, sz = 0, sv = 0, st0 = 1)
drd0_t0_1_st0_1 <- drd(1.1, a=1, z=0.5, v=2, t0= 1, d=0, sz = 0, sv = 0, st0 = 1)

rts <- seq (0, 10, by=0.01)
y_st0 <- prd(rts, a=1, z=0.5, v=2, t0= 1, d=0, sz = 0, sv = 0, st0 = 0)
y_st1 <- prd(rts, a=1, z=0.5, v=2, t0= 1, d=0, sz = 0, sv = 0, st0 = 1)
plot(rts, y_st0)
plot(rts, y_st1)

prd(1, a=1, z=0.5, v=2, t0=0.1, d=0.1, sz = 0, sv = 0, st0 = 10)
