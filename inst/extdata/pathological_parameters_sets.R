
require(rtdists)

rlba_norm(10, A=.89, b=.69, t0 = .52, mean_v=c(.87, .84), sd_v=c(.43, .18), posdrift = TRUE)

source("lba-math.r")

rlba(100, A=.89, b=.69, t0 = .52, vs=c(.87, .84), s=c(.43, .18), posdrift = TRUE)
