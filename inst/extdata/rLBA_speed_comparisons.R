
set.seed(1)
rt1 <- rLBA(500, A=0.5, b=1, t0 = 0.5, mean_v=
              list(c(2.7, seq(1.6, 2.6, length.out = 29)), 
                   seq(0.6, 1.6, length.out = 30)), sd_v=c(1,1.2))
head(rt1)

rt2 <- rLBA(500, A=0.5, b=1, t0 = 0.5, mean_v=list(2.6, 1.6), sd_v=c(1,1.2))

library(profvis)

profvis({
replicate(50, rLBA(500, A=0.5, b=1, t0 = 0.5, mean_v=list(seq(1.6, 2.6, length.out = 500), seq(0.6, 1.6, length.out = 500)), sd_v=c(1,1.2), silent = TRUE))
})

str(rt3)

head(rt1)
str(rt1)

require(microbenchmark)
microbenchmark(rLBA(500, A=0.5, b=1, t0 = 0.5, mean_v= list(seq(1.6, 2.6, length.out = 20), seq(0.6, 1.6, length.out = 20)), sd_v=c(1,1.2), silent = TRUE),
              rLBA(500, A=0.5, b=1, t0 = 0.5, mean_v=list(seq(1.6, 2.6, length.out = 500), seq(0.6, 1.6, length.out = 500)), sd_v=c(1,1.2), silent = TRUE),
              rLBA(500, A=0.5, b=1, t0 = 0.5, mean_v=list(2.6, 1.6), sd_v=c(1,1.2), silent = TRUE))
