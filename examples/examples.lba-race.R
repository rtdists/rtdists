

## check random generated values against race functions:

## 1. Without st0:
r_lba <- rLBA(1e4, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=0.2)
x <- seq(0.5, 4, length.out = 100) # for plotting
# PDF
y <- n1PDF(x, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1.0), sd_v=0.2) # PDF
hist(r_lba$rt[r_lba$response==1],probability = TRUE, breaks = "FD")
lines(x=x,y=y/mean(r_lba$response == 1))
# CDF
plot(ecdf(r_lba$rt[r_lba$response==1]))
y <- n1CDF(x, A=0.5, b=1, t0 = 0.5, st0 = 0, mean_v=c(1.2, 1.0), sd_v=0.2)
lines(x=x,y=y/mean(r_lba$response == 1), col = "red", lwd = 4.5, lty = 2)
# KS test
\dontrun{
normalised_n1CDF = function(rt,...) n1CDF(rt,...)/n1CDF(rt=Inf,...) 
ks.test(r_lba$rt[r_lba$response==1], normalised_n1CDF, A=0.5, b=1, t0 = 0.5, 
        mean_v=c(1.2, 1.0), sd_v=0.2)
}

\dontrun{
## Other examples (don't run to save time):
  
## 2. With st0 = 0.2:
r_lba <- rLBA(1e4, A=0.5, b=1, t0 = 0.5, st0 = 0.2, mean_v=c(1.2, 1), sd_v=0.2)
x <- seq(0.5, 4, length.out = 100) # for plotting
# PDF
y <- n1PDF(x, A=0.5, b=1, t0 = 0.5, st0 = 0.2, mean_v=c(1.2, 1.0), sd_v=0.2) # PDF
hist(r_lba$rt[r_lba$response==1],probability = TRUE, breaks = "FD")
lines(x=x,y=y/mean(r_lba$response == 1))
# CDF
plot(ecdf(r_lba$rt[r_lba$response==1]))
y <- n1CDF(x, A=0.5, b=1, t0 = 0.5, st0 = 0.2, mean_v=c(1.2, 1.0), sd_v=0.2)
lines(x=x,y=y/mean(r_lba$response == 1), col = "red", lwd = 4.5, lty = 2)
# KS test
normalised_n1CDF = function(rt,...) n1CDF(rt,...)/n1CDF(rt=Inf,...) 
ks.test(r_lba$rt[r_lba$response==1], normalised_n1CDF, A=0.5, b=1, t0 = 0.5, 
        st0 = 0.2, mean_v=c(1.2, 1.0), sd_v=0.2)


xx <- rLBA(10, A=0.5, b=1, t0 = 0.5, mean_v=1.2, sd_v=0.2)

# default uses normal distribution for drift rate:
n1PDF(xx$rt, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1.0), sd_v=0.2)

# other distributions:
n1PDF(xx$rt, A=0.5, b=1, t0 = 0.5, shape_v=c(1.2, 1), scale_v=c(0.2,0.3), distribution = "gamma")
n1PDF(xx$rt, A=0.5, b=1, t0 = 0.5, shape_v=c(1.2, 1), scale_v=c(0.2,0.3), distribution = "frechet")
n1PDF(xx$rt, A=0.5, b=1, t0 = 0.5, meanlog_v = c(0.5, 0.8), sdlog_v = 0.5, distribution = "lnorm")

# add st0:
n1PDF(xx$rt, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1.0), sd_v=0.2, st0 = 0.4)


# use different A parameters for each RT:
n1PDF(xx$rt, A=runif(10, 0.4, 0.6), 
      b=1, t0 = 0.5, mean_v=c(1.2, 1.0), sd_v=0.2)

# use different A parameters for each RT and each accumulator:
n1PDF(xx$rt, A=list(runif(10, 0.4, 0.6), runif(10, 0.2, 0.4)), 
      b=1, t0 = 0.5, mean_v=c(1.2, 1.0), sd_v=0.2)


### vectorize drift rates:

# vector versus list:
v1 <- n1PDF(xx$rt, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1.0), sd_v=0.2)
v2 <- n1PDF(xx$rt, A=0.5, b=1, t0 = 0.5, mean_v=list(1.2, 1.0), sd_v=0.2)
identical(v1, v2)  # TRUE

# drift rate per trial:
n1PDF(xx$rt, A=0.5, b=1, t0 = 0.5, mean_v=list(rnorm(10, 1.2), rnorm(10, 1)), sd_v=0.2)

# combine list with vector:
n1PDF(xx$rt, A=0.5, b=1, t0 = 0.5, mean_v=list(rnorm(10, 1.2), rnorm(10, 1)), sd_v=c(0.2, 0.1))

# t0 per trial and accumulator:
n1PDF(xx$rt, A=0.5, b=1, t0 = c(0.5), mean_v=c(1.2, 1.0), sd_v=0.2)
n1PDF(xx$rt, A=0.5, b=1, t0 = c(0.5, 0.6), mean_v=c(1.2, 1.0), sd_v=0.2) # per trial only
n1PDF(xx$rt, A=0.5, b=1, t0 = list(0.5, 0.6), mean_v=c(1.2, 1.0), sd_v=0.2) # per drift rate only
n1PDF(xx$rt, A=0.5, b=1, t0 = list(c(0.4, 0.5), c(0.5, 0.6)), mean_v=c(1.2, 1.0), sd_v=0.2)
}
