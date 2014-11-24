

## check random generated values against race functions:

## 1. Without st0:
r_lba <- rlba_norm(1e4, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=0.2)
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
normalised_n1CDF = function(t,...) n1CDF(t,...)/n1CDF(t=Inf,...) 
ks.test(r_lba$rt[r_lba$response==1], normalised_n1CDF, A=0.5, b=1, t0 = 0.5, 
        mean_v=c(1.2, 1.0), sd_v=0.2)
}

## 2. With st0 = 0.2:
r_lba <- rlba_norm(1e4, A=0.5, b=1, t0 = 0.5, st0 = 0.2, mean_v=c(1.2, 1), sd_v=0.2)
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
\dontrun{
normalised_n1CDF = function(t,...) n1CDF(t,...)/n1CDF(t=Inf,...) 
ks.test(r_lba$rt[r_lba$response==1], normalised_n1CDF, A=0.5, b=1, t0 = 0.5, 
        st0 = 0.2, mean_v=c(1.2, 1.0), sd_v=0.2)
}

## Other examples:

xx <- rlba_norm(10, A=0.5, b=1, t0 = 0.5, mean_v=1.2, sd_v=0.2)
n1PDF(x, A=0.5, b=1, t0 = 0.5, st0 = 0.2, mean_v=c(1.2, 1.0), sd_v=0.2)

n1PDF(xx$rt, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1.0), sd_v=0.2, st0 = 0.4)

n1PDF(xx$rt, A=0.5, b=1, t0 = 0.5, 
      meanlog_v = c(0.5, 0.8), sdlog_v = 0.5, distribution = "lnorm")
n1PDF(xx$rt, A=0.5, b=1, t0 = 0.5, 
      meanlog_v = c(0.5, 0.8), sdlog_v = 0.5, st0 = 0.2, distribution = "lnorm")

n1CDF(xx$rt, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1.0), sd_v=0.2)
n1CDF(xx$rt, A=0.5, b=1, t0 = 0.5, mean_v=c(1.0, 1.2), sd_v=0.2)

n1CDF(xx$rt, A=0.5, b=1, t0 = 0.5, meanlog_v = c(0.5, 0.8), sdlog_v = 0.5, distribution = "lnorm")
n1CDF(xx$rt, A=0.5, b=1, t0 = 0.5, meanlog_v = c(0.8, 0.5), sdlog_v = 0.5, distribution = "lnorm")