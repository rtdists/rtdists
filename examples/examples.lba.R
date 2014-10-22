
## random number generation using different distributions for v:
rlba(10, A=0.5, b=1, v=c(1.2, 1), sv=c(1,1.5), t0 = 0.5)  # default uses "normal"
rlba(10, A=0.5, b=1, v=c(1.2, 1), sv=c(1, 1.5), t0 = 0.5, v_distribution = "gamma")
rlba(10, A=0.5, b=1, v=c(1.2, 1), sv=c(1, 1.5), t0 = 0.5, v_distribution = "frechet")
rlba(10, A=0.5, b=1, v=c(1.2, 1), sv=c(1, 1.5), t0 = 0.5, v_distribution = "lognormal")


## make nice distribution plots
# use somewhat plausible values for plotting:
A <- 0.2
b <- 0.5
v <- 1.0
t0 <- 0.3
sv <- 1.1  # note that frechet needs sv > 1.

# get all possible distributions
(v_distributions <- eval(formals(dlba)$v_distribution))

# plot density (pdf):
curve(x-10, ylim = c(0, 5), xlim=c(0,3), 
      main="Density of LBA versions", ylab="density", xlab="response time")
for (i in seq_along(v_distributions)) 
  curve(dlba(x, A=A, b=b, v=v, t0=t0, sv=sv,v_distribution = v_distributions[i]), add=TRUE, lty = i)
legend("topright", legend=v_distributions, title = expression("Distribution of"~~italic(v)), 
       lty = seq_along(v_distributions))

# plot cdf:
curve(x-10, xlim = c(0, 3), ylim=c(0,1), 
      main="CDF of LBA versions", ylab="cumulative probability", xlab="response time")
for (i in seq_along(v_distributions)) 
  curve(plba(x, A=A, b=b, v=v, t0=t0, sv=sv,v_distribution = v_distributions[i]), add=TRUE, lty = i)
legend("bottomright", legend=v_distributions, title = expression("Distribution of"~~italic(v)), 
       lty = seq_along(v_distributions))

