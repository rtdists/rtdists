
## the shifted Wald distribution is the default (i.e., A = 0), which
## corresponds to a single evidence accumulation process without start point
## variability:
rwald(10, b = 1, t0 = 0.3, v = 2)

# use somewhat plausible values for plotting:
b <- 0.5
t0 <- 0.3

# plot density of the shifted Wald and the effect of start point variability:
curve(dwald(x, b = b, t0 = t0, v = 2), ylim = c(0, 4), xlim = c(0, 3),
      main = "Density/PDF of the Wald accumulator", ylab = "density",
      xlab = "response time", lty = 2)
curve(dwald(x, A = 0.2, b = b, t0 = t0, v = 2), add = TRUE)
legend("topright", legend = c("0.2", "0"),
       title = expression("Start point variability"~~italic(A)), lty = c(1, 2))

# plot cdf:
curve(pwald(x, b = b, t0 = t0, v = 2), xlim = c(0, 3), ylim = c(0, 1),
      ylab = "cumulative probability", xlab = "response time", lty = 2,
      main = "Distribution/CDF of the Wald accumulator")
curve(pwald(x, A = 0.2, b = b, t0 = t0, v = 2), add = TRUE)
legend("bottomright", legend = c("0.2", "0"),
       title = expression("Start point variability"~~italic(A)), lty = c(1, 2))

## the quantile function inverts the CDF:
q <- qwald(c(0.1, 0.3, 0.5, 0.7, 0.9), b = b, t0 = t0, v = 2)
q
pwald(q, b = b, t0 = t0, v = 2)  # 0.1, 0.3, 0.5, 0.7, 0.9

## variability in non-decision time (only available for A = 0):
curve(dwald(x, b = b, t0 = t0, v = 2), xlim = c(0, 3), ylim = c(0, 4),
      main = "Shifted Wald with variable non-decision time",
      ylab = "density", xlab = "response time")
curve(dwald(x, b = b, t0 = t0, v = 2, st0 = 0.3), add = TRUE, lty = 2)
legend("topright", legend = c("0", "0.3"),
       title = expression("Non-decision time variability"~~italic(s)[t0]),
       lty = 1:2)

## log densities are calculated on the log scale throughout, so likelihoods
## remain usable far into the tail:
dwald(c(0.31, 5), b = b, t0 = t0, v = 2, log = TRUE)

\dontrun{
## parameter recovery for the shifted Wald:
set.seed(1)
rts <- rwald(500, b = 0.5, t0 = 0.3, v = 2)

objective_fun <- function(pars, rt) {
  -sum(dwald(rt, b = pars[1], t0 = pars[2], v = pars[3], log = TRUE))
}
nlminb(c(1, 0.1, 1), objective_fun, rt = rts,
       lower = c(0.01, 0.01, 0.01), upper = c(10, min(rts), 10))$par
}
