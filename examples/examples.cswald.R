## censored shifted Wald: two shifted Wald accumulators race, the loser
## censors the winner. Unbiased start point, single bound b:
curve(dcswald(x, response = "upper", b = 1, t0 = 0.3, v = 1.5),
      xlim = c(0, 3), ylim = c(0, 1.5), ylab = "defective density",
      xlab = "t", main = "censored shifted Wald (b = 1, v = 1.5, t0 = 0.3)")
curve(dcswald(x, response = "lower", b = 1, t0 = 0.3, v = 1.5),
      add = TRUE, lty = 2)
legend("topright", legend = c("upper", "lower"), lty = 1:2)

## the two response probabilities sum to 1 (the race is always proper):
pu <- pcswald(Inf, response = "upper", b = 1, t0 = 0.3, v = 1.5)
pl <- pcswald(Inf, response = "lower", b = 1, t0 = 0.3, v = 1.5)
c(upper = pu, lower = pl, sum = pu + pl)

## the same model in the Wiener-style parameterization (a = 2 * b), here
## with a start point biased towards the upper bound:
dcswald(0.9, response = "upper", a = 2, t0 = 0.3, v = 1.5, w = 0.6)

## when drift rate and the disfavored bound are high the lower accumulator
## (drift -v) practically never finishes and the censored shifted Wald
## reduces to the plain shifted Wald (a single Wald accumulator with A = 0):
dcswald(0.9, response = "upper", b = 3, t0 = 0.3, v = 6.5) -
  dwald(0.9, A = 0, b = 3, t0 = 0.3, v = 6.5)

## random generation simulates the race; response proportions follow the
## sign of the drift rate:
rts <- rcswald(1e4, b = 1, t0 = 0.3, v = 1.5)
prop.table(table(rts$response))
pcswald(Inf, response = "lower", b = 1, t0 = 0.3, v = 1.5)

\dontrun{
## ---------------------------------------------------------------------------
## maximum likelihood estimation, and the bias from ignoring the censoring:
## fitting only the correct (upper) RTs with a plain shifted Wald treats the
## error-censored trials as if they never happened and biases the estimates
## (Miller et al., 2018, Fig. 2); the censored shifted Wald likelihood
## accounts for them.
set.seed(1)
true <- c(b = 0.8, t0 = 0.25, v = 1.2)  # ~7% errors
dat <- rcswald(5e3, b = true["b"], t0 = true["t0"], v = true["v"])

ll_csw <- function(par, dat) {
  d <- dcswald(dat, b = par[1], t0 = par[2], v = par[3], log = TRUE)
  if (any(!is.finite(d))) return(1e10)
  -sum(d)
}
ll_sw <- function(par, rt) {
  d <- log(dwald(rt, A = 0, b = par[1], t0 = par[2], v = par[3]))
  if (any(!is.finite(d))) return(1e10)
  -sum(d)
}
start <- c(1, 0.1, 1)
fit_csw <- nlminb(start, ll_csw, dat = dat,
                  lower = c(0.05, 0.01, 0.05))
fit_sw <- nlminb(start, ll_sw, rt = dat$rt[dat$response == "upper"],
                 lower = c(0.05, 0.01, 0.05))
rbind(true = true, cswald = fit_csw$par, plain_sw = fit_sw$par)
}
