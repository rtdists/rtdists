context("Censored shifted Wald distribution functions work correctly")

## The censored shifted Wald (Miller, Scherbaum, Heck, Goschke, & Enge, 2018)
## is a competing-risks race of two shifted Wald accumulators with drift rates
## v and -v. Its functions are checked against (a) the closed-form equations
## of the paper evaluated in its notation, (b) statmod plus the defective
## factorization as an independent reference, and (c) each other (integration,
## inversion, simulation). The internal switch to the "simple" censored
## shifted Wald (Eq. 4) is checked against its analytic error bound.

## shifted Wald density and survival in Miller et al.'s notation (their
## Eqs. 2 and 3 with gamma = threshold, delta = drift, theta = 0), valid for
## both signs of delta
f_sw <- function(t, gamma, delta)
  gamma/sqrt(2*pi*t^3) * exp(-(gamma - delta*t)^2/(2*t))
S_sw <- function(t, gamma, delta)
  1 - (pnorm((delta*t - gamma)/sqrt(t)) +
         exp(2*delta*gamma)*pnorm((-delta*t - gamma)/sqrt(t)))

test_that("dcswald agrees with Miller et al. (2018) Eq. 5 for both responses", {
  t <- seq(0.05, 4, length.out = 60)
  b <- 0.9; v <- 1.4; t0 <- 0.25

  ## unbiased start: both thresholds equal b
  expect_equal(dcswald(t + t0, "upper", b = b, t0 = t0, v = v),
               f_sw(t, b, v) * S_sw(t, b, -v))
  expect_equal(dcswald(t + t0, "lower", b = b, t0 = t0, v = v),
               f_sw(t, b, -v) * S_sw(t, b, v))

  ## biased start: b_upper = 2*b*(1 - w), b_lower = 2*b*w
  w <- 0.65
  expect_equal(dcswald(t + t0, "upper", b = b, t0 = t0, v = v, w = w),
               f_sw(t, 2*b*(1 - w), v) * S_sw(t, 2*b*w, -v))
  expect_equal(dcswald(t + t0, "lower", b = b, t0 = t0, v = v, w = w),
               f_sw(t, 2*b*w, -v) * S_sw(t, 2*b*(1 - w), v))

  ## negative drift rate favors the lower response
  expect_equal(dcswald(t + t0, "lower", b = b, t0 = t0, v = -v),
               f_sw(t, b, v) * S_sw(t, b, -v))

  ## s rescales the evidence scale (Ratcliff convention s = 0.1)
  s <- 0.1
  expect_equal(dcswald(t + t0, "upper", b = b*s, t0 = t0, v = v*s, s = s),
               f_sw(t, b, v) * S_sw(t, b, -v))

  ## log = TRUE does not lose the far tail to underflow
  ld <- dcswald(2000, "upper", b = b, t0 = t0, v = v, log = TRUE)
  expect_true(is.finite(ld) && ld < -1000)
  expect_identical(dcswald(2000, "upper", b = b, t0 = t0, v = v), 0)
})

test_that("the b and a parameterizations are identical with a = 2*b", {
  t <- c(0.4, 0.8, 1.5)
  for (resp in c("upper", "lower")) {
    expect_identical(dcswald(t, resp, b = 0.8, t0 = 0.2, v = 1.2, w = 0.6),
                     dcswald(t, resp, a = 1.6, t0 = 0.2, v = 1.2, w = 0.6))
    expect_identical(pcswald(t, resp, b = 0.8, t0 = 0.2, v = 1.2, w = 0.6),
                     pcswald(t, resp, a = 1.6, t0 = 0.2, v = 1.2, w = 0.6))
  }
  expect_identical(
    qcswald(c(0.2, 0.6), "upper", b = 0.8, t0 = 0.2, v = 1.2),
    qcswald(c(0.2, 0.6), "upper", a = 1.6, t0 = 0.2, v = 1.2))
})

test_that("the race is proper: the response probabilities sum to one", {
  ## the grid of Miller et al.-style parameter sets used when verifying the
  ## branch: gamma_u x gamma_l x delta, realized via the a parameterization
  for (gu in c(0.5, 1.2, 3)) for (gl in c(0.5, 1, 2.5))
    for (delta in c(0.5, 2, 4)) {
      a <- gu + gl
      w <- gl/a
      p <- pcswald(Inf, "upper", a = a, t0 = 0.2, v = delta, w = w) +
        pcswald(Inf, "lower", a = a, t0 = 0.2, v = delta, w = w)
      expect_equal(p, 1, tolerance = 1e-8)
    }
})

test_that("the censored shifted Wald reduces to the shifted Wald when the losing accumulator never finishes", {
  ## Miller et al. (2018, p. 123): with nothing censored, Eq. 4 is the plain
  ## shifted Wald likelihood. Internally, once the losing accumulator's total
  ## finishing probability exp(-2*k*|v|/s^2) is below one ulp its survival
  ## factor is dropped, and d/p/qcswald are *identical* to d/p/qwald (A = 0).
  b <- 3; v <- 6.5; t0 <- 0.3  # 2*b*v = 39 > -log(eps/2) = 36.74
  rt <- c(0.35, 0.5, 0.8, 1.5, 5)
  expect_identical(dcswald(rt, "upper", b = b, t0 = t0, v = v),
                   dwald(rt, b = b, t0 = t0, v = v))
  expect_identical(pcswald(rt, "upper", b = b, t0 = t0, v = v),
                   pwald(rt, b = b, t0 = t0, v = v))
  expect_identical(qcswald(c(0.1, 0.5, 0.9), "upper", b = b, t0 = t0, v = v),
                   qwald(c(0.1, 0.5, 0.9), b = b, t0 = t0, v = v))

  ## the same identity holds through the s scaling
  expect_identical(
    dcswald(rt, "upper", b = b/10, t0 = t0, v = v/10, s = 0.1),
    dwald(rt, b = b/10, t0 = t0, v = v/10, s = 0.1))
})

test_that("dropping the survival factor is exact at the switch point", {
  ## the relative error of the simple branch equals the losing accumulator's
  ## total mass exp(-delta), delta = 2*k*|v|/s^2 (see the implementation
  ## notes); directly above the cutoff the two branches agree to < 1 ulp.
  t <- exp(seq(log(0.01), log(20), length.out = 100))
  for (delta in c(10, 20, 30, 36, 37, 39)) {
    v <- delta/2  # k_lower = 1 at b = 1, w = 0.5... realized via a = 2, w = 0.5
    crisk <- f_sw(t, 1, v) * S_sw(t, 1, -v)
    got <- dcswald(t + 0.2, "upper", a = 2, t0 = 0.2, v = v)
    expect_equal(got, crisk, tolerance = max(4e-15, 2 * exp(-delta)),
                 info = paste("delta =", delta))
  }
})

test_that("the defective functions match statmod via the factorization", {
  ## F_defective(t | k, v < 0) = exp(2*k*v) * F_IG(t | mu = k/|v|, k^2):
  ## conditional on finishing, the losing accumulator is the ordinary inverse
  ## Gaussian of its mirrored drift. This checks the negative-drift survival
  ## factor against statmod, which never sees a negative drift.
  skip_if_not_installed("statmod")
  t <- seq(0.05, 20, length.out = 80)
  for (k in c(0.5, 1, 2.5)) for (v in c(0.5, 2, 4)) {
    S_def <- 1 - exp(-2*k*v) *
      statmod::pinvgauss(t, mean = k/v, shape = k^2)
    expect_equal(
      dcswald(t + 0.1, "upper", a = 2*k, t0 = 0.1, v = v) /
        statmod::dinvgauss(t, mean = k/v, shape = k^2),
      S_def, tolerance = 1e-12)
  }
})

test_that("pcswald is the integral of dcswald, over both responses", {
  b <- 0.8; v <- 1.2; t0 <- 0.25; w <- 0.4
  for (resp in c("upper", "lower")) {
    for (rt in c(0.6, 1.2, Inf)) {
      num <- integrate(dcswald, lower = t0, upper = rt, response = resp,
                       b = b, t0 = t0, v = v, w = w, rel.tol = 1e-10)$value
      expect_equal(pcswald(rt, resp, b = b, t0 = t0, v = v, w = w), num,
                   tolerance = 1e-8)
    }
  }
})

test_that("qcswald inverts pcswald", {
  b <- 0.8; v <- 1.2; t0 <- 0.25
  ## defective probabilities up to the response probability
  for (resp in c("upper", "lower")) {
    mass <- pcswald(Inf, resp, b = b, t0 = t0, v = v)
    p <- c(0.1, 0.5, 0.9) * mass
    q <- qcswald(p, resp, b = b, t0 = t0, v = v)
    expect_equal(pcswald(q, resp, b = b, t0 = t0, v = v), p,
                 tolerance = 1e-6)
    ## scale_p rescales by the response probability
    expect_equal(qcswald(c(0.1, 0.5, 0.9), resp, b = b, t0 = t0, v = v,
                         scale_p = TRUE), q, tolerance = 1e-6)
  }
  ## p = 0 gives t0, probabilities above the response probability give NA
  expect_identical(qcswald(0, "upper", b = b, t0 = t0, v = v), t0)
  expect_warning(
    out <- qcswald(0.5, "lower", b = b, t0 = t0, v = v),
    "predicted response probability")
  expect_identical(out, NA_real_)
})

test_that("the sign of the drift rate mirrors the race", {
  t <- c(0.5, 0.9, 2)
  ## flipping the sign of v and the bias swaps the two responses (w chosen
  ## so that 1 - w is exactly representable)
  expect_identical(dcswald(t, "upper", b = 1, t0 = 0.2, v = 1.5, w = 0.25),
                   dcswald(t, "lower", b = 1, t0 = 0.2, v = -1.5, w = 0.75))
  expect_identical(pcswald(t, "upper", b = 1, t0 = 0.2, v = 1.5, w = 0.25),
                   pcswald(t, "lower", b = 1, t0 = 0.2, v = -1.5, w = 0.75))
  ## v = 0 with an unbiased start gives equal response probabilities
  expect_equal(pcswald(Inf, "upper", b = 1, t0 = 0.2, v = 0), 0.5,
               tolerance = 1e-8)
})

test_that("rcswald simulates the race", {
  skip_on_cran()
  set.seed(4)
  b <- 0.8; v <- 1.2; t0 <- 0.25
  r <- rcswald(4e4, b = b, t0 = t0, v = v)
  expect_true(is.data.frame(r) && all(c("rt", "response") %in% names(r)))
  expect_identical(levels(r$response), c("lower", "upper"))
  expect_true(all(is.finite(r$rt)) && all(r$rt > t0))

  ## response proportions follow the sign of the drift rate
  p_up <- pcswald(Inf, "upper", b = b, t0 = t0, v = v)
  expect_equal(mean(r$response == "upper"), p_up, tolerance = 0.01)
  set.seed(5)
  r_neg <- rcswald(4e4, b = b, t0 = t0, v = -v)
  expect_equal(mean(r_neg$response == "lower"), p_up, tolerance = 0.01)

  ## the RTs of each response follow the conditional defective distribution
  ## (subsampled: below the switch point every pcswald value is a quadrature)
  for (resp in c("upper", "lower")) {
    mass <- pcswald(Inf, resp, b = b, t0 = t0, v = v)
    rts <- r$rt[r$response == resp]
    rts <- rts[seq_len(min(length(rts), 2000))]
    ks <- suppressWarnings(ks.test(
      rts, function(q) pcswald(q, resp, b = b, t0 = t0, v = v)/mass))
    expect_true(ks$statistic < 0.04)
  }
})

test_that("the defective density is the censored likelihood of a data set", {
  ## sum(dcswald(data, ..., log = TRUE)) is Eq. 5's log-likelihood; the
  ## data.frame convenience input matches the explicit call
  skip_on_cran()
  set.seed(6)
  r <- rcswald(500, b = 0.8, t0 = 0.25, v = 1.2)
  ll_df <- dcswald(r, b = 0.8, t0 = 0.25, v = 1.2, log = TRUE)
  ll_expl <- dcswald(r$rt, r$response, b = 0.8, t0 = 0.25, v = 1.2,
                     log = TRUE)
  expect_identical(ll_df, ll_expl)
  expect_true(all(is.finite(ll_df)))
})

test_that("all parameters are recycled and responses can be mixed", {
  rt <- c(0.5, 0.9, 1.4, 2)
  resp <- c("upper", "lower", "upper", "lower")
  one_by_one <- vapply(seq_along(rt), function(i)
    dcswald(rt[i], resp[i], b = 0.9, t0 = 0.2, v = 1.1, w = 0.55),
    numeric(1))
  expect_identical(dcswald(rt, resp, b = 0.9, t0 = 0.2, v = 1.1, w = 0.55),
                   one_by_one)
  ## numeric responses: 1 = lower, 2 = upper
  expect_identical(dcswald(rt, c(2, 1, 2, 1), b = 0.9, t0 = 0.2, v = 1.1,
                           w = 0.55), one_by_one)
  ## vectorized parameters
  expect_identical(
    dcswald(rep(0.9, 3), "upper", b = c(0.7, 0.9, 1.1), t0 = 0.2, v = 1.1),
    vapply(c(0.7, 0.9, 1.1), function(b)
      dcswald(0.9, "upper", b = b, t0 = 0.2, v = 1.1), numeric(1)))
})

test_that("dcswald and pcswald are 0 below t0 and dcswald is 0 at infinity", {
  expect_identical(dcswald(c(0.1, 0.2), "upper", b = 1, t0 = 0.2, v = 2),
                   c(0, 0))
  expect_identical(pcswald(c(0.1, 0.2), "upper", b = 1, t0 = 0.2, v = 2),
                   c(0, 0))
  expect_identical(dcswald(Inf, "upper", b = 1, t0 = 0.2, v = 2), 0)
  ## deep into the tail nothing overflows, even in the competing-risks branch
  ld <- dcswald(c(1e3, 1e5), "upper", b = 1, t0 = 0.2, v = 2, log = TRUE)
  expect_true(all(is.finite(ld)))
})

test_that("the negative-drift CDF no longer overflows (pig_log regression)", {
  ## before the log-sum-exp fix, the lower tail of pig_log overflowed to Inf
  ## for a negative mu once rt was large; the total mass is exp(2*k*v)
  lp <- rtdists:::pig_log(c(20, 1e3, 1e5, 1e8), rep(-0.5, 4), rep(1, 4))
  expect_equal(exp(lp), rep(exp(-4), 4))
  ## and at infinity both tails give the defective mass
  expect_equal(exp(rtdists:::pig_log(Inf, -0.5, 1)), exp(-4))
  expect_equal(exp(rtdists:::pig_log(Inf, -0.5, 1, lower.tail = FALSE)),
               1 - exp(-4))
})

test_that("input is checked", {
  expect_error(dcswald(0.5, "upper", t0 = 0.2, v = 1),
               "exactly one of b")
  expect_error(dcswald(0.5, "upper", b = 1, a = 2, t0 = 0.2, v = 1),
               "exactly one of b")
  expect_error(dcswald(0.5, "upper", b = 1, t0 = 0.2, v = 1, w = 1),
               "strictly between")
  expect_error(dcswald(0.5, "upper", b = 1, t0 = 0.2, v = 1, w = 0),
               "strictly between")
  expect_error(dcswald(0.5, "upper", b = -1, t0 = 0.2, v = 1),
               "positive")
  expect_error(dcswald(0.5, "upper", b = 1, t0 = 0.2, v = 1, s = 0),
               "positive")
  expect_error(dcswald(0.5, "middle", b = 1, t0 = 0.2, v = 1))
  expect_error(dcswald(0.5, 3, b = 1, t0 = 0.2, v = 1), "response")
  expect_error(dcswald(-0.5, "upper", b = 1, t0 = 0.2, v = 1),
               "positive values")
  expect_error(rcswald(c(10, 10), b = 1, t0 = 0.2, v = 1), "length 1")
})
