## The censored shifted Wald (Miller, Scherbaum, Heck, Goschke, & Enge, 2018)
## is a competing-risks race of two shifted Wald accumulators with drift rates
## v and -v. Its functions are checked against (a) the closed-form equations
## of the paper evaluated in its notation, (b) numerical integration of those
## equations, and (c) each other (integration, inversion, simulation).

## shifted Wald density and CDF in Miller et al.'s notation (their Eqs. 2 and
## 3 with gamma = threshold, delta = drift, theta = 0), valid for both signs
## of delta
f_sw <- function(t, gamma, delta)
  gamma/sqrt(2*pi*t^3) * exp(-(gamma - delta*t)^2/(2*t))
F_sw <- function(t, gamma, delta)
  pnorm((delta*t - gamma)/sqrt(t)) +
    exp(2*delta*gamma)*pnorm((-delta*t - gamma)/sqrt(t))
S_sw <- function(t, gamma, delta) 1 - F_sw(t, gamma, delta)
## Eq. 5: defective densities of the upper (threshold gu, drift delta) and
## lower (threshold gl, drift -delta) response on the decision-time scale
f_up <- function(t, gu, gl, delta) f_sw(t, gu, delta) * S_sw(t, gl, -delta)
f_low <- function(t, gu, gl, delta) f_sw(t, gl, -delta) * S_sw(t, gu, delta)
## the corresponding defective CDFs by quadrature (independent of pcswald)
P_num <- function(x, resp, gu, gl, delta) {
  fun <- if (resp == "upper") f_up else f_low
  vapply(x, function(u)
    integrate(fun, lower = 0, upper = u, gu = gu, gl = gl, delta = delta,
              rel.tol = 1e-10)$value, numeric(1))
}

test_that("dcswald agrees with Miller et al. (2018) Eq. 5 for both responses", {
  t <- seq(0.05, 4, length.out = 60)
  b <- 0.9; v <- 1.4; t0 <- 0.25

  ## unbiased start: both thresholds equal b
  expect_equal(dcswald(t + t0, "upper", b = b, t0 = t0, v = v),
               f_up(t, b, b, v))
  expect_equal(dcswald(t + t0, "lower", b = b, t0 = t0, v = v),
               f_low(t, b, b, v))

  ## biased start: b_upper = 2*b*(1 - w), b_lower = 2*b*w
  w <- 0.65
  expect_equal(dcswald(t + t0, "upper", b = b, t0 = t0, v = v, w = w),
               f_up(t, 2*b*(1 - w), 2*b*w, v))
  expect_equal(dcswald(t + t0, "lower", b = b, t0 = t0, v = v, w = w),
               f_low(t, 2*b*(1 - w), 2*b*w, v))

  ## negative drift rate favors the lower response
  expect_equal(dcswald(t + t0, "lower", b = b, t0 = t0, v = -v),
               f_up(t, b, b, v))

  ## s rescales the evidence scale (Ratcliff convention s = 0.1)
  s <- 0.1
  expect_equal(dcswald(t + t0, "upper", b = b*s, t0 = t0, v = v*s, s = s),
               f_up(t, b, b, v))

  ## log = TRUE returns the log density
  expect_equal(dcswald(t + t0, "upper", b = b, t0 = t0, v = v, log = TRUE),
               log(f_up(t, b, b, v)))
})

test_that("the b and a parameterizations are identical with a = 2*b", {
  t <- c(0.4, 0.8, 1.5)
  for (resp in c("upper", "lower")) {
    expect_identical(dcswald(t, resp, b = 0.8, t0 = 0.2, v = 1.2, w = 0.6),
                     dcswald(t, resp, a = 1.6, t0 = 0.2, v = 1.2, w = 0.6))
    expect_identical(pcswald(t, resp, b = 0.8, t0 = 0.2, v = 1.2, w = 0.6),
                     pcswald(t, resp, a = 1.6, t0 = 0.2, v = 1.2, w = 0.6))
  }
})

test_that("pcswald has a closed form for an unbiased start point", {
  ## With equal thresholds k, the losing accumulator's CDF is c*F with
  ## c = exp(-2*k*|v|) and F the CDF of the favored accumulator, so Eq. 5
  ## integrates to F - c*F^2/2 (favored) and c*(F - F^2/2) (disfavored
  ## response). Checked against quadrature of the paper's equations.
  b <- 0.9; v <- 1.4; t0 <- 0.25
  t <- c(0.1, 0.4, 0.8, 1.5, 3)
  for (resp in c("upper", "lower")) {
    num <- P_num(t, resp, b, b, v)
    tot <- P_num(Inf, resp, b, b, v)
    expect_equal(pcswald(t + t0, resp, b = b, t0 = t0, v = v), num,
                 tolerance = 1e-8)
    expect_equal(pcswald(Inf, resp, b = b, t0 = t0, v = v), tot,
                 tolerance = 1e-8)
    expect_equal(pcswald(t + t0, resp, b = b, t0 = t0, v = v,
                         lower.tail = FALSE), tot - num, tolerance = 1e-8)
    expect_identical(pcswald(Inf, resp, b = b, t0 = t0, v = v,
                             lower.tail = FALSE), 0)
  }
  ## the response probabilities are exact
  cc <- exp(-2*b*v)
  expect_equal(pcswald(Inf, "upper", b = b, t0 = t0, v = v), 1 - cc/2)
  expect_equal(pcswald(Inf, "lower", b = b, t0 = t0, v = v), cc/2)
  ## the same through the s scaling
  s <- 0.1
  expect_equal(pcswald(t + t0, "lower", b = b*s, t0 = t0, v = v*s, s = s),
               P_num(t, "lower", b, b, v), tolerance = 1e-8)
})

test_that("pcswald integrates dcswald for a biased start point", {
  ## unequal thresholds have no closed form; the defective density is
  ## integrated numerically and checked against quadrature of Eq. 5
  b <- 0.8; v <- 1.2; t0 <- 0.25; w <- 0.4
  gu <- 2*b*(1 - w); gl <- 2*b*w
  for (resp in c("upper", "lower")) {
    tot <- P_num(Inf, resp, gu, gl, v)
    for (rt in c(0.6, 1.2, Inf)) {
      num <- P_num(rt - t0, resp, gu, gl, v)
      expect_equal(pcswald(rt, resp, b = b, t0 = t0, v = v, w = w), num,
                   tolerance = 1e-8)
      expect_equal(pcswald(rt, resp, b = b, t0 = t0, v = v, w = w,
                           lower.tail = FALSE), tot - num, tolerance = 1e-8)
    }
  }
})

test_that("the race is proper: the response probabilities sum to one", {
  ## grid of Miller et al.-style parameter sets, gamma_u x gamma_l x delta,
  ## realized via the a parameterization (w = 0.5 takes the closed form,
  ## everything else the numerical integration)
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
  ## shifted Wald likelihood. Once the losing accumulator's total finishing
  ## probability exp(-2*k*|v|/s^2) is below one ulp, its survival factor is
  ## exactly 1 and dcswald/pcswald equal dwald/pwald with A = 0.
  b <- 3; v <- 6.5; t0 <- 0.3  # exp(-2*b*v) is about 1e-17
  rt <- c(0.35, 0.5, 0.8, 1.5, 5)
  expect_equal(dcswald(rt, "upper", b = b, t0 = t0, v = v),
               dwald(rt, A = 0, b = b, t0 = t0, v = v), tolerance = 1e-15)
  expect_equal(pcswald(rt, "upper", b = b, t0 = t0, v = v),
               pwald(rt, A = 0, b = b, t0 = t0, v = v), tolerance = 1e-15)
  expect_equal(pcswald(Inf, "upper", b = b, t0 = t0, v = v), 1,
               tolerance = 1e-15)

  ## the same identity holds through the s scaling
  expect_equal(
    dcswald(rt, "upper", b = b/10, t0 = t0, v = v/10, s = 0.1),
    dwald(rt, A = 0, b = b/10, t0 = t0, v = v/10, s = 0.1), tolerance = 1e-15)
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
  expect_equal(pcswald(Inf, "lower", b = 1, t0 = 0.2, v = 0), 0.5,
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
  ## (KS statistic below its 1% critical value, 1.63/sqrt(n))
  for (resp in c("upper", "lower")) {
    mass <- pcswald(Inf, resp, b = b, t0 = t0, v = v)
    rts <- r$rt[r$response == resp]
    ks <- suppressWarnings(ks.test(
      rts, function(q) pcswald(q, resp, b = b, t0 = t0, v = v)/mass))
    expect_true(ks$statistic < 1.63/sqrt(length(rts)))
  }

  ## biased start point: unequal thresholds for the two accumulators
  set.seed(6)
  w <- 0.35
  r_w <- rcswald(4e4, b = b, t0 = t0, v = v, w = w)
  expect_equal(mean(r_w$response == "upper"),
               pcswald(Inf, "upper", b = b, t0 = t0, v = v, w = w),
               tolerance = 0.01)
  ## (subsampled: every pcswald value is a quadrature here)
  rts <- r_w$rt[r_w$response == "lower"]
  rts <- rts[seq_len(min(length(rts), 1000))]
  mass <- pcswald(Inf, "lower", b = b, t0 = t0, v = v, w = w)
  ks <- suppressWarnings(ks.test(
    rts, function(q) pcswald(q, "lower", b = b, t0 = t0, v = v, w = w)/mass))
  expect_true(ks$statistic < 1.63/sqrt(length(rts)))

  ## trialwise parameters are recycled to n
  set.seed(7)
  r_vec <- rcswald(2e4, b = b, t0 = t0, v = c(v, -v))
  expect_equal(mean(r_vec$response[c(TRUE, FALSE)] == "upper"), p_up,
               tolerance = 0.02)
  expect_equal(mean(r_vec$response[c(FALSE, TRUE)] == "lower"), p_up,
               tolerance = 0.02)
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
  expect_identical(pcswald(r, b = 0.8, t0 = 0.25, v = 1.2),
                   pcswald(r$rt, r$response, b = 0.8, t0 = 0.25, v = 1.2))
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
  ## mixed closed-form and numerical CDFs in one call
  expect_identical(
    pcswald(rep(0.9, 3), "upper", b = 0.9, t0 = 0.2, v = 1.1,
            w = c(0.5, 0.4, 0.5)),
    vapply(c(0.5, 0.4, 0.5), function(w)
      pcswald(0.9, "upper", b = 0.9, t0 = 0.2, v = 1.1, w = w), numeric(1)))
})

test_that("dcswald and pcswald are 0 below t0 and dcswald is 0 at infinity", {
  expect_identical(dcswald(c(0.1, 0.2), "upper", b = 1, t0 = 0.2, v = 2),
                   c(0, 0))
  expect_identical(pcswald(c(0.1, 0.2), "upper", b = 1, t0 = 0.2, v = 2),
                   c(0, 0))
  expect_identical(pcswald(c(0.1, 0.2), "upper", b = 1, t0 = 0.2, v = 2,
                           w = 0.4), c(0, 0))
  expect_identical(dcswald(Inf, "upper", b = 1, t0 = 0.2, v = 2), 0)
  ## the survival function at t0 is the response probability
  for (w in c(0.5, 0.4)) {
    expect_equal(
      pcswald(0.2, "upper", b = 1, t0 = 0.2, v = 2, w = w, lower.tail = FALSE),
      pcswald(Inf, "upper", b = 1, t0 = 0.2, v = 2, w = w))
    expect_identical(
      pcswald(Inf, "upper", b = 1, t0 = 0.2, v = 2, w = w, lower.tail = FALSE),
      0)
  }
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

## ---------------------------------------------------------------------------
## quantile function
## ---------------------------------------------------------------------------

test_that("qcswald inverts pcswald", {
  b <- 0.8; v <- 1.2; t0 <- 0.25
  ## defective probabilities up to the response probability, for the
  ## closed-form (w = 0.5) and the numerical (w = 0.4) CDF
  for (w in c(0.5, 0.4)) for (resp in c("upper", "lower")) {
    mass <- pcswald(Inf, resp, b = b, t0 = t0, v = v, w = w)
    p <- c(0.1, 0.5, 0.9) * mass
    q <- qcswald(p, resp, b = b, t0 = t0, v = v, w = w)
    expect_equal(pcswald(q, resp, b = b, t0 = t0, v = v, w = w), p,
                 tolerance = 1e-6)
    ## scale_p rescales by the response probability
    expect_equal(qcswald(c(0.1, 0.5, 0.9), resp, b = b, t0 = t0, v = v,
                         w = w, scale_p = TRUE), q, tolerance = 1e-6)
  }
  ## p = 0 gives t0, probabilities above the response probability give NA
  expect_identical(qcswald(0, "upper", b = b, t0 = t0, v = v), t0)
  expect_warning(
    out <- qcswald(0.5, "lower", b = b, t0 = t0, v = v),
    "predicted response probability")
  expect_identical(out, NA_real_)
  ## the b and a parameterizations agree, and data.frame input works
  expect_identical(
    qcswald(c(0.2, 0.6), "upper", b = 0.8, t0 = 0.2, v = 1.2),
    qcswald(c(0.2, 0.6), "upper", a = 1.6, t0 = 0.2, v = 1.2))
  expect_identical(
    qcswald(data.frame(p = c(0.2, 0.05), response = c("upper", "lower")),
            b = 0.8, t0 = 0.2, v = 1.2),
    qcswald(c(0.2, 0.05), c("upper", "lower"), b = 0.8, t0 = 0.2, v = 1.2))
})
