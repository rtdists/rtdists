context("Shifted Wald distribution functions work correctly")

## The shifted Wald is the A = 0 case of the Wald accumulator. Its four
## distribution functions are checked against (a) the closed-form equations in
## the literature, (b) statmod as an independent reference implementation, and
## (c) each other (integration, inversion, simulation).

test_that("dwald and pwald agree with the shifted Wald equations", {
  rt <- seq(0.05, 3, length.out = 50)
  b <- 1
  v <- 2

  ## Anders, Alario, & Van Maanen (2016), Eq. 3 (= Tillman et al., 2020, Eq. 2)
  eq <- b*(2*pi*rt^3)^(-1/2) * exp(-(1/(2*rt))*(v*rt - b)^2)
  expect_equal(dwald(rt, b = b, t0 = 0, v = v), eq)
  ## A = 0 is the default
  expect_identical(dwald(rt, b = b, t0 = 0, v = v),
                   dwald(rt, A = 0, b = b, t0 = 0, v = v))

  ## shifted by t0
  expect_equal(dwald(rt + 0.3, b = b, t0 = 0.3, v = v), eq)

  ## Tillman et al. (2020), Appendix A, Eq. 11 with a zero start point
  cdf <- pnorm((v*rt - b)/sqrt(rt)) + exp(2*v*b)*pnorm((-v*rt - b)/sqrt(rt))
  expect_equal(pwald(rt, b = b, t0 = 0, v = v), cdf)

  ## Anders, Alario, & Van Maanen (2016), Eq. 1, in their notation
  ## (alpha = b, gamma = v, theta = t0)
  alpha <- b; gamma <- v; theta <- 0.3
  X <- rt + theta
  eq1 <- alpha/sqrt(2*pi*(X - theta)^3) *
    exp(-(alpha - gamma*(X - theta))^2/(2*(X - theta)))
  expect_equal(dwald(X, b = alpha, t0 = theta, v = gamma), eq1)

  ## Giner & Smyth (2016), Eq. 5: the density for an infinite mean, which is
  ## the zero drift rate case, with shape lambda = b^2
  eq5 <- (2*pi*rt^3/b^2)^(-1/2) * exp(-b^2/(2*rt))
  expect_equal(dwald(rt, b = b, t0 = 0, v = 0), eq5)
})

test_that("the moments match Anders et al. (2016)", {
  skip_on_cran()
  ## Eq. 3: E(X) = alpha/gamma + theta;  Eq. 4: SD(X) = sqrt(alpha/gamma^3)
  alpha <- 1; gamma <- 2; theta <- 0.3
  set.seed(9)
  x <- rwald(4e5, b = alpha, t0 = theta, v = gamma)
  expect_equal(mean(x), alpha/gamma + theta, tolerance = 0.01)
  expect_equal(sd(x), sqrt(alpha/gamma^3), tolerance = 0.01)

  ## Eq. 2 of both Anders et al. (2016) and Giner & Smyth (2016): the mode,
  ## with kappa = 3/(2*alpha*gamma), which is the starting value qwald uses
  kappa <- 3/(2*alpha*gamma)
  mode <- (alpha/gamma)*(sqrt(1 + kappa^2) - kappa)
  numerical <- optimize(function(z) dwald(z, b = alpha, t0 = theta, v = gamma),
                        c(theta, theta + 3), maximum = TRUE)$maximum - theta
  expect_equal(numerical, mode, tolerance = 1e-4)
})

test_that("the chi-square identity of Giner & Smyth (2016) holds in both tails", {
  ## Their Eq. 4 gives (X - mu)^2/(phi*X*mu^2) ~ chisq(1), so by their Eq. 6 the
  ## upper chi-square tail equals the sum of the two Wald tails at the roots
  ## q1 < mu < q2 of that quadratic. This checks both tails of pwald at once,
  ## against a function computed by an entirely different route.
  for (pars in list(c(1, 2), c(3, 4), c(2, 0.5))) {
    b <- pars[1]; v <- pars[2]
    mu <- b/v; phi <- 1/b^2
    for (q1 in c(0.1, 0.01)*mu) {
      z <- (q1 - mu)^2/(phi*mu^2*q1)
      q2 <- max(Re(polyroot(c(mu^2, -2*mu - phi*mu^2*z, 1))))
      expect_equal(pwald(q1, b = b, t0 = 0, v = v) +
                     pwald(q2, b = b, t0 = 0, v = v, lower.tail = FALSE),
                   pchisq(z, df = 1, lower.tail = FALSE))
    }
  }
})

test_that("the diffusion constant s only rescales b and v", {
  rt <- seq(0.3, 3, length.out = 20)
  for (s in c(0.5, 0.1)) {
    expect_equal(dwald(rt, b = 1*s, t0 = 0.2, v = 2*s, s = s),
                 dwald(rt, b = 1, t0 = 0.2, v = 2))
    expect_equal(pwald(rt, b = 1*s, t0 = 0.2, v = 2*s, s = s),
                 pwald(rt, b = 1, t0 = 0.2, v = 2))
    expect_equal(qwald(c(0.1, 0.5, 0.9), b = 1*s, t0 = 0.2, v = 2*s, s = s),
                 qwald(c(0.1, 0.5, 0.9), b = 1, t0 = 0.2, v = 2))
  }
})

test_that("d/p/q/rwald agree with statmod", {
  skip_if_not_installed("statmod")
  rt <- c(0.201, 0.25, 0.3, 0.5, 0.9, 1.5, 3, 10, 100)
  t0 <- 0.2
  p <- c(1e-12, 1e-6, 0.01, 0.1, 0.5, 0.9, 0.99, 1 - 1e-6)
  grid <- expand.grid(b = c(0.5, 1, 2, 3), v = c(0.01, 0.5, 2, 4),
                      s = c(1, 0.5, 0.1))
  for (i in seq_len(nrow(grid))) {
    b <- grid$b[i]; v <- grid$v[i]; s <- grid$s[i]
    mu <- b/v
    lambda <- (b/s)^2
    x <- rt - t0
    expect_equal(dwald(rt, b = b, t0 = t0, v = v, s = s),
                 statmod::dinvgauss(x, mean = mu, shape = lambda))
    expect_equal(dwald(rt, b = b, t0 = t0, v = v, s = s, log = TRUE),
                 statmod::dinvgauss(x, mean = mu, shape = lambda, log = TRUE))
    expect_equal(pwald(rt, b = b, t0 = t0, v = v, s = s),
                 statmod::pinvgauss(x, mean = mu, shape = lambda))
    expect_equal(pwald(rt, b = b, t0 = t0, v = v, s = s, lower.tail = FALSE),
                 statmod::pinvgauss(x, mean = mu, shape = lambda,
                                    lower.tail = FALSE))
    expect_equal(pwald(rt, b = b, t0 = t0, v = v, s = s, log.p = TRUE),
                 statmod::pinvgauss(x, mean = mu, shape = lambda,
                                    log.p = TRUE))
    ## statmod::qinvgauss() errors for some of these parameter values (its
    ## left-tail starting value overshoots when the CV is small), so only
    ## compare where it returns.
    ref <- suppressWarnings(
      tryCatch(statmod::qinvgauss(p, mean = mu, shape = lambda),
               error = function(e) NULL))
    if (!is.null(ref))
      expect_equal(qwald(p, b = b, t0 = t0, v = v, s = s) - t0, ref)
  }
})

test_that("pwald is accurate when 2*b*v/s^2 is large", {
  skip_if_not_installed("statmod")
  ## the previous natural-scale implementation lost up to 28% relative accuracy
  ## in this regime, which is reached with the common convention s = 0.1
  for (g in list(c(2, 3, 0.1, 0.5), c(2, 3, 0.1, 0.9),
                 c(3, 4, 0.1, 0.9), c(5, 5, 0.1, 1.2))) {
    b <- g[1]; v <- g[2]; s <- g[3]; rt <- g[4]
    expect_equal(pwald(rt, b = b, t0 = 0.2, v = v, s = s),
                 statmod::pinvgauss(rt - 0.2, mean = b/v, shape = (b/s)^2),
                 tolerance = 1e-10)
  }
})

test_that("dwald integrates to one and pwald is its integral", {
  for (st0 in c(0, 0.15)) {
    expect_equal(integrate(dwald, lower = 0, upper = Inf,
                           b = 1, t0 = 0.2, v = 2, st0 = st0)$value, 1,
                 tolerance = 1e-4)
    for (rt in c(0.4, 0.8, 2)) {
      expect_equal(integrate(dwald, lower = 0, upper = rt,
                             b = 1, t0 = 0.2, v = 2, st0 = st0)$value,
                   pwald(rt, b = 1, t0 = 0.2, v = 2, st0 = st0),
                   tolerance = 1e-5)
    }
  }
})

test_that("qwald inverts pwald", {
  p <- c(1e-12, 1e-6, seq(0.01, 0.99, length.out = 40), 1 - 1e-6)
  grid <- expand.grid(b = c(0.5, 1, 3), v = c(0.5, 2, 4), s = c(1, 0.1))
  for (i in seq_len(nrow(grid))) {
    b <- grid$b[i]; v <- grid$v[i]; s <- grid$s[i]
    q <- qwald(p, b = b, t0 = 0.2, v = v, s = s)
    expect_equal(pwald(q, b = b, t0 = 0.2, v = v, s = s), p,
                 tolerance = 1e-10)
  }
})

test_that("qwald honours lower.tail and log.p", {
  p <- c(1e-9, 0.05, 0.3, 0.5, 0.8, 0.999)
  qu <- qwald(p, b = 1, t0 = 0.2, v = 2, lower.tail = FALSE)
  expect_equal(pwald(qu, b = 1, t0 = 0.2, v = 2, lower.tail = FALSE), p,
               tolerance = 1e-10)
  ql <- qwald(log(p), b = 1, t0 = 0.2, v = 2, log.p = TRUE)
  expect_equal(pwald(ql, b = 1, t0 = 0.2, v = 2), p, tolerance = 1e-10)
  ## boundaries
  expect_equal(qwald(c(0, 1), b = 1, t0 = 0.2, v = 2), c(0.2, Inf))
  expect_true(all(is.na(qwald(c(-0.1, 1.1), b = 1, t0 = 0.2, v = 2))))
})

test_that("log and lower.tail do not lose the tails", {
  ## on the natural scale these underflow to 0 and 1, respectively
  expect_true(is.finite(dwald(0.21, b = 5, t0 = 0.2, v = 2, log = TRUE)))
  expect_true(dwald(0.21, b = 5, t0 = 0.2, v = 2) == 0)
  expect_true(pwald(50, b = 1, t0 = 0.2, v = 4, lower.tail = FALSE) > 0)
  expect_equal(exp(pwald(1.5, b = 1, t0 = 0.2, v = 2, log.p = TRUE)),
               pwald(1.5, b = 1, t0 = 0.2, v = 2))
})

test_that("st0 gives the convolution with a uniform non-decision time", {
  b <- 1; v <- 2; t0 <- 0.2; st0 <- 0.1
  rt <- t0 + c(0.02, 0.05, 0.12, 0.3, 0.8, 2)
  ## density: (F(rt - t0) - F(rt - t0 - st0)) / st0
  ref_d <- (pwald(rt, b = b, t0 = t0, v = v) -
              pwald(rt, b = b, t0 = t0 + st0, v = v))/st0
  expect_equal(dwald(rt, b = b, t0 = t0, v = v, st0 = st0), ref_d)
  ## CDF against numerical integration of that density
  ref_p <- vapply(rt, function(t)
    integrate(dwald, lower = 0, upper = t, b = b, t0 = t0, v = v,
              st0 = st0)$value, 0)
  expect_equal(pwald(rt, b = b, t0 = t0, v = v, st0 = st0), ref_p,
               tolerance = 1e-6)
  ## upper tail and quantile function stay consistent
  expect_equal(pwald(rt, b = b, t0 = t0, v = v, st0 = st0, lower.tail = FALSE),
               1 - pwald(rt, b = b, t0 = t0, v = v, st0 = st0))
  p <- c(0.1, 0.3, 0.5, 0.9)
  expect_equal(pwald(qwald(p, b = b, t0 = t0, v = v, st0 = st0),
                     b = b, t0 = t0, v = v, st0 = st0), p, tolerance = 1e-6)
  ## st0 = 0 must not change the default path
  expect_identical(dwald(rt, b = b, t0 = t0, v = v, st0 = 0),
                   dwald(rt, b = b, t0 = t0, v = v))
  expect_identical(pwald(rt, b = b, t0 = t0, v = v, st0 = 0),
                   pwald(rt, b = b, t0 = t0, v = v))
})

test_that("st0 is rejected for A > 0", {
  expect_error(dwald(0.5, A = 0.2, b = 1, t0 = 0.2, v = 2, st0 = 0.1),
               "only supported for A = 0")
  expect_error(pwald(0.5, A = 0.2, b = 1, t0 = 0.2, v = 2, st0 = 0.1),
               "only supported for A = 0")
  expect_error(qwald(0.5, A = 0.2, b = 1, t0 = 0.2, v = 2, st0 = 0.1),
               "only supported for A = 0")
})

test_that("a drift rate of zero gives the Levy distribution", {
  rt <- c(0.3, 0.5, 1, 5)
  b <- 1; t0 <- 0.2
  x <- rt - t0
  ## Levy: F(x) = 2 * Phi(-sqrt(lambda/x)) with lambda = b^2
  expect_equal(pwald(rt, b = b, t0 = t0, v = 0), 2*pnorm(-sqrt(b^2/x)))
  expect_equal(dwald(rt, b = b, t0 = t0, v = 0),
               b/sqrt(2*pi*x^3)*exp(-b^2/(2*x)))
  ## quantile: lambda / qchisq(1 - p, 1)
  p <- c(0.1, 0.5, 0.9)
  expect_equal(qwald(p, b = b, t0 = t0, v = 0),
               t0 + b^2/qchisq(p, df = 1, lower.tail = FALSE))
})

test_that("dwald and pwald are 0 below t0 and for negative drift rates", {
  expect_equal(dwald(c(0.1, 0.2), b = 1, t0 = 0.2, v = 2), c(0, 0))
  expect_equal(pwald(c(0.1, 0.2), b = 1, t0 = 0.2, v = 2), c(0, 0))
  expect_equal(dwald(c(0.5, 1), b = 1, t0 = 0.2, v = -2), c(0, 0))
  expect_equal(pwald(c(0.5, 1), b = 1, t0 = 0.2, v = -2), c(0, 0))
  ## the upper tail of a process that has not (or will never) finish is 1
  expect_equal(pwald(c(0.1, 0.2), b = 1, t0 = 0.2, v = 2, lower.tail = FALSE),
               c(1, 1))
  ## NA/NaN parameters give 0 in the cores, which is what the race relies on
  ## (the exported functions error on them, as the LBA ones do)
  for (par in list(list(b = NA_real_, v = 2), list(b = NaN, v = 2),
                   list(b = 1, v = NA_real_), list(b = 1, v = NaN))) {
    expect_equal(dwald_core(0.5, A = 0, b = par$b, t0 = 0.2, v = par$v,
                            s = 1, nn = 1), 0)
    expect_equal(pwald_core(0.5, A = 0, b = par$b, t0 = 0.2, v = par$v,
                            s = 1, nn = 1), 0)
  }
})

test_that("the distribution is proper, i.e. pwald(Inf) is 1", {
  for (A in c(0, 0.3)) {
    for (v in c(2, 0)) {
      expect_equal(pwald(Inf, A = A, b = 1, t0 = 0.2, v = v), 1)
      expect_equal(pwald(Inf, A = A, b = 1, t0 = 0.2, v = v,
                         lower.tail = FALSE), 0)
    }
    ## an accumulator with a negative drift rate never finishes
    expect_equal(pwald(Inf, A = A, b = 1, t0 = 0.2, v = -2), 0)
  }
  expect_equal(pwald(Inf, b = 1, t0 = 0.2, v = 2, st0 = 0.15), 1)
  expect_equal(dwald(Inf, b = 1, t0 = 0.2, v = 2), 0)
  ## Inf mixed into a vector does not affect the other elements
  expect_equal(pwald(c(0.5, Inf, 1), b = 1, t0 = 0.2, v = 2),
               c(pwald(0.5, b = 1, t0 = 0.2, v = 2), 1,
                 pwald(1, b = 1, t0 = 0.2, v = 2)))
})

test_that("all parameters are recycled and rt order does not matter", {
  rt <- c(0.5, 0.9, 1.4, 2.1)
  vs <- c(1, 2, 3, 4)
  expect_equal(dwald(rt, b = 1, t0 = 0.2, v = vs),
               vapply(seq_along(rt), function(i)
                 dwald(rt[i], b = 1, t0 = 0.2, v = vs[i]), 0))
  expect_equal(pwald(rt, b = 1, t0 = 0.2, v = vs),
               vapply(seq_along(rt), function(i)
                 pwald(rt[i], b = 1, t0 = 0.2, v = vs[i]), 0))
  expect_equal(dwald(rt, b = c(0.8, 1.2), t0 = c(0.1, 0.2), v = 2),
               dwald(rt, b = c(0.8, 1.2, 0.8, 1.2), t0 = c(0.1, 0.2, 0.1, 0.2),
                     v = c(2, 2, 2, 2)))
  ## order invariance
  o <- c(3, 1, 4, 2)
  expect_equal(pwald(rt[o], b = 1, t0 = 0.2, v = vs[o]),
               pwald(rt, b = 1, t0 = 0.2, v = vs)[o])
})

test_that("input is checked", {
  expect_error(dwald(0.5, A = 2, b = 1, t0 = 0.2, v = 2),
               "b cannot be smaller than A")
  expect_error(pwald(0.5, A = 2, b = 1, t0 = 0.2, v = 2),
               "b cannot be smaller than A")
  expect_error(dwald(-0.5, b = 1, t0 = 0.2, v = 2),
               "rt needs to contain only positive values")
  expect_error(dwald(0.5, b = "1", t0 = 0.2, v = 2),
               "b needs to be a numeric vector")
  expect_error(rwald(c(10, 10), b = 1, t0 = 0.2, v = 2),
               "n needs to be of length 1")
})

test_that("rwald returns RTs from the right distribution", {
  skip_on_cran()
  set.seed(4)
  n <- 1e5
  x <- rwald(n, b = 1, t0 = 0.2, v = 2)
  expect_true(is.numeric(x))
  expect_false(is.matrix(x))
  expect_equal(length(x), n)
  ## mean of a shifted Wald is t0 + b/v
  expect_equal(mean(x), 0.2 + 1/2, tolerance = 0.01)
  p <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  expect_equal(unname(quantile(x, p)), qwald(p, b = 1, t0 = 0.2, v = 2),
               tolerance = 0.01)
  ## with variability in non-decision time
  set.seed(5)
  y <- rwald(n, b = 1, t0 = 0.2, v = 2, st0 = 0.15)
  expect_equal(mean(y), 0.2 + 0.15/2 + 1/2, tolerance = 0.01)
  expect_equal(unname(quantile(y, p)),
               qwald(p, b = 1, t0 = 0.2, v = 2, st0 = 0.15), tolerance = 0.01)
  ## a race of two accumulators still returns rt and response
  set.seed(6)
  z <- rwald(100, A = 0.2, b = 1, t0 = 0.2, v = c(2, 1))
  expect_true(is.matrix(z))
  expect_equal(colnames(z), c("rt", "response"))
})

test_that("the Wald accumulator still works as part of a race", {
  ## the shifted Wald branch is what dRDM/n1PDF call with A = 0
  rt <- c(0.3, 0.5, 0.9, 2)
  expect_equal(
    dRDM(rt, response = 1, A = 0, b = 1, t0 = 0.2, v = c(2, 1), silent = TRUE),
    dwald(rt, b = 1, t0 = 0.2, v = 2) *
      (1 - pwald(rt, b = 1, t0 = 0.2, v = 1)))
  expect_equal(
    pRDM(Inf, response = 1, A = 0, b = 1, t0 = 0.2, v = c(2, 1),
         silent = TRUE) +
      pRDM(Inf, response = 2, A = 0, b = 1, t0 = 0.2, v = c(2, 1),
           silent = TRUE), 1, tolerance = 1e-4)
})
