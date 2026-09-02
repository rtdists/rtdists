
context("RDM works correctly")

test_that("dwald agrees with the published density equations", {
  rt <- seq(0.05, 3, length.out = 50)
  b <- 1
  v <- 2
  
  # Tillman et al. (2020), Equation 2: shifted Wald, no start point variability
  eq2 <- b*(2*pi*rt^3)^(-1/2) * exp(-(1/(2*rt))*(v*rt - b)^2)
  expect_equal(dwald(rt, A = 0, b = b, t0 = 0, v = v), eq2)
  
  # Equation 5: with start point variability
  A <- 0.4
  alpha <- (b - A - rt*v)/sqrt(rt)
  beta <- (b - rt*v)/sqrt(rt)
  eq5 <- (1/A)*(-v*pnorm(alpha) + dnorm(alpha)/sqrt(rt) + 
                  v*pnorm(beta) - dnorm(beta)/sqrt(rt))
  expect_equal(dwald(rt, A = A, b = b, t0 = 0, v = v), eq5)
  
  # Equation 6: with start point variability and a drift rate of zero
  eq6 <- (1/A)*(dnorm((b - A)/sqrt(rt))/sqrt(rt) - dnorm(b/sqrt(rt))/sqrt(rt))
  expect_equal(dwald(rt, A = A, b = b, t0 = 0, v = 0), eq6)
})

test_that("dwald integrates to one and pwald is its integral", {
  for (A in c(0, 0.4, 1)) {
    expect_equal(integrate(dwald, lower = 0, upper = Inf, 
                           A = A, b = 1, t0 = 0.2, v = 2)$value, 1, 
                 tolerance = 1e-4)
    for (rt in c(0.4, 0.8, 2)) {
      expect_equal(integrate(dwald, lower = 0, upper = rt, 
                             A = A, b = 1, t0 = 0.2, v = 2)$value,
                   pwald(rt, A = A, b = 1, t0 = 0.2, v = 2), 
                   tolerance = 1e-5)
    }
  }
})

test_that("pwald keeps its precision when 2*b*v/s^2 is large", {
  # the terms inside exp() underflowed, e.g. with s = 0.1 (Ratcliff convention)
  for (s in c(0.05, 0.1)) {
    for (A in c(0, 0.2)) {
      for (rt in c(0.7, 0.85, 1, 1.5, 3)) {
        expect_equal(integrate(dwald, lower = 0, upper = rt, 
                               A = A, b = 2, t0 = 0, v = 2, s = s, 
                               rel.tol = 1e-10, subdivisions = 1000L)$value,
                     pwald(rt, A = A, b = 2, t0 = 0, v = 2, s = s), 
                     tolerance = 1e-8)
      }
    }
  }
})

test_that("pwald is 1 at rt = Inf", {
  expect_equal(pwald(Inf, A = 0, b = 1, t0 = 0.2, v = 2), 1)
  expect_equal(pwald(Inf, A = 0.5, b = 1, t0 = 0.2, v = 2), 1)
  expect_equal(pwald(c(0.5, Inf, 1), A = 0.5, b = 1, t0 = 0.2, v = 2)[2], 1)
})

test_that("dwald and pwald are 0 below t0 and for negative drift rates", {
  expect_equal(dwald(c(0.1, 0.3), A = 0.4, b = 1, t0 = 0.5, v = 2), c(0, 0))
  expect_equal(pwald(c(0.1, 0.3), A = 0.4, b = 1, t0 = 0.5, v = 2), c(0, 0))
  expect_equal(dwald(1, A = 0.4, b = 1, t0 = 0.2, v = -1), 0)
  expect_equal(pwald(1, A = 0.4, b = 1, t0 = 0.2, v = -1), 0)
})

test_that("accumulators with a negative drift rate never win", {
  set.seed(3)
  x <- rRDM(200, A = 0.4, b = 1, t0 = 0.25, v = c(2.5, -1), silent = TRUE)
  expect_true(all(x$response == 1))
  expect_equal(dRDM(x$rt, rep(2, nrow(x)), A = 0.4, b = 1, t0 = 0.25, 
                    v = c(2.5, -1), silent = TRUE), rep(0, nrow(x)))
})

test_that("missing parameters give 0 as they do for the LBA", {
  for (par in list(list(A = NA_real_, b = 1), list(A = NaN, b = 1), 
                   list(A = 0.4, b = NA_real_), list(A = 0.4, b = NaN))) {
    expect_equal(n1PDF(1, A = par$A, b = par$b, t0 = 0.2, v = c(2, 1), 
                       distribution = "wald", silent = TRUE), 
                 n1PDF(1, A = par$A, b = par$b, t0 = 0.2, mean_v = c(2, 1), 
                       sd_v = c(1, 1), silent = TRUE))
  }
})

test_that("dRDM with wald works as expected", {
  rt <- 2
  vs <- c(0.8, 0.4)
  
  f <- dwald(rt, A = 0.5, b = 1, t0 = 0.5, v = vs[1])
  F <- pwald(rt, A = 0.5, b = 1, t0 = 0.5, v = vs[2])
  
  expect_equivalent(n1PDF(rt, A = 0.5, b = 1, t0 = 0.5, v = vs, 
                          distribution = "wald", silent = TRUE), f*(1-F))
  expect_equivalent(dRDM(rt, 1, A = 0.5, b = 1, t0 = 0.5, v = vs, 
                         silent = TRUE), f*(1-F))
})

test_that("dRDM is identical to n1PDF", {
  n <- 100
  x <- rRDM(n, A = 0.5, b = 1, t0 = 0.5, v = c(1.2, 1), silent = TRUE)
  
  n1 <- vector("numeric", n)
  n1[x$response == 1] <- n1PDF(x$rt[x$response == 1], A = 0.5, b = 1, t0 = 0.5, 
                               v = c(1.2, 1), distribution = "wald", 
                               silent = TRUE)
  n1[x$response == 2] <- n1PDF(x$rt[x$response == 2], A = 0.5, b = 1, t0 = 0.5, 
                               v = c(1, 1.2), distribution = "wald", 
                               silent = TRUE)
  expect_identical(n1, dRDM(x$rt, x$response, A = 0.5, b = 1, t0 = 0.5, 
                            v = c(1.2, 1), silent = TRUE))
})

test_that("dRDM is the defective density of a race with three accumulators", {
  rt <- seq(0.4, 2.5, length.out = 20)
  v <- c(2.5, 1.5, 1)
  A <- 0.4
  b <- 1
  t0 <- 0.25
  
  for (i in 1:3) {
    others <- setdiff(1:3, i)
    manual <- dwald(rt, A = A, b = b, t0 = t0, v = v[i])
    for (j in others) 
      manual <- manual * (1 - pwald(rt, A = A, b = b, t0 = t0, v = v[j]))
    expect_equal(dRDM(rt, i, A = A, b = b, t0 = t0, v = v, silent = TRUE), 
                 manual)
  }
})

test_that("dRDM sums to one across responses and pRDM is its integral", {
  A <- 0.4; b <- 1; t0 <- 0.25; v <- c(2.5, 1.5)
  mass <- 0
  for (i in 1:2) 
    mass <- mass + integrate(dRDM, lower = 0, upper = Inf, response = i, 
                             A = A, b = b, t0 = t0, v = v, silent = TRUE)$value
  expect_equal(mass, 1, tolerance = 1e-4)
  
  for (i in 1:2) for (rt in c(0.5, 1, 2)) 
    expect_equal(integrate(dRDM, lower = 0, upper = rt, response = i, 
                           A = A, b = b, t0 = t0, v = v, silent = TRUE)$value,
                 pRDM(rt, i, A = A, b = b, t0 = t0, v = v, silent = TRUE),
                 tolerance = 1e-5)
})

test_that("qRDM inverts pRDM", {
  A <- 0.4; b <- 1; t0 <- 0.25; v <- c(2.5, 1.5)
  p <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  q <- qRDM(p, 1, A = A, b = b, t0 = t0, v = v, scale_p = TRUE, silent = TRUE)
  max_p <- pRDM(Inf, 1, A = A, b = b, t0 = t0, v = v, silent = TRUE)
  expect_equal(pRDM(q, 1, A = A, b = b, t0 = t0, v = v, silent = TRUE)/max_p, p,
               tolerance = 1e-4)
})

test_that("s only rescales A, b, and v", {
  rt <- seq(0.4, 2.5, length.out = 20)
  expect_equal(dRDM(rt, 1, A = 0.4, b = 1, t0 = 0.25, v = c(2.5, 1.5), 
                    silent = TRUE),
               dRDM(rt, 1, A = 0.8, b = 2, t0 = 0.25, v = c(5, 3), s = 2, 
                    silent = TRUE))
  expect_equal(dwald(rt, A = 0.4, b = 1, t0 = 0.25, v = 2), 
               dwald(rt, A = 0.04, b = 0.1, t0 = 0.25, v = 0.2, s = 0.1))
})

test_that("RDM parameters can vary across trials", {
  set.seed(1)
  n <- 20
  rt <- runif(n, 0.4, 2)
  response <- sample(1:2, n, replace = TRUE)
  A <- runif(n, 0.2, 0.5)
  b <- A + runif(n, 0.4, 0.8)
  t0 <- runif(n, 0.15, 0.3)
  v1 <- runif(n, 1.5, 3)
  v2 <- runif(n, 0.8, 2)
  
  single <- vapply(seq_len(n), function(i) 
    dRDM(rt[i], response[i], A = A[i], b = b[i], t0 = t0[i], 
         v = c(v1[i], v2[i]), silent = TRUE), 0)
  expect_equal(dRDM(rt, response, A = list(A, A), b = list(b, b), 
                    t0 = list(t0, t0), v = list(v1, v2), silent = TRUE), single)
  
  single <- vapply(seq_len(n), function(i) 
    pRDM(rt[i], response[i], A = A[i], b = b[i], t0 = t0[i], 
         v = c(v1[i], v2[i]), silent = TRUE), 0)
  expect_equal(pRDM(rt, response, A = list(A, A), b = list(b, b), 
                    t0 = list(t0, t0), v = list(v1, v2), silent = TRUE), single)
})

test_that("dRDM and pRDM do not depend on the order of rt", {
  set.seed(1)
  rt <- runif(30, 0.35, 2)
  response <- sample(1:2, 30, replace = TRUE)
  o <- order(rt)
  expect_equal(dRDM(rt, response, A = 0.4, b = 1, t0 = 0.25, v = c(2.5, 1.5), 
                    silent = TRUE)[o],
               dRDM(rt[o], response[o], A = 0.4, b = 1, t0 = 0.25, 
                    v = c(2.5, 1.5), silent = TRUE))
  expect_equal(pRDM(rt, response, A = 0.4, b = 1, t0 = 0.25, v = c(2.5, 1.5), 
                    silent = TRUE)[o],
               pRDM(rt[o], response[o], A = 0.4, b = 1, t0 = 0.25, 
                    v = c(2.5, 1.5), silent = TRUE))
})

test_that("qRDM uses the t0 of each trial", {
  p <- c(0.2, 0.5, 0.8)
  t0 <- c(0.2, 0.9, 0.2)
  vectorised <- qRDM(p, 1, A = 0.4, b = 1, t0 = list(t0, t0), v = c(2.5, 1.5), 
                     scale_p = TRUE, silent = TRUE)
  single <- vapply(seq_along(p), function(i) 
    qRDM(p[i], 1, A = 0.4, b = 1, t0 = t0[i], v = c(2.5, 1.5), 
         scale_p = TRUE, silent = TRUE), 0)
  expect_false(anyNA(vectorised))
  expect_equal(vectorised, single)
})

test_that("dRDM and pRDM accept a data.frame", {
  df <- data.frame(rt = seq(0.4, 2.5, length.out = 10), response = rep(1:2, 5))
  expect_equal(dRDM(df, A = 0.4, b = 1, t0 = 0.25, v = c(2.5, 1.5), 
                    silent = TRUE),
               dRDM(df$rt, df$response, A = 0.4, b = 1, t0 = 0.25, 
                    v = c(2.5, 1.5), silent = TRUE))
  expect_equal(pRDM(df, A = 0.4, b = 1, t0 = 0.25, v = c(2.5, 1.5), 
                    silent = TRUE),
               pRDM(df$rt, df$response, A = 0.4, b = 1, t0 = 0.25, 
                    v = c(2.5, 1.5), silent = TRUE))
})

test_that("RDM functions check their input", {
  expect_error(dwald(1, A = 1, b = 0.5, t0 = 0.2, v = 2), "b cannot be smaller")
  expect_error(rwald(10, A = 1, b = 0.5, t0 = 0.2, v = c(2, 1)), 
               "b cannot be smaller")
  expect_error(n1PDF(1, A = 0.5, b = 1, t0 = 0.2, mean_v = c(2, 1), 
                     distribution = "wald", silent = TRUE), 
               "v needs to be passed")
  expect_error(dRDM(1, 3, A = 0.5, b = 1, t0 = 0.2, v = c(2, 1), silent = TRUE), 
               "response needs to be")
})

test_that("rRDM recovers the predicted RT distribution", {
  skip_on_cran()
  A <- 0.4; b <- 1; t0 <- 0.25; v <- c(2.5, 1.5)
  set.seed(2)
  x <- rRDM(1e5, A = A, b = b, t0 = t0, v = v, silent = TRUE)
  
  for (i in 1:2) {
    expect_equal(mean(x$response == i), 
                 pRDM(Inf, i, A = A, b = b, t0 = t0, v = v, silent = TRUE), 
                 tolerance = 0.01)
    q <- as.numeric(quantile(x$rt[x$response == i], c(0.1, 0.3, 0.5, 0.7, 0.9)))
    empirical <- vapply(q, function(z) mean(x$rt <= z & x$response == i), 0)
    expect_equal(empirical, 
                 pRDM(q, i, A = A, b = b, t0 = t0, v = v, silent = TRUE), 
                 tolerance = 0.01)
  }
})
