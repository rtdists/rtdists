context("LBA: A_small branch bugs")

## Regression tests for bugs in the A_small (< 1e-10) branches of the
## lognormal and gamma single-accumulator functions. When any element of
## A is below 1e-10, the A_small code path is triggered, which previously
## contained several bugs:
##   1. dlba_lnorm_core: used `max` instead of `min` in zlognorm, Gmin, u, v
##   2. plba_lnorm_core: same max/min issue; also missing rem_t0 call
##   3. dlba_gamma_core: diffG function double-subsetted shape_v; term1-3
##      used full vectors instead of subsetted ones
##   4. plba_gamma_core: term1/term2 used full vectors instead of subsetted

test_that("dlba_lnorm with mixed A values matches individual calls", {
  d_mixed <- dlba_lnorm(
    rt = rep(1, 4), A = c(1e-12, 0.5, 1e-12, 0.5), b = rep(1, 4),
    t0 = rep(0, 4), meanlog_v = rep(0, 4), sdlog_v = rep(0.5, 4)
  )
  d_small  <- dlba_lnorm(rt = 1, A = 1e-12, b = 1, t0 = 0, meanlog_v = 0, sdlog_v = 0.5)
  d_normal <- dlba_lnorm(rt = 1, A = 0.5,   b = 1, t0 = 0, meanlog_v = 0, sdlog_v = 0.5)
  expect_equal(d_mixed, c(d_small, d_normal, d_small, d_normal))
})

test_that("plba_lnorm with mixed A values matches individual calls", {
  p_mixed <- plba_lnorm(
    rt = rep(1, 4), A = c(1e-12, 0.5, 1e-12, 0.5), b = rep(1, 4),
    t0 = rep(0, 4), meanlog_v = rep(0, 4), sdlog_v = rep(0.5, 4)
  )
  p_small  <- plba_lnorm(rt = 1, A = 1e-12, b = 1, t0 = 0, meanlog_v = 0, sdlog_v = 0.5)
  p_normal <- plba_lnorm(rt = 1, A = 0.5,   b = 1, t0 = 0, meanlog_v = 0, sdlog_v = 0.5)
  expect_equal(p_mixed, c(p_small, p_normal, p_small, p_normal))
})

test_that("plba_lnorm with A small and t0 > 0 applies rem_t0", {
  p_with_t0  <- plba_lnorm(rt = 2, A = 1e-12, b = 1, t0 = 0.5, meanlog_v = 0, sdlog_v = 0.5)
  p_expected <- plba_lnorm(rt = 1.5, A = 1e-12, b = 1, t0 = 0,   meanlog_v = 0, sdlog_v = 0.5)
  expect_equal(p_with_t0, p_expected)
})

test_that("plba_lnorm with mixed A and t0 > 0 applies rem_t0 to all entries", {
  p_mixed <- plba_lnorm(
    rt = rep(2, 4), A = c(1e-12, 0.5, 1e-12, 0.5), b = rep(1, 4),
    t0 = rep(0.5, 4), meanlog_v = rep(0, 4), sdlog_v = rep(0.5, 4)
  )
  p_small  <- plba_lnorm(rt = 1.5, A = 1e-12, b = 1, t0 = 0, meanlog_v = 0, sdlog_v = 0.5)
  p_normal <- plba_lnorm(rt = 1.5, A = 0.5,   b = 1, t0 = 0, meanlog_v = 0, sdlog_v = 0.5)
  expect_equal(p_mixed, c(p_small, p_normal, p_small, p_normal))
})

test_that("dlba_gamma with mixed A values matches individual calls", {
  d_mixed <- dlba_gamma(
    rt = rep(1, 4), A = c(1e-12, 0.5, 1e-12, 0.5), b = rep(1, 4),
    t0 = rep(0, 4), shape_v = rep(1, 4), rate_v = rep(1, 4)
  )
  d_small  <- dlba_gamma(rt = 1, A = 1e-12, b = 1, t0 = 0, shape_v = 1, rate_v = 1)
  d_normal <- dlba_gamma(rt = 1, A = 0.5,   b = 1, t0 = 0, shape_v = 1, rate_v = 1)
  expect_equal(d_mixed, c(d_small, d_normal, d_small, d_normal))
})

test_that("plba_gamma with mixed A values matches individual calls", {
  p_mixed <- plba_gamma(
    rt = rep(1, 4), A = c(1e-12, 0.5, 1e-12, 0.5), b = rep(1, 4),
    t0 = rep(0, 4), shape_v = rep(1, 4), rate_v = rep(1, 4)
  )
  p_small  <- plba_gamma(rt = 1, A = 1e-12, b = 1, t0 = 0, shape_v = 1, rate_v = 1)
  p_normal <- plba_gamma(rt = 1, A = 0.5,   b = 1, t0 = 0, shape_v = 1, rate_v = 1)
  expect_equal(p_mixed, c(p_small, p_normal, p_small, p_normal))
})

test_that("dlba_lnorm mixed A produces no warnings", {
  expect_silent(dlba_lnorm(
    rt = rep(1, 4), A = c(1e-12, 0.5, 1e-12, 0.5), b = rep(1, 4),
    t0 = rep(0, 4), meanlog_v = rep(0, 4), sdlog_v = rep(0.5, 4)
  ))
})

test_that("dlba_gamma mixed A produces no warnings", {
  expect_silent(dlba_gamma(
    rt = rep(1, 4), A = c(1e-12, 0.5, 1e-12, 0.5), b = rep(1, 4),
    t0 = rep(0, 4), shape_v = rep(1, 4), rate_v = rep(1, 4)
  ))
})

test_that("plba_lnorm mixed A produces no warnings", {
  expect_silent(plba_lnorm(
    rt = rep(1, 4), A = c(1e-12, 0.5, 1e-12, 0.5), b = rep(1, 4),
    t0 = rep(0, 4), meanlog_v = rep(0, 4), sdlog_v = rep(0.5, 4)
  ))
})

test_that("plba_gamma mixed A produces no warnings", {
  expect_silent(plba_gamma(
    rt = rep(1, 4), A = c(1e-12, 0.5, 1e-12, 0.5), b = rep(1, 4),
    t0 = rep(0, 4), shape_v = rep(1, 4), rate_v = rep(1, 4)
  ))
})

## Frechet functions do not have an A_small branch; they use a ps_below_zero
## filter that subsets all variables upfront. These tests confirm that mixed
## A values and rem_t0 are handled correctly.

test_that("dlba_frechet with mixed A values matches individual calls", {
  d_mixed <- dlba_frechet(
    rt = rep(1, 4), A = c(1e-12, 0.5, 1e-12, 0.5), b = rep(1, 4),
    t0 = rep(0, 4), shape_v = rep(2, 4), scale_v = rep(1, 4)
  )
  d_small  <- dlba_frechet(rt = 1, A = 1e-12, b = 1, t0 = 0, shape_v = 2, scale_v = 1)
  d_normal <- dlba_frechet(rt = 1, A = 0.5,   b = 1, t0 = 0, shape_v = 2, scale_v = 1)
  expect_equal(d_mixed, c(d_small, d_normal, d_small, d_normal))
})

test_that("plba_frechet with mixed A values matches individual calls", {
  p_mixed <- plba_frechet(
    rt = rep(1, 4), A = c(1e-12, 0.5, 1e-12, 0.5), b = rep(1, 4),
    t0 = rep(0, 4), shape_v = rep(2, 4), scale_v = rep(1, 4)
  )
  p_small  <- plba_frechet(rt = 1, A = 1e-12, b = 1, t0 = 0, shape_v = 2, scale_v = 1)
  p_normal <- plba_frechet(rt = 1, A = 0.5,   b = 1, t0 = 0, shape_v = 2, scale_v = 1)
  expect_equal(p_mixed, c(p_small, p_normal, p_small, p_normal))
})

test_that("plba_frechet with A small and t0 > 0 applies rem_t0", {
  p_with_t0  <- plba_frechet(rt = 2, A = 1e-12, b = 1, t0 = 0.5, shape_v = 2, scale_v = 1)
  p_expected <- plba_frechet(rt = 1.5, A = 1e-12, b = 1, t0 = 0,   shape_v = 2, scale_v = 1)
  expect_equal(p_with_t0, p_expected)
})
