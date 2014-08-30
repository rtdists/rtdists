
context("rrd: random number generation for the diffusion model")

test_that("rrd works", {
  expect_that(rrd(10, c(rep(1, 4), rep(0.1, 4))), is_a("matrix"))
})
