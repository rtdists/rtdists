

context("LBA vectorization")

test_that("_norm vectorizes", {
  n <- 10
  x <- rlba_norm(n, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=c(0.2,0.3))
  
  as <- seq(0.2, 0.5, length.out = n)
  bs <- seq(0.5, 1.5, length.out = n)
  t0s <- seq(0, 0.5, length.out = n/2)
  vs <- seq(0.8, 1.2, length.out = 10)
  
  o1 <- plba_norm(x$response, A=as, b=bs, t0 = t0s, mean_v=vs, sd_v=0.2)
  o2 <- mapply(plba_norm, t = x$response, A = as, b = bs, t0 = t0s, mean_v=vs, MoreArgs = list(sd_v=0.2))
  expect_identical(o1, o2)
  p1 <- dlba_norm(x$response, A=as, b=bs, t0 = t0s, mean_v=1.2, sd_v=vs)
  p2 <- mapply(dlba_norm, t = x$response, A = as, b = bs, t0 = t0s, sd_v=vs, MoreArgs = list(mean_v=1.2))
  expect_identical(p1, p2)
})

test_that("_gamma vectorizes", {
  n <- 10
  x <- rlba_norm(n, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=c(0.2,0.3))
  
  as <- seq(0.2, 0.5, length.out = n)
  bs <- seq(0.5, 1.5, length.out = n)
  t0s <- seq(0, 0.5, length.out = n/2)
  vs <- seq(0.8, 1.2, length.out = 10)
  
  o1 <- dlba_gamma(x$response, A=as, b=bs, t0 = t0s, shape_v=vs, scale_v=0.2)
  o2 <- mapply(dlba_gamma, t = x$response, A = as, b = bs, t0 = t0s, shape_v=vs, MoreArgs = list(scale_v=0.2))
  expect_identical(o1, o2)
  p1 <- plba_gamma(x$response, A=as, b=bs, t0 = t0s, shape_v=1.2, scale_v=vs)
  p2 <- mapply(plba_gamma, t = x$response, A = as, b = bs, t0 = t0s, scale_v=vs, MoreArgs = list(shape_v=1.2))
  expect_identical(p1, p2)
})

test_that("_frechet vectorizes", {
  n <- 10
  x <- rlba_norm(n, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=c(0.2,0.3))
  
  as <- seq(0.2, 0.5, length.out = n)
  bs <- seq(0.5, 1.5, length.out = n)
  t0s <- seq(0, 0.5, length.out = n/2)
  vs <- seq(0.8, 1.2, length.out = 10)
  
  o1 <- dlba_frechet(x$response, A=as, b=bs, t0 = t0s, shape_v=vs, scale_v=0.2)
  o2 <- mapply(dlba_frechet, t = x$response, A = as, b = bs, t0 = t0s, shape_v=vs, MoreArgs = list(scale_v=0.2))
  expect_identical(o1, o2)
  p1 <- plba_frechet(x$response, A=as, b=bs, t0 = t0s, shape_v=1.2, scale_v=vs)
  p2 <- mapply(plba_frechet, t = x$response, A = as, b = bs, t0 = t0s, scale_v=vs, MoreArgs = list(shape_v=1.2))
  expect_identical(p1, p2)
})

test_that("_lnorm vectorizes", {
  n <- 10
  x <- rlba_norm(n, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=c(0.2,0.3))
  
  as <- seq(0.2, 0.5, length.out = n)
  bs <- seq(0.5, 1.5, length.out = n)
  t0s <- seq(0, 0.5, length.out = n/2)
  vs <- seq(0.8, 1.2, length.out = 10)
  
  o1 <- dlba_lnorm(x$response, A=as, b=bs, t0 = t0s, meanlog_v=vs, sdlog_v=0.2)
  o2 <- mapply(dlba_lnorm, t = x$response, A = as, b = bs, t0 = t0s, meanlog_v=vs, MoreArgs = list(sdlog_v=0.2))
  expect_identical(o1, o2)
  p1 <- plba_lnorm(x$response, A=as, b=bs, t0 = t0s, meanlog_v=1.2, sdlog_v=vs)
  p2 <- mapply(plba_lnorm, t = x$response, A = as, b = bs, t0 = t0s, sdlog_v=vs, MoreArgs = list(meanlog_v=1.2))
  expect_identical(p1, p2)
})

context("LBA input")

test_that("_norm input works as they should", {
  expect_error(rlba_norm("10", A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=c(0.2,0.3)))
  expect_error(rlba_norm(10, A=c(0.5, 0.6), b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=c(0.2,0.3)))    
  x <- rlba_norm(10, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=c(0.2,0.3))
  expect_error(plba_norm(x$response, A=as.character(seq(0.2, 0.5, length.out = 10)), b=0.5, t0 = 0.5, mean_v=1.2, sd_v=0.2))
  expect_is(plba_norm(x$response, A=seq(0.2, 0.5, length.out = 10), b=0.5, t0 = seq(0.5, 0.7, 0.1), mean_v=1.2, sd_v=0.2), "numeric")
})

test_that("_gamma input works as they should", {
  expect_error(rlba_gamma("10", A=0.5, b=1, t0 = 0.5, shape_v =c(1.2, 1), scale_v =c(0.2,0.3)), regexp = "numeric and finite")
  expect_error(rlba_gamma(10, A=c(0.5, 0.6), b=1, t0 = 0.5, shape_v=c(1.2, 1), scale_v=c(0.2,0.3)), "length") 
  x <- rlba_gamma(10, A=0.5, b=1, t0 = 0.5, shape_v=c(1.2, 1), scale_v=c(0.2,0.3))
  expect_error(plba_gamma(x$response, A=as.character(seq(0.2, 0.5, length.out = 10)), b=0.5, t0 = 0.5, shape_v=1.2, scale_v=0.2), "numeric")
  expect_is(plba_gamma(x$response, A=seq(0.2, 0.5, length.out = 10), b=0.5, t0 = seq(0.5, 0.7, 0.1), shape_v=1.2, scale_v=0.2), "numeric")
})

test_that("_frechet input works as they should", {
  expect_error(rlba_frechet("10", A=0.5, b=1, t0 = 0.5, shape_v =c(1.2, 1), scale_v =c(0.2,0.3)), regexp = "numeric and finite")
  expect_error(rlba_frechet(10, A=c(0.5, 0.6), b=1, t0 = 0.5, shape_v=c(1.2, 1), scale_v=c(0.2,0.3)), "length") 
  x <- rlba_frechet(10, A=0.5, b=1, t0 = 0.5, shape_v=c(1.2, 1), scale_v=c(0.2,0.3))
  expect_error(plba_frechet(x$response, A=as.character(seq(0.2, 0.5, length.out = 10)), b=0.5, t0 = 0.5, shape_v=1.2, scale_v=0.2), "numeric")
  expect_is(plba_frechet(x$response, A=seq(0.2, 0.5, length.out = 10), b=0.5, t0 = seq(0.5, 0.7, 0.1), shape_v=1.2, scale_v=0.2), "numeric")
})

test_that("_lnorm input works as they should", {
  expect_error(rlba_lnorm("10", A=0.5, b=1, t0 = 0.5, meanlog_v =c(1.2, 1), sdlog_v =c(0.2,0.3)), regexp = "numeric and finite")
  expect_error(rlba_lnorm(10, A=c(0.5, 0.6), b=1, t0 = 0.5, meanlog_v=c(1.2, 1), sdlog_v=c(0.2,0.3)), "length") 
  x <- rlba_lnorm(10, A=0.5, b=1, t0 = 0.5, meanlog_v=c(1.2, 1), sdlog_v=c(0.2,0.3))
  expect_error(plba_lnorm(x$response, A=as.character(seq(0.2, 0.5, length.out = 10)), b=0.5, t0 = 0.5, meanlog_v=1.2, sdlog_v=0.2), "numeric")
  expect_is(plba_lnorm(x$response, A=seq(0.2, 0.5, length.out = 10), b=0.5, t0 = seq(0.5, 0.7, 0.1), meanlog_v=1.2, sdlog_v=0.2), "numeric")
})