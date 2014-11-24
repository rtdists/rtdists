

context("LBA vectorization")

test_that("_norm", {
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


context("LBA vectorization")

test_that("_norm", {
  expect_error(rlba_norm("10", A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=c(0.2,0.3)))
  expect_error(rlba_norm(10, A=c(0.5, 0.6), b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=c(0.2,0.3)))    
  x <- rlba_norm(10, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=c(0.2,0.3))
  expect_error(plba_norm(x$response, A=as.character(seq(0.2, 0.5, length.out = 10)), b=0.5, t0 = 0.5, mean_v=1.2, sd_v=0.2))
  expect_is(plba_norm(x$response, A=seq(0.2, 0.5, length.out = 10), b=0.5, t0 = 0.5, mean_v=1.2, sd_v=0.2), "numeric")
})
