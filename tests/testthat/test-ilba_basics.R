

context("iLBA works correctly")

test_that("diLBA norm is identical to n1PDF", {
  n <- 100
  x <- riLBA(n, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=c(0.2,0.3))
  
  ex1_n1pdf <- vector("numeric", n)
  ex1_n1pdf[x$response == 1] <- n1PDF(x$rt[x$response == 1], A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=c(0.2,0.3), distribution = "norm", silent = TRUE)
  ex1_n1pdf[x$response == 2] <- n1PDF(x$rt[x$response == 2], A=0.5, b=1, t0 = 0.5, mean_v=c(1, 1.2), sd_v=c(0.3,0.2), distribution = "norm", silent = TRUE)
  ex1_dilba <- diLBA(x$rt, x$response, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=c(0.2,0.3), distribution = "norm")
  expect_identical(ex1_n1pdf, ex1_dilba)
  
  ex2_n1pdf <- vector("numeric", n)
  ex2_n1pdf[x$response == 1] <- n1PDF(x$rt[x$response == 1], A=list(0.5, 0.6), b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=c(0.2,0.3), distribution = "norm", silent = TRUE)
  ex2_n1pdf[x$response == 2] <- n1PDF(x$rt[x$response == 2], A=list(0.6, 0.5), b=1, t0 = 0.5, mean_v=c(1, 1.2), sd_v=c(0.3,0.2), distribution = "norm", silent = TRUE)
  ex2_dilba <- diLBA(x$rt, x$response, A=list(0.5, 0.6), b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=c(0.2,0.3), distribution = "norm")
  expect_identical(ex2_n1pdf, ex2_dilba)
  
  ex3_n1pdf <- vector("numeric", n)
  ex3_n1pdf[x$response == 1] <- n1PDF(x$rt[x$response == 1], A=c(0.5, 0.6), b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=c(0.2,0.3), distribution = "norm", silent = TRUE)
  ex3_n1pdf[x$response == 2] <- n1PDF(x$rt[x$response == 2], A=c(0.5, 0.6), b=1, t0 = 0.5, mean_v=c(1, 1.2), sd_v=c(0.3,0.2), distribution = "norm", silent = TRUE)
  ex3_dilba <- diLBA(x$rt, x$response, A=c(0.5, 0.6), b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=c(0.2,0.3), distribution = "norm")
  expect_identical(ex3_n1pdf, ex3_dilba)
  
  ex4_n1pdf <- vector("numeric", n)
  ex4_n1pdf[x$response == 1] <- n1PDF(x$rt[x$response == 1], A=list(c(0.5, 0.6), c(0.6, 0.5)), b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=c(0.2,0.3), distribution = "norm", silent = TRUE)
  ex4_n1pdf[x$response == 2] <- n1PDF(x$rt[x$response == 2], A=list(c(0.6, 0.5), c(0.5, 0.6)), b=1, t0 = 0.5, mean_v=c(1, 1.2), sd_v=c(0.3,0.2), distribution = "norm", silent = TRUE)
  ex4_dilba <- diLBA(x$rt, x$response, A=list(c(0.5, 0.6), c(0.6, 0.5)), b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=c(0.2,0.3), distribution = "norm")
  expect_identical(ex4_n1pdf, ex4_dilba)
  
})


test_that("piLBA norm is identical to n1CDF", {
  x <- seq(0, 3, by =0.1)
  
  o1a <- n1CDF(x, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=c(0.2,0.3), distribution = "norm", silent = TRUE)
  o1b <- n1CDF(x, A=0.5, b=1, t0 = 0.5, mean_v=c(1, 1.2), sd_v=c(0.3,0.2), distribution = "norm", silent = TRUE)
  o1c <- piLBA(c(x, x), rep(1:2, each = length(x)), A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=c(0.2,0.3), distribution = "norm", silent = TRUE)
  expect_identical(c(o1a, o1b), o1c)
  
  o2a <- n1CDF(x, A=0.5, b=1, t0 = list(0.5, 0.2), mean_v=c(1.2, 1), sd_v=c(0.2,0.3), distribution = "norm", silent = TRUE)
  o2b <- n1CDF(x, A=0.5, b=1, t0 = list(0.2, 0.5), mean_v=c(1, 1.2), sd_v=c(0.3,0.2), distribution = "norm", silent = TRUE)
  o2c <- piLBA(c(x, x), rep(1:2, each = length(x)), A=0.5, b=1, t0 = list(0.5, 0.2), mean_v=c(1.2, 1), sd_v=c(0.2,0.3), distribution = "norm", silent = TRUE)
  expect_identical(c(o2a, o2b), o2c)
  
  o3a <- n1CDF(x, A=list(seq(0.5, 0.6, length.out = length(x)), 0.2), b=1, t0 = list(0.5, 0.2), mean_v=c(1.2, 1), sd_v=c(0.2,0.3), distribution = "norm", silent = TRUE)
  o3b <- n1CDF(x, A=list(0.2, seq(0.5, 0.6, length.out = length(x))), b=1, t0 = list(0.2, 0.5), mean_v=c(1, 1.2), sd_v=c(0.3,0.2), distribution = "norm", silent = TRUE)
  o3c <- piLBA(c(x, x), rep(1:2, each = length(x)), A=list(seq(0.5, 0.6, length.out = length(x)), 0.2), b=1, t0 = list(0.5, 0.2), mean_v=c(1.2, 1), sd_v=c(0.2,0.3), distribution = "norm", silent = TRUE)
  expect_identical(c(o2a, o2b), o2c)
  
})

