

context("iLBA works correctly")

test_that("norm", {
  n <- 100
  x <- rlba_norm(n, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=c(0.2,0.3))
  
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
