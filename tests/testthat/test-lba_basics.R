

context("LBA works correctly")

test_that("dLBA norm is identical to n1PDF", {
  n <- 100
  x <- rLBA(n, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=c(0.2,0.3))
  
  ex1_n1pdf <- vector("numeric", n)
  ex1_n1pdf[x$response == 1] <- n1PDF(x$rt[x$response == 1], A=0.5, b=1, 
                                      t0 = 0.5, 
                                      mean_v=c(1.2, 1), 
                                      sd_v=c(0.2,0.3), 
                                      distribution = "norm", silent = TRUE)
  ex1_n1pdf[x$response == 2] <- n1PDF(x$rt[x$response == 2], 
                                      A=0.5, b=1, t0 = 0.5, 
                                      mean_v=c(1, 1.2), 
                                      sd_v=c(0.3,0.2), 
                                      distribution = "norm", silent = TRUE)
  ex1_dLBA <- dLBA(x$rt, x$response, A=0.5, b=1, t0 = 0.5, 
                   mean_v=c(1.2, 1), sd_v=c(0.2,0.3), distribution = "norm")
  expect_identical(ex1_n1pdf, ex1_dLBA)
  
  ex2_n1pdf <- vector("numeric", n)
  ex2_n1pdf[x$response == 1] <- n1PDF(x$rt[x$response == 1], 
                                      A=list(0.5, 0.6), 
                                      b=1, t0 = 0.5, 
                                      mean_v=c(1.2, 1), 
                                      sd_v=c(0.2,0.3), 
                                      distribution = "norm", silent = TRUE)
  ex2_n1pdf[x$response == 2] <- n1PDF(x$rt[x$response == 2], 
                                      A=list(0.6, 0.5), 
                                      b=1, t0 = 0.5, 
                                      mean_v=c(1, 1.2), 
                                      sd_v=c(0.3,0.2), 
                                      distribution = "norm", silent = TRUE)
  ex2_dLBA <- dLBA(x$rt, x$response, 
                   A=list(0.5, 0.6), b=1, t0 = 0.5, 
                   mean_v=c(1.2, 1), sd_v=c(0.2,0.3), 
                   distribution = "norm")
  expect_identical(ex2_n1pdf, ex2_dLBA)
  
  ex3_n1pdf <- vector("numeric", n)
  ex3_n1pdf[x$response == 1] <- n1PDF(x$rt[x$response == 1], 
                                      A=rep(c(0.5, 0.6), n/2)[x$response == 1], 
                                      b=1, t0 = 0.5, 
                                      mean_v=c(1.2, 1), 
                                      sd_v=c(0.2,0.3), 
                                      distribution = "norm", silent = TRUE)
  ex3_n1pdf[x$response == 2] <- n1PDF(x$rt[x$response == 2], 
                                      A=rep(c(0.5, 0.6), n/2)[x$response == 2], 
                                      b=1, t0 = 0.5, 
                                      mean_v=c(1, 1.2), 
                                      sd_v=c(0.3,0.2), 
                                      distribution = "norm", silent = TRUE)
  ex3_dLBA <- dLBA(x$rt, x$response, A=c(0.5, 0.6), b=1, t0 = 0.5, 
                   mean_v=c(1.2, 1), sd_v=c(0.2,0.3), 
                   distribution = "norm")
  expect_identical(ex3_n1pdf, ex3_dLBA)
  
  ex4_n1pdf <- vector("numeric", n)
  ex4_n1pdf[x$response == 1] <- n1PDF(x$rt[x$response == 1], 
                                      A=list(
                                        rep(c(0.5, 0.6), n/2)[x$response == 1], 
                                        rep(c(0.6, 0.5), n/2)[x$response == 1]
                                      ), 
                                      b=1, t0 = 0.5, 
                                      mean_v=c(1.2, 1), 
                                      sd_v=c(0.2,0.3), 
                                      distribution = "norm", silent = TRUE)
  ex4_n1pdf[x$response == 2] <- n1PDF(x$rt[x$response == 2], 
                                      A=list(
                                        rep(c(0.6, 0.5), n/2)[x$response == 2], 
                                        rep(c(0.5, 0.6), n/2)[x$response == 2]
                                      ),  
                                      b=1, t0 = 0.5, 
                                      mean_v=c(1, 1.2), 
                                      sd_v=c(0.3,0.2), 
                                      distribution = "norm", silent = TRUE)
  ex4_dLBA <- dLBA(x$rt, x$response, 
                   A=list(c(0.5, 0.6), c(0.6, 0.5)), b=1, t0 = 0.5, 
                   mean_v=c(1.2, 1), sd_v=c(0.2,0.3), 
                   distribution = "norm")
  expect_identical(ex4_n1pdf, ex4_dLBA)
  
  ex5_n1pdf <- vector("numeric", n)
  ex5_n1pdf[x$response == 1] <- n1PDF(x$rt[x$response == 1], 
                                      A=list(
                                        rep(c(0.5, 0.6), n/2)[x$response == 1], 
                                        rep(c(0.6, 0.5), n/2)[x$response == 1]
                                      ),
                                      b=1, 
                                      t0 = list(0.5, 0.3), 
                                      mean_v=c(1.2, 1), 
                                      sd_v=c(0.2,0.3), 
                                      distribution = "norm", silent = TRUE)
  ex5_n1pdf[x$response == 2] <- n1PDF(x$rt[x$response == 2], 
                                      A=list(
                                        rep(c(0.6, 0.5), n/2)[x$response == 2], 
                                        rep(c(0.5, 0.6), n/2)[x$response == 2]
                                      ), 
                                      b=1, 
                                      t0 = list(0.3, 0.5), 
                                      mean_v=c(1, 1.2), 
                                      sd_v=c(0.3,0.2), 
                                      distribution = "norm", silent = TRUE)
  ex5_dLBA <- dLBA(x$rt, x$response, 
                   A=list(c(0.5, 0.6), c(0.6, 0.5)), 
                   b=1, t0 = list(0.5, 0.3), 
                   mean_v=c(1.2, 1), 
                   sd_v=c(0.2,0.3), 
                   distribution = "norm")
  expect_identical(ex5_n1pdf, ex5_dLBA)
  
  ex6_n1pdf <- vector("numeric", n)
  ex6_n1pdf[x$response == 1] <- n1PDF(x$rt[x$response == 1], 
                                      A=list(
                                        rep(c(0.5, 0.6), n/2)[x$response == 1], 
                                        rep(c(0.6, 0.5), n/2)[x$response == 1]
                                      ), 
                                      b=1, t0 = list(0.5, 0.3), 
                                      mean_v=c(1.2, 1), 
                                      sd_v=c(0.2,0.3), 
                                      distribution = "norm", 
                                      silent = TRUE, st0 = 0.2)
  ex6_n1pdf[x$response == 2] <- n1PDF(x$rt[x$response == 2], 
                                      A=list(
                                        rep(c(0.6, 0.5), n/2)[x$response == 2], 
                                        rep(c(0.5, 0.6), n/2)[x$response == 2]
                                      ), 
                                      b=1, t0 = list(0.3, 0.5), 
                                      mean_v=c(1, 1.2), 
                                      sd_v=c(0.3,0.2), 
                                      distribution = "norm", silent = TRUE, 
                                      st0 = 0.2)
  ex6_dLBA <- dLBA(x$rt, x$response, 
                   A=list(c(0.5, 0.6), c(0.6, 0.5)), 
                   b=1, t0 = list(0.5, 0.3), 
                   mean_v=c(1.2, 1), 
                   sd_v=c(0.2,0.3), 
                   distribution = "norm", st0 = 0.2)
  expect_identical(ex6_n1pdf, ex6_dLBA)
  
  
  drifts1 <- rnorm(n/2, 1, 0.1)
  drifts2 <- rnorm(n/2, 1.5, 0.2)
  
  sddrift1 <- runif(n/5, 0.1, 0.3)
  sddrift2 <- runif(n/5, 0.2, 0.4)
    
  ex7_n1pdf <- vector("numeric", n)
  ex7_n1pdf[x$response == 1] <- n1PDF(x$rt[x$response == 1], 
                                      A=list(
                                        rep(c(0.5, 0.6), n/2)[x$response == 1], 
                                        rep(c(0.6, 0.5), n/2)[x$response == 1]
                                      ), 
                                      b = list(
                                        rep(c(1, 1.1, 1.2, 1.3, 1.4), n/5)[x$response == 1], 
                                        rep(c(1, 1.1, 1.2, 1.3, 1.4) - 0.2, n/5)[x$response == 1]
                                      ),
                                      t0 = list(
                                        seq(0.1, 0.5, length.out = n)[x$response == 1], 
                                        seq(0.3, 0.8, length.out = n)[x$response == 1]), 
                                      mean_v = list(
                                        rep(drifts1, 2)[x$response == 1], 
                                        rep(drifts2, 2)[x$response == 1]
                                      ), 
                                      sd_v = list(
                                        rep(sddrift1, 5)[x$response == 1], 
                                        rep(sddrift2, 5)[x$response == 1]
                                      ), 
                                      distribution = "norm", 
                                      silent = TRUE, st0 = 0.2)
  ex7_n1pdf[x$response == 2] <- n1PDF(x$rt[x$response == 2], 
                                      A=list(
                                        rep(c(0.6, 0.5), n/2)[x$response == 2], 
                                        rep(c(0.5, 0.6), n/2)[x$response == 2]
                                      ), 
                                      b = list(
                                        rep(c(1, 1.1, 1.2, 1.3, 1.4) - 0.2, n/5)[x$response == 2],
                                        rep(c(1, 1.1, 1.2, 1.3, 1.4), n/5)[x$response == 2]
                                      ),
                                      t0 = list(
                                        seq(0.3, 0.8, length.out = n)[x$response == 2],
                                        seq(0.1, 0.5, length.out = n)[x$response == 2]
                                      ), 
                                      mean_v = list(
                                        rep(drifts2, 2)[x$response == 2], 
                                        rep(drifts1, 2)[x$response == 2]
                                      ), 
                                      sd_v = list(
                                        rep(sddrift2, 5)[x$response == 2], 
                                        rep(sddrift1, 5)[x$response == 2]
                                      ), 
                                      distribution = "norm", silent = TRUE, 
                                      st0 = 0.2)
  ex7_dLBA <- dLBA(x$rt, x$response, 
                   A = list(
                     c(0.5, 0.6), 
                     c(0.6, 0.5)
                   ),
                   b = list(
                     c(1, 1.1, 1.2, 1.3, 1.4), 
                     c(1, 1.1, 1.2, 1.3, 1.4) - 0.2
                   ),
                   t0 = list(
                     seq(0.1, 0.5, length.out = n), 
                     seq(0.3, 0.8, length.out = n)), 
                   mean_v = list(
                     drifts1, 
                     drifts2
                   ), 
                   sd_v = list(
                     sddrift1, 
                     sddrift2
                   ), 
                   distribution = "norm", st0 = 0.2)
  expect_identical(ex7_n1pdf, ex7_dLBA)
  
})

test_that("pLBA norm is identical to n1CDF", {
  n <- 100
  x <- rLBA(n, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1), sd_v=c(0.2,0.3))
  
  ex1_n1cdf <- vector("numeric", n)
  ex1_n1cdf[x$response == 1] <- n1CDF(x$rt[x$response == 1], A=0.5, b=1, 
                                      t0 = 0.5, 
                                      mean_v=c(1.2, 1), 
                                      sd_v=c(0.2,0.3), 
                                      distribution = "norm", silent = TRUE)
  ex1_n1cdf[x$response == 2] <- n1CDF(x$rt[x$response == 2], 
                                      A=0.5, b=1, t0 = 0.5, 
                                      mean_v=c(1, 1.2), 
                                      sd_v=c(0.3,0.2), 
                                      distribution = "norm", silent = TRUE)
  ex1_pLBA <- pLBA(x$rt, x$response, A=0.5, b=1, t0 = 0.5, 
                   mean_v=c(1.2, 1), sd_v=c(0.2,0.3), distribution = "norm")
  expect_identical(ex1_n1cdf, ex1_pLBA)
  
  ex2_n1cdf <- vector("numeric", n)
  ex2_n1cdf[x$response == 1] <- n1CDF(x$rt[x$response == 1], 
                                      A=list(0.5, 0.6), 
                                      b=1, t0 = 0.5, 
                                      mean_v=c(1.2, 1), 
                                      sd_v=c(0.2,0.3), 
                                      distribution = "norm", silent = TRUE)
  ex2_n1cdf[x$response == 2] <- n1CDF(x$rt[x$response == 2], 
                                      A=list(0.6, 0.5), 
                                      b=1, t0 = 0.5, 
                                      mean_v=c(1, 1.2), 
                                      sd_v=c(0.3,0.2), 
                                      distribution = "norm", silent = TRUE)
  ex2_pLBA <- pLBA(x$rt, x$response, 
                   A=list(0.5, 0.6), b=1, t0 = 0.5, 
                   mean_v=c(1.2, 1), sd_v=c(0.2,0.3), 
                   distribution = "norm")
  expect_identical(ex2_n1cdf, ex2_pLBA)
  
  ex3_n1cdf <- vector("numeric", n)
  ex3_n1cdf[x$response == 1] <- n1CDF(x$rt[x$response == 1], 
                                      A=rep(c(0.5, 0.6), n/2)[x$response == 1], 
                                      b=1, t0 = 0.5, 
                                      mean_v=c(1.2, 1), 
                                      sd_v=c(0.2,0.3), 
                                      distribution = "norm", silent = TRUE)
  ex3_n1cdf[x$response == 2] <- n1CDF(x$rt[x$response == 2], 
                                      A=rep(c(0.5, 0.6), n/2)[x$response == 2], 
                                      b=1, t0 = 0.5, 
                                      mean_v=c(1, 1.2), 
                                      sd_v=c(0.3,0.2), 
                                      distribution = "norm", silent = TRUE)
  ex3_pLBA <- pLBA(x$rt, x$response, A=c(0.5, 0.6), b=1, t0 = 0.5, 
                   mean_v=c(1.2, 1), sd_v=c(0.2,0.3), 
                   distribution = "norm")
  expect_identical(ex3_n1cdf, ex3_pLBA)
  
  ex4_n1cdf <- vector("numeric", n)
  ex4_n1cdf[x$response == 1] <- n1CDF(x$rt[x$response == 1], 
                                      A=list(
                                        rep(c(0.5, 0.6), n/2)[x$response == 1], 
                                        rep(c(0.6, 0.5), n/2)[x$response == 1]
                                      ), 
                                      b=1, t0 = 0.5, 
                                      mean_v=c(1.2, 1), 
                                      sd_v=c(0.2,0.3), 
                                      distribution = "norm", silent = TRUE)
  ex4_n1cdf[x$response == 2] <- n1CDF(x$rt[x$response == 2], 
                                      A=list(
                                        rep(c(0.6, 0.5), n/2)[x$response == 2], 
                                        rep(c(0.5, 0.6), n/2)[x$response == 2]
                                      ),  
                                      b=1, t0 = 0.5, 
                                      mean_v=c(1, 1.2), 
                                      sd_v=c(0.3,0.2), 
                                      distribution = "norm", silent = TRUE)
  ex4_pLBA <- pLBA(x$rt, x$response, 
                   A=list(c(0.5, 0.6), c(0.6, 0.5)), b=1, t0 = 0.5, 
                   mean_v=c(1.2, 1), sd_v=c(0.2,0.3), 
                   distribution = "norm")
  expect_identical(ex4_n1cdf, ex4_pLBA)
  
  ex5_n1cdf <- vector("numeric", n)
  ex5_n1cdf[x$response == 1] <- n1CDF(x$rt[x$response == 1], 
                                      A=list(
                                        rep(c(0.5, 0.6), n/2)[x$response == 1], 
                                        rep(c(0.6, 0.5), n/2)[x$response == 1]
                                      ),
                                      b=1, 
                                      t0 = list(0.5, 0.3), 
                                      mean_v=c(1.2, 1), 
                                      sd_v=c(0.2,0.3), 
                                      distribution = "norm", silent = TRUE)
  ex5_n1cdf[x$response == 2] <- n1CDF(x$rt[x$response == 2], 
                                      A=list(
                                        rep(c(0.6, 0.5), n/2)[x$response == 2], 
                                        rep(c(0.5, 0.6), n/2)[x$response == 2]
                                      ), 
                                      b=1, 
                                      t0 = list(0.3, 0.5), 
                                      mean_v=c(1, 1.2), 
                                      sd_v=c(0.3,0.2), 
                                      distribution = "norm", silent = TRUE)
  ex5_pLBA <- pLBA(x$rt, x$response, 
                   A=list(c(0.5, 0.6), c(0.6, 0.5)), 
                   b=1, t0 = list(0.5, 0.3), 
                   mean_v=c(1.2, 1), 
                   sd_v=c(0.2,0.3), 
                   distribution = "norm")
  expect_identical(ex5_n1cdf, ex5_pLBA)
  
  ex6_n1cdf <- vector("numeric", n)
  ex6_n1cdf[x$response == 1] <- n1CDF(x$rt[x$response == 1], 
                                      A=list(
                                        rep(c(0.5, 0.6), n/2)[x$response == 1], 
                                        rep(c(0.6, 0.5), n/2)[x$response == 1]
                                      ), 
                                      b=1, t0 = list(0.5, 0.3), 
                                      mean_v=c(1.2, 1), 
                                      sd_v=c(0.2,0.3), 
                                      distribution = "norm", 
                                      silent = TRUE, st0 = 0.2)
  ex6_n1cdf[x$response == 2] <- n1CDF(x$rt[x$response == 2], 
                                      A=list(
                                        rep(c(0.6, 0.5), n/2)[x$response == 2], 
                                        rep(c(0.5, 0.6), n/2)[x$response == 2]
                                      ), 
                                      b=1, t0 = list(0.3, 0.5), 
                                      mean_v=c(1, 1.2), 
                                      sd_v=c(0.3,0.2), 
                                      distribution = "norm", silent = TRUE, 
                                      st0 = 0.2)
  ex6_pLBA <- pLBA(x$rt, x$response, 
                   A=list(c(0.5, 0.6), c(0.6, 0.5)), 
                   b=1, t0 = list(0.5, 0.3), 
                   mean_v=c(1.2, 1), 
                   sd_v=c(0.2,0.3), 
                   distribution = "norm", st0 = 0.2)
  expect_identical(ex6_n1cdf, ex6_pLBA)
  
  
  drifts1 <- rnorm(n/2, 1, 0.1)
  drifts2 <- rnorm(n/2, 1.5, 0.2)
  
  sddrift1 <- runif(n/5, 0.1, 0.3)
  sddrift2 <- runif(n/5, 0.2, 0.4)
    
  ex7_n1cdf <- vector("numeric", n)
  ex7_n1cdf[x$response == 1] <- n1CDF(x$rt[x$response == 1], 
                                      A=list(
                                        rep(c(0.5, 0.6), n/2)[x$response == 1], 
                                        rep(c(0.6, 0.5), n/2)[x$response == 1]
                                      ), 
                                      b = list(
                                        rep(c(1, 1.1, 1.2, 1.3, 1.4), n/5)[x$response == 1], 
                                        rep(c(1, 1.1, 1.2, 1.3, 1.4) - 0.2, n/5)[x$response == 1]
                                      ),
                                      t0 = list(
                                        seq(0.1, 0.5, length.out = n)[x$response == 1], 
                                        seq(0.3, 0.8, length.out = n)[x$response == 1]), 
                                      mean_v = list(
                                        rep(drifts1, 2)[x$response == 1], 
                                        rep(drifts2, 2)[x$response == 1]
                                      ), 
                                      sd_v = list(
                                        rep(sddrift1, 5)[x$response == 1], 
                                        rep(sddrift2, 5)[x$response == 1]
                                      ), 
                                      distribution = "norm", 
                                      silent = TRUE, st0 = 0.2)
  ex7_n1cdf[x$response == 2] <- n1CDF(x$rt[x$response == 2], 
                                      A=list(
                                        rep(c(0.6, 0.5), n/2)[x$response == 2], 
                                        rep(c(0.5, 0.6), n/2)[x$response == 2]
                                      ), 
                                      b = list(
                                        rep(c(1, 1.1, 1.2, 1.3, 1.4) - 0.2, n/5)[x$response == 2],
                                        rep(c(1, 1.1, 1.2, 1.3, 1.4), n/5)[x$response == 2]
                                      ),
                                      t0 = list(
                                        seq(0.3, 0.8, length.out = n)[x$response == 2],
                                        seq(0.1, 0.5, length.out = n)[x$response == 2]
                                      ), 
                                      mean_v = list(
                                        rep(drifts2, 2)[x$response == 2], 
                                        rep(drifts1, 2)[x$response == 2]
                                      ), 
                                      sd_v = list(
                                        rep(sddrift2, 5)[x$response == 2], 
                                        rep(sddrift1, 5)[x$response == 2]
                                      ), 
                                      distribution = "norm", silent = TRUE, 
                                      st0 = 0.2)
  ex7_pLBA <- pLBA(x$rt, x$response, 
                   A = list(
                     c(0.5, 0.6), 
                     c(0.6, 0.5)
                   ),
                   b = list(
                     c(1, 1.1, 1.2, 1.3, 1.4), 
                     c(1, 1.1, 1.2, 1.3, 1.4) - 0.2
                   ),
                   t0 = list(
                     seq(0.1, 0.5, length.out = n), 
                     seq(0.3, 0.8, length.out = n)), 
                   mean_v = list(
                     drifts1, 
                     drifts2
                   ), 
                   sd_v = list(
                     sddrift1, 
                     sddrift2
                   ), 
                   distribution = "norm", st0 = 0.2)
  expect_identical(ex7_n1cdf, ex7_pLBA)
  
})


test_that("pLBA norm is identical to n1CDF", {
  x <- seq(0, 3, by =0.1)
  
  o1a <- n1CDF(x, A=0.5, b=1, t0 = 0.5, 
               mean_v=c(1.2, 1), 
               sd_v=c(0.2,0.3), 
               distribution = "norm", silent = TRUE)
  o1b <- n1CDF(x, A=0.5, b=1, t0 = 0.5, 
               mean_v=c(1, 1.2), 
               sd_v=c(0.3,0.2), 
               distribution = "norm", silent = TRUE)
  o1c <- pLBA(c(x, x), rep(1:2, each = length(x)), 
              A=0.5, b=1, t0 = 0.5, 
              mean_v=c(1.2, 1), 
              sd_v=c(0.2,0.3), 
              distribution = "norm", silent = TRUE)
  expect_identical(c(o1a, o1b), o1c)
  
  o2a <- n1CDF(x, A=0.5, b=1, t0 = list(0.5, 0.2), 
               mean_v=c(1.2, 1), 
               sd_v=c(0.2,0.3), 
               distribution = "norm", silent = TRUE)
  o2b <- n1CDF(x, A=0.5, b=1, t0 = list(0.2, 0.5), 
               mean_v=c(1, 1.2), 
               sd_v=c(0.3,0.2), 
               distribution = "norm", silent = TRUE)
  o2c <- pLBA(c(x, x), rep(1:2, each = length(x)), 
              A=0.5, b=1, t0 = list(0.5, 0.2), 
              mean_v=c(1.2, 1), 
              sd_v=c(0.2,0.3), 
              distribution = "norm", silent = TRUE)
  expect_identical(c(o2a, o2b), o2c)
  
  o3a <- n1CDF(x, A=list(seq(0.5, 0.6, length.out = length(x)), 0.2), 
               b=1, t0 = list(0.5, 0.2), 
               mean_v=c(1.2, 1), 
               sd_v=c(0.2,0.3), 
               distribution = "norm", silent = TRUE)
  o3b <- n1CDF(x, A=list(0.2, seq(0.5, 0.6, length.out = length(x))), 
               b=1, t0 = list(0.2, 0.5), 
               mean_v=c(1, 1.2), 
               sd_v=c(0.3,0.2), 
               distribution = "norm", silent = TRUE)
  o3c <- pLBA(c(x, x), rep(1:2, each = length(x)), 
              A=list(seq(0.5, 0.6, length.out = length(x)), 0.2), 
              b=1, t0 = list(0.5, 0.2), 
              mean_v=c(1.2, 1), 
              sd_v=c(0.2,0.3), 
              distribution = "norm", silent = TRUE)
  expect_identical(c(o2a, o2b), o2c)
  
})

test_that("qLBA is equivalent to manual calculation",{
  p11_norm <- structure(list(par =
                               structure(
                                 c(
                                   0.11829564514883,
                                   -2.74097628458775,
                                   1.04498350371894,
                                   0.451351158199909,
                                   0.124346087574482,
                                   0.260994169758609
                                 ),
                                 .Names = c("A", "v1", "v2", "b", "t0", "sv")
                               )))
  q <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  pred_prop_correct <- pLBA(Inf, 2, 
                            A = p11_norm$par["A"], 
                            b = p11_norm$par["A"]+p11_norm$par["b"], 
                            t0 = p11_norm$par["t0"], 
                            mean_v = c(p11_norm$par["v1"], p11_norm$par["v2"]), 
                            sd_v = c(1, p11_norm$par["sv"]), 
                            silent = TRUE)

  expect_equal(qLBA(pred_prop_correct*q, 2, 
                    A = p11_norm$par["A"], 
                    b = p11_norm$par["A"]+p11_norm$par["b"], 
                    t0 = p11_norm$par["t0"], 
                    mean_v = c(p11_norm$par["v1"], p11_norm$par["v2"]), 
                    sd_v = c(1, p11_norm$par["sv"]), silent = TRUE),
               c(0.487170939752454, 0.551026400837336, 0.608185370581083, 
                 0.680979476696082, 0.830128589908231))
  
  expect_equal(suppressWarnings(qLBA(pred_prop_correct*q, 1, 
                                     A = p11_norm$par["A"], 
                                     b = p11_norm$par["A"]+p11_norm$par["b"], 
                                     t0 = p11_norm$par["t0"], 
                                     mean_v = c(p11_norm$par["v1"], 
                                                p11_norm$par["v2"]), 
                                     sd_v = c(1, p11_norm$par["sv"]), 
                                     silent = TRUE)),
               as.numeric(rep(NA, 5)))
})

test_that("rLBA works wit all distributions", {
  expect_is(rLBA(100, A=c(0.5, 0.6), b=1, t0 = 0.5, 
                 meanlog_v=c(1.2, 1), sdlog_v=c(0.2,0.3), 
                 distribution = "lnorm"), 
            "data.frame")
  expect_is(rLBA(100, A=c(0.5, 0.6), b=1, t0 = 0.5, 
                 shape_v=c(1.2, 1), scale_v=c(0.2,0.3), 
                 distribution = "frechet"), 
            "data.frame")
  expect_is(rLBA(100, A=c(0.5, 0.6), b=1, t0 = 0.5, 
                 shape_v=c(1.2, 1), scale_v=c(0.2,0.3), 
                 distribution = "gamma"), 
            "data.frame")
})
