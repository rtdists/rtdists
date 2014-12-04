
context("LBA race functions: Input")

test_that("n1PDF: List input for A and b", {
  samples <- 100
  A <- runif(2, 0.3, 0.9)
  b <- A+runif(2, 0, 0.5)
  t0 <- runif(2, 0.1, 0.7)
  v1 <- runif(4, 0.5, 1.5)
  v2 <- runif(4, 0.1, 0.5)
  st0 <- runif(1, 0.1, 0.5)
  r_lba <- rlba_norm(samples, A=A[1], b=b[1], t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2])
  
  p1 <- n1PDF(r_lba$rt[r_lba$response==1], A=A[1], b=b[1], t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2]) # PDF
  p2 <- n1PDF(r_lba$rt[r_lba$response==1], A=list(A[1], A[1]), b=b[1], t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2]) # PDF
  expect_equal(p2, p1)
  p3 <- n1PDF(r_lba$rt[r_lba$response==1], A=A[1], b=list(b[1], b[1]), t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2]) # PDF
  expect_equal(p3, p1)
  p4 <- n1PDF(r_lba$rt[r_lba$response==1], A=list(A[1], A[1]), b=list(b[1], b[1]), t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2]) # PDF
  expect_equal(p4, p1)
  
  p5 <- n1PDF(r_lba$rt[r_lba$response==1], A=rep(A[1], 13), b=list(rep(b[1], 13), rep(b[1], 5)) , t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2]) # PDF
  expect_equal(p5, p1)
  
})


test_that("n1CDF: List input for A and b", {
  samples <- 100
  A <- runif(2, 0.3, 0.9)
  b <- A+runif(2, 0, 0.5)
  t0 <- runif(2, 0.1, 0.7)
  v1 <- runif(4, 0.5, 1.5)
  v2 <- runif(4, 0.1, 0.5)
  st0 <- runif(1, 0.1, 0.5)
  r_lba <- rlba_norm(samples, A=A[1], b=b[1], t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2])
  
  p1 <- n1CDF(r_lba$rt[r_lba$response==1], A=A[1], b=b[1], t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2]) # PDF
  p2 <- n1CDF(r_lba$rt[r_lba$response==1], A=list(A[1], A[1]), b=b[1], t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2]) # PDF
  expect_equal(p2, p1)
  p3 <- n1CDF(r_lba$rt[r_lba$response==1], A=A[1], b=list(b[1], b[1]), t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2]) # PDF
  expect_equal(p3, p1)
  p4 <- n1CDF(r_lba$rt[r_lba$response==1], A=list(A[1], A[1]), b=list(b[1], b[1]), t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2]) # PDF
  expect_equal(p4, p1)
  
  p5 <- n1CDF(r_lba$rt[r_lba$response==1], A=rep(A[1], 13), b=list(rep(b[1], 13), rep(b[1], 5)) , t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2]) # PDF
  expect_equal(p5, p1)
  
})

test_that("n1CDF: List input for drift rate", {
  samples <- 100
  A <- runif(2, 0.3, 0.9)
  b <- A+runif(2, 0, 0.5)
  t0 <- runif(2, 0.1, 0.7)
  v1 <- runif(4, 0.5, 1.5)
  v2 <- runif(4, 0.1, 0.5)
  st0 <- runif(1, 0.1, 0.5)
  r_lba <- rlba_norm(samples, A=A[1], b=b[1], t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2])
  
  v1 <- n1PDF(r_lba$rt, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1.0), sd_v=0.2)
  v2 <- n1PDF(r_lba$rt, A=0.5, b=1, t0 = 0.5, mean_v=list(1.2, 1.0), sd_v=0.2)
  expect_identical(v1, v2)
  
})


