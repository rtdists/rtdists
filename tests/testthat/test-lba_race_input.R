
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
  
  p1 <- n1PDF(r_lba$rt, A=A[1], b=b[1], t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2], silent = TRUE) # PDF
  p2 <- n1PDF(r_lba$rt, A=list(A[1], A[1]), b=b[1], t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2], silent = TRUE) # PDF
  expect_equal(p2, p1)
  p3 <- n1PDF(r_lba$rt, A=A[1], b=list(b[1], b[1]), t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2], silent = TRUE) # PDF
  expect_equal(p3, p1)
  p4 <- n1PDF(r_lba$rt, A=list(A[1], A[1]), b=list(b[1], b[1]), t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2], silent = TRUE) # PDF
  expect_equal(p4, p1)
  
  p5 <- n1PDF(r_lba$rt, A=rep(A[1], 13), b=list(rep(b[1], 13), rep(b[1], 5)) , t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2], silent = TRUE) # PDF
  expect_equal(p5, p1)
  
  p6_a <- n1PDF(r_lba$rt[seq_len(samples/2)], A=A[1], b=b[1], t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2], silent = TRUE)
  p6_b <- n1PDF(r_lba$rt[seq_len(samples/2)+samples/2], A=A[2], b=b[2], t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2], silent = TRUE)
  p7 <- n1PDF(r_lba$rt, A=rep(A, each = samples/2), b=rep(b, each = samples/2), t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2], silent = TRUE)
  expect_identical(c(p6_a, p6_b), p7)
  
  p6x <- n1PDF(r_lba$rt, A=A[1], b=b[1], t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2], st0 = 0.2, silent = TRUE) # PDF
  p7x <- n1PDF(r_lba$rt, A=list(A[1], A[1]), b=b[1], t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2], st0 = 0.2, silent = TRUE) # PDF
  expect_equal(p7x, p6x)
  p8x <- n1PDF(r_lba$rt, A=A[1], b=list(b[1], b[1]), t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2], st0 = 0.2, silent = TRUE) # PDF
  expect_equal(p8x, p6x)
  p9x <- n1PDF(r_lba$rt, A=list(A[1], A[1]), b=list(b[1], b[1]), t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2], st0 = 0.2, silent = TRUE) # PDF
  expect_equal(p9x, p6x)
  
  p10x <- n1PDF(r_lba$rt, A=rep(A[1], 13), b=list(rep(b[1], 13), rep(b[1], 5)) , t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2], st0 = 0.2, silent = TRUE) # PDF
  expect_equal(p10x, p9x)
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
  
  p1 <- n1CDF(r_lba$rt, A=A[1], b=b[1], t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2], silent = TRUE) # PDF
  p2 <- n1CDF(r_lba$rt, A=list(A[1], A[1]), b=b[1], t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2], silent = TRUE) # PDF
  expect_equal(p2, p1)
  p3 <- n1CDF(r_lba$rt, A=A[1], b=list(b[1], b[1]), t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2], silent = TRUE) # PDF
  expect_equal(p3, p1)
  p4 <- n1CDF(r_lba$rt, A=list(A[1], A[1]), b=list(b[1], b[1]), t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2], silent = TRUE) # PDF
  expect_equal(p4, p1)
  
  p5 <- n1CDF(r_lba$rt, A=rep(A[1], 13), b=list(rep(b[1], 13), rep(b[1], 5)) , t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2], silent = TRUE) # PDF
  expect_equal(p5, p1)
  
  
  p6 <- n1CDF(r_lba$rt, A=A[1], b=b[1], t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2], st0 = 0.2, silent = TRUE) # PDF
  p7 <- n1CDF(r_lba$rt, A=list(A[1], A[1]), b=b[1], t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2], st0 = 0.2, silent = TRUE) # PDF
  expect_equal(p7, p6)
  p8 <- n1CDF(r_lba$rt, A=A[1], b=list(b[1], b[1]), t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2], st0 = 0.2, silent = TRUE) # PDF
  expect_equal(p8, p6)
  p9 <- n1CDF(r_lba$rt, A=list(A[1], A[1]), b=list(b[1], b[1]), t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2], st0 = 0.2, silent = TRUE) # PDF
  expect_equal(p9, p6)
  
  p10 <- n1CDF(r_lba$rt, A=rep(A[1], 13), b=list(rep(b[1], 13), rep(b[1], 5)) , t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2], st0 = 0.2, silent = TRUE) # PDF
  expect_equal(p10, p9)

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


test_that("n1PDF: Trialwise input for t0", {
  samples <- 100
  A <- runif(2, 0.3, 0.9)
  b <- A+runif(2, 0, 0.5)
  t0 <- runif(2, 0.1, 0.7)
  v1 <- runif(4, 0.5, 1.5)
  v2 <- runif(4, 0.1, 0.5)
  st0 <- runif(1, 0.1, 0.5)
  r_lba <- rlba_norm(samples, A=A[1], b=b[1], t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2])
  
  v1_a <- n1PDF(r_lba$rt[seq_len(samples/2)], A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1.0), sd_v=0.2)
  v1_b <- n1PDF(r_lba$rt[seq_len(samples/2)+samples/2], A=0.5, b=1, t0 = 0.2, mean_v=c(1.2, 1.0), sd_v=0.2)
  v2 <- n1PDF(r_lba$rt, A=0.5, b=1, t0 = rep(c(0.5, 0.2), each = 50), mean_v=c(1.2, 1.0), sd_v=0.2)
  expect_identical(c(v1_a, v1_b), v2)
  
  v1 <- runif(4, 0.5, 1.5)
  v2 <- runif(4, 0.5, 1.5)
  st0 <- runif(1, 0.1, 0.5)
  
  r_lba <- rlba_gamma(samples, A=A[1], b=b[1], t0 = t0[1], shape_v =v1[1:2], scale_v = v2[1:2])
  
  v3_a <- n1PDF(r_lba$rt[seq_len(samples/2)], A=0.5, b=1, t0 = 0.5, shape_v=v1[3:4], scale_v = v2[3:4], distribution = "gamma")
  v3_b <- n1PDF(r_lba$rt[seq_len(samples/2)+samples/2], A=0.5, b=1, t0 = 0.2, shape_v=v1[3:4], scale_v = v2[3:4], distribution = "gamma")
  v3 <- n1PDF(r_lba$rt, A=0.5, b=1, t0 = rep(c(0.5, 0.2), each = 50), shape_v=v1[3:4], scale_v = v2[3:4], distribution = "gamma")
  expect_identical(c(v3_a, v3_b), v3)
  
})


test_that("n1PDF: Trialwise input for st0", {
  samples <- 100
  A <- runif(2, 0.3, 0.9)
  b <- A+runif(2, 0, 0.5)
  t0 <- runif(2, 0.1, 0.7)
  v1 <- runif(4, 0.5, 1.5)
  v2 <- runif(4, 0.1, 0.5)
  st0 <- runif(2, 0.1, 0.2)
  r_lba <- rlba_norm(samples, A=A[1], b=b[1], t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2], st0 = st0[1])
  
  v1_a <- n1PDF(r_lba$rt[seq_len(samples/2)], A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1.0), sd_v=0.2, st0 = st0[1])
  v1_b <- n1PDF(r_lba$rt[seq_len(samples/2)+samples/2], A=0.5, b=1, t0 = 0.2, mean_v=c(1.2, 1.0), sd_v=0.2, st0 = st0[2])
  v2 <- n1PDF(r_lba$rt, A=0.5, b=1, t0 = rep(c(0.5, 0.2), each = 50), mean_v=c(1.2, 1.0), sd_v=0.2, st0 = rep(st0, each = 50))
  expect_identical(c(v1_a, v1_b), v2)
  
})

test_that("n1PDF: Trialwise input for drift rate", {
  samples <- 100
  A <- runif(2, 0.3, 0.9)
  b <- A+runif(2, 0, 0.5)
  t0 <- runif(2, 0.1, 0.7)
  v1 <- runif(4, 0.5, 1.5)
  v2 <- runif(4, 0.1, 0.5)
  st0 <- runif(1, 0.1, 0.5)
  r_lba <- rlba_norm(samples, A=A[1], b=b[1], t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2])
  
  v1_a <- n1PDF(r_lba$rt[seq_len(samples/2)], A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1.0), sd_v=0.2)
  v1_b <- n1PDF(r_lba$rt[seq_len(samples/2)+samples/2], A=0.5, b=1, t0 = 0.5, mean_v=c(2.2, 0.5), sd_v=0.5)
  v2 <- n1PDF(r_lba$rt, A=0.5, b=1, t0 = 0.5, mean_v=list(rep(c(1.2, 2.2), each = 50), rep(c(1.0, 0.5), each = 50)), sd_v=list(rep(c(0.2, 0.5), each = 50)))
  expect_identical(c(v1_a, v1_b), v2)
  
  v3_a <- n1PDF(r_lba$rt[seq_len(samples/2)], A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1.0, 0.6), sd_v=0.2)
  v3_b <- n1PDF(r_lba$rt[seq_len(samples/2)+samples/2], A=0.5, b=1, t0 = 0.5, mean_v=c(2.2, 0.5, 1.2), sd_v=0.5)
  v4 <- n1PDF(r_lba$rt, A=0.5, b=1, t0 = 0.5, mean_v=list(rep(c(1.2, 2.2), each = 50), rep(c(1.0, 0.5), each = 50), rep(c(0.6, 1.2), each = 50)), sd_v=list(rep(c(0.2, 0.5), each = 50)))
  expect_identical(c(v3_a, v3_b), v4)
  
  v5_a <- n1PDF(r_lba$rt[seq_len(samples/2)], A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1.0), sd_v=0.2, st0 = 0.2)
  v5_b <- n1PDF(r_lba$rt[seq_len(samples/2)+samples/2], A=0.5, b=1, t0 = 0.5, mean_v=c(2.2, 0.5), sd_v=0.5, st0 = 0.2)
  v6 <- n1PDF(r_lba$rt, A=0.5, b=1, t0 = 0.5, mean_v=list(rep(c(1.2, 2.2), each = 50), rep(c(1.0, 0.5), each = 50)), sd_v=list(rep(c(0.2, 0.5), each = 50)), st0 = 0.2)
  expect_identical(c(v5_a, v5_b), v6)
  
  v7_a <- n1PDF(r_lba$rt[seq_len(samples/2)], A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1.0, 0.6), sd_v=0.2, st0 = 0.2)
  v7_b <- n1PDF(r_lba$rt[seq_len(samples/2)+samples/2], A=0.5, b=1, t0 = 0.5, mean_v=c(2.2, 0.5, 1.2), sd_v=0.5, st0 = 0.2)
  v8 <- n1PDF(r_lba$rt, A=0.5, b=1, t0 = 0.5, mean_v=list(rep(c(1.2, 2.2), each = 50), rep(c(1.0, 0.5), each = 50), rep(c(0.6, 1.2), each = 50)), sd_v=list(rep(c(0.2, 0.5), each = 50)), st0 = 0.2)
  expect_identical(c(v7_a, v7_b), v8)
  
})
