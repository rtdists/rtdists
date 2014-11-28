
context("LBA race functions")

test_that("Norm: n1CDF corresponds to random derivates", {
  normalised_n1CDF = function(t,...) n1CDF(t,...)/n1CDF(t=Inf,...) 
  samples <- 2e3
  p_min <- 0.0001
  p_max <- 0.05
  A <- runif(2, 0.3, 0.9)
  b <- runif(2, 0.5, 1.5)
  t0 <- runif(2, 0.1, 0.7)
  v1 <- runif(4, 0.5, 1.5)
  v2 <- runif(4, 0.1, 0.5)
  st0 <- runif(1, 0.1, 0.5)
  r_lba1 <- rlba_norm(samples, A=A[1], b=b[1], t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2])
  r_lba2 <- rlba_norm(samples, A=A[2], b=b[2], t0 = t0[2], mean_v=v1[3:4], sd_v=v2[3:4], st0 = st0[1])
  expect_less_than(ks.test(r_lba1$rt[r_lba1$response==1], "n1CDF",A=A[1], b=b[1], t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2])$p.value, p_min)
  expect_less_than(ks.test(r_lba1$rt[r_lba1$response==1], normalised_n1CDF, A=A[1]+0.1, b=b[1], t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2])$p.value, p_min)
  expect_less_than(ks.test(r_lba2$rt[r_lba2$response==1], normalised_n1CDF, A=A[2], b=b[2], t0 = t0[2], mean_v=v1[3:4], sd_v=v2[3:4], st0 = 0)$p.value, p_min)
  expect_less_than(ks.test(r_lba2$rt[r_lba2$response==1], normalised_n1CDF, A=A[2], b=b[2], t0 = t0[2], mean_v=v1[3:4], sd_v=v2[3:4], st0 = st0[1]+0.1)$p.value, p_min)
  
  expect_more_than(ks.test(r_lba1$rt[r_lba1$response==1], normalised_n1CDF, A=A[1], b=b[1], t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2])$p.value, p_max)
  expect_more_than(ks.test(r_lba2$rt[r_lba2$response==1], normalised_n1CDF, A=A[2], b=b[2], t0 = t0[2], mean_v=v1[3:4], sd_v=v2[3:4], st0 = st0[1])$p.value, p_max)
})

test_that("Gamma: n1CDF corresponds to random derivates", {
  normalised_n1CDF = function(t,...) n1CDF(t,...)/n1CDF(t=Inf,...) 
  samples <- 2e3
  p_min <- 0.001
  p_max <- 0.05
  A <- runif(2, 0.3, 0.9)
  b <- runif(2, 0.5, 1.5)
  t0 <- runif(2, 0.1, 0.7)
  v1 <- runif(4, 0.5, 1.5)
  v2 <- runif(4, 0.1, 0.5)
  st0 <- runif(1, 0.1, 0.5)
  r_lba1 <- rlba_gamma(samples, A=A[1], b=b[1], t0 = t0[1], shape_v=v1[1:2], scale_v=v2[1:2])
  r_lba2 <- rlba_gamma(samples, A=A[2], b=b[2], t0 = t0[2], shape_v=v1[3:4], scale_v=v2[3:4], st0 = st0[1])
  expect_less_than(ks.test(r_lba1$rt[r_lba1$response==1], "n1CDF",A=A[1], b=b[1], t0 = t0[1], shape_v=v1[1:2], scale_v=v2[1:2], distribution = "gamma")$p.value, p_min)
  expect_less_than(ks.test(r_lba1$rt[r_lba1$response==1], normalised_n1CDF, A=A[1]+0.1, b=b[1], t0 = t0[1], shape_v=v1[1:2], scale_v=v2[1:2], distribution = "gamma")$p.value, p_min)
  expect_less_than(ks.test(r_lba2$rt[r_lba2$response==1], normalised_n1CDF, A=A[2], b=b[2], t0 = t0[2], shape_v=v1[3:4], scale_v=v2[3:4], st0 = 0, distribution = "gamma")$p.value, p_min)
  expect_less_than(ks.test(r_lba2$rt[r_lba2$response==1], normalised_n1CDF, A=A[2], b=b[2], t0 = t0[2], shape_v=v1[3:4], scale_v=v2[3:4], st0 = st0[1]+0.1, distribution = "gamma")$p.value, p_min)
  
  expect_more_than(ks.test(r_lba1$rt[r_lba1$response==1], normalised_n1CDF, A=A[1], b=b[1], t0 = t0[1], shape_v=v1[1:2], scale_v=v2[1:2], distribution = "gamma")$p.value, p_max)
  expect_more_than(ks.test(r_lba2$rt[r_lba2$response==1], normalised_n1CDF, A=A[2], b=b[2], t0 = t0[2], shape_v=v1[3:4], scale_v=v2[3:4], st0 = st0[1], distribution = "gamma")$p.value, p_max)
})

test_that("Frechet: n1CDF corresponds to random derivates", {
  normalised_n1CDF = function(t,...) n1CDF(t,...)/n1CDF(t=Inf,...) 
  samples <- 2e3
  p_min <- 0.001
  p_max <- 0.05
  A <- runif(2, 0.3, 0.9)
  b <- runif(2, 0.5, 1.5)
  t0 <- runif(2, 0.1, 0.7)
  v1 <- runif(4, 0.5, 1.5)
  v2 <- runif(4, 0.5, 1.5)
  st0 <- runif(1, 0.1, 0.5)
  r_lba1 <- rlba_frechet(samples, A=A[1], b=b[1], t0 = t0[1], shape_v=v1[1:2], scale_v=v2[1:2])
  r_lba2 <- rlba_frechet(samples, A=A[2], b=b[2], t0 = t0[2], shape_v=v1[3:4], scale_v=v2[3:4], st0 = st0[1])
  expect_less_than(ks.test(r_lba1$rt[r_lba1$response==1], "n1CDF",A=A[1], b=b[1], t0 = t0[1], shape_v=v1[1:2], scale_v=v2[1:2], distribution = "frechet")$p.value, p_min)
  expect_less_than(ks.test(r_lba1$rt[r_lba1$response==1], normalised_n1CDF, A=A[1]+0.2, b=b[1], t0 = t0[1], shape_v=v1[1:2], scale_v=v2[1:2], distribution = "frechet")$p.value, p_min)
  expect_less_than(ks.test(r_lba2$rt[r_lba2$response==1], normalised_n1CDF, A=A[2], b=b[2], t0 = t0[2], shape_v=v1[3:4], scale_v=v2[3:4], st0 = 0, distribution = "frechet")$p.value, p_min)
  expect_less_than(ks.test(r_lba2$rt[r_lba2$response==1], normalised_n1CDF, A=A[2], b=b[2], t0 = t0[2], shape_v=v1[3:4], scale_v=v2[3:4], st0 = st0[1]+0.1, distribution = "frechet")$p.value, p_min)
  
  expect_more_than(ks.test(r_lba1$rt[r_lba1$response==1], normalised_n1CDF, A=A[1], b=b[1], t0 = t0[1], shape_v=v1[1:2], scale_v=v2[1:2], distribution = "frechet")$p.value, p_max)
  expect_more_than(ks.test(r_lba2$rt[r_lba2$response==1], normalised_n1CDF, A=A[2], b=b[2], t0 = t0[2], shape_v=v1[3:4], scale_v=v2[3:4], st0 = st0[1], distribution = "frechet")$p.value, p_max)
})
