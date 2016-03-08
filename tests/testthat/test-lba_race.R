
context("LBA race functions: RNG is equivalent to n1")

x <- .Random.seed
set.seed(3)

tryCatch.W.E <- function(expr)
{
  mc <- match.call()
  mc2 <- match.call(definition = ks.test, call =  as.call(mc[[2]]))
  mc2[[1]] <- list
  
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e), warning = w.handler),warning = W, data = eval(mc2, envir = parent.frame()))
}

conditional_save_t <- function(t, distribution) {
  mc <- match.call()
  ex_data <- t$data
  #if (!is.null(t$warning)) save(ex_data, file = paste0(mc[[2]], "_", distribution, "_problem.Rdata"))
  #browser()
  #str(t)  
}

test_that("Norm: n1CDF corresponds to random derivates", {
  testthat::skip_on_cran()
  normalised_n1CDF = function(rt,...) n1CDF(rt,...)/n1CDF(rt=Inf,...) 
  samples <- 1e3
  p_min <- 0.001
  p_max <- 0.05
  A <- runif(2, 0.3, 0.9)
  b <- A+runif(2, 0, 0.5)
  t0 <- runif(2, 0.1, 0.7)
  v1 <- runif(4, 0.5, 1.5)
  v2 <- runif(4, 0.1, 0.5)
  st0 <- runif(1, 0.1, 0.5)
  r_lba1 <- rlba_norm(samples, A=A[1], b=b[1], t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2])
  r_lba2 <- rlba_norm(samples, A=A[2], b=b[2], t0 = t0[2], mean_v=v1[3:4], sd_v=v2[3:4], st0 = st0[1])
  t1 <- tryCatch.W.E(ks.test(r_lba1$rt[r_lba1$response==1], normalised_n1CDF, A=A[1], b=b[1]+0.1, t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2]))
  expect_less_than(t1$value$p.value, p_min)
  conditional_save_t(t1, "norm")
  
  t2 <- tryCatch.W.E(ks.test(r_lba2$rt[r_lba2$response==1], normalised_n1CDF, A=A[2], b=b[2], t0 = t0[2], mean_v=v1[3:4], sd_v=v2[3:4], st0 = 0))
  expect_less_than(t2$value$p.value, p_min)
  conditional_save_t(t2, "norm")
  
  t3 <- tryCatch.W.E(ks.test(r_lba2$rt[r_lba2$response==1], normalised_n1CDF, A=A[2], b=b[2], t0 = t0[2], mean_v=v1[3:4], sd_v=v2[3:4], st0 = st0[1]+0.1))
  expect_less_than(t3$value$p.value, p_min+0.003)
  conditional_save_t(t3, "norm")
  
  t4 <- tryCatch.W.E(ks.test(r_lba1$rt[r_lba1$response==1], normalised_n1CDF, A=A[1], b=b[1], t0 = t0[1], mean_v=v1[1:2], sd_v=v2[1:2]))
  expect_more_than(t4$value$p.value, p_max)
  conditional_save_t(t4, "norm")
  
  t5 <- tryCatch.W.E(ks.test(r_lba2$rt[r_lba2$response==1], normalised_n1CDF, A=A[2], b=b[2], t0 = t0[2], mean_v=v1[3:4], sd_v=v2[3:4], st0 = st0[1]))
  conditional_save_t(t5, "norm")
  
  expect_more_than(t5$value$p.value, p_max)
  
  #if (any(sapply(list(t1, t2, t3, t4, t5), function(x) !is.null(x$warning)))) browser()
})

test_that("Gamma: n1CDF corresponds to random derivates", {
  testthat::skip_on_cran()
  normalised_n1CDF = function(rt,...) n1CDF(rt,...)/n1CDF(rt=Inf,...) 
  samples <- 1e3
  p_min <- 0.01
  p_max <- 0.05
  A <- runif(2, 0.3, 0.9)
  b <- A+runif(2, 0, 0.5)
  t0 <- runif(2, 0.1, 0.7)
  v1 <- runif(4, 0.5, 1.5)
  v2 <- runif(4, 0.1, 0.5)
  st0 <- runif(1, 0.25, 0.75)
  r_lba1 <- rlba_gamma(samples, A=A[1], b=b[1], t0 = t0[1], shape_v=v1[1:2], scale_v=v2[1:2])
  r_lba2 <- rlba_gamma(samples, A=A[2], b=b[2], t0 = t0[2], shape_v=v1[3:4], scale_v=v2[3:4], st0 = st0[1])
  
  t1 <- tryCatch.W.E(ks.test(r_lba1$rt[r_lba1$response==1], normalised_n1CDF, A=A[1]+0.1, b=b[1]+0.2, t0 = t0[1], shape_v=v1[1:2], scale_v=v2[1:2], distribution = "gamma"))
  expect_less_than(t1$value$p.value, p_min)
  conditional_save_t(t1, "gamma")
  
  t2 <- tryCatch.W.E(ks.test(r_lba2$rt[r_lba2$response==1], normalised_n1CDF, A=A[2], b=b[2], t0 = t0[2], shape_v=v1[3:4], scale_v=v2[3:4], st0 = 0, distribution = "gamma"))
  expect_less_than(t2$value$p.value, p_min)
  conditional_save_t(t2, "gamma")
  
  t3 <- tryCatch.W.E(ks.test(r_lba2$rt[r_lba2$response==1], normalised_n1CDF, A=A[2], b=b[2], t0 = t0[2], shape_v=v1[3:4], scale_v=v2[3:4], st0 = st0[1]+0.5, distribution = "gamma"))
  expect_less_than(t3$value$p.value, p_min)
  conditional_save_t(t3, "gamma")
  
  t4 <- tryCatch.W.E(ks.test(r_lba1$rt[r_lba1$response==1], normalised_n1CDF, A=A[1], b=b[1], t0 = t0[1], shape_v=v1[1:2], scale_v=v2[1:2], distribution = "gamma"))
  expect_more_than(t4$value$p.value, p_max)
  conditional_save_t(t4, "gamma")
  
  t5 <- tryCatch.W.E(ks.test(r_lba2$rt[r_lba2$response==1], normalised_n1CDF, A=A[2], b=b[2], t0 = t0[2], shape_v=v1[3:4], scale_v=v2[3:4], st0 = st0[1], distribution = "gamma"))
  expect_more_than(t5$value$p.value, p_max)
  conditional_save_t(t5, "gamma")
  
  #if (any(sapply(list(t1, t2, t3, t4, t5), function(x) !is.null(x$warning)))) browser()
})

test_that("Frechet: n1CDF corresponds to random derivates", {
  testthat::skip_on_cran()
  normalised_n1CDF = function(rt,...) n1CDF(rt,...)/n1CDF(rt=Inf,...) 
  samples <- 2e2
  p_min <- 0.001
  p_max <- 0.05
  A <- runif(2, 0.3, 0.9)
  b <- A+runif(2, 0, 0.5)
  t0 <- runif(2, 0.1, 0.7)
  v1 <- runif(4, 0.5, 1.5)
  v2 <- runif(4, 0.5, 1.5)
  st0 <- runif(1, 0.25, 0.5)
  r_lba1 <- rlba_frechet(samples, A=A[1], b=b[1], t0 = t0[1], shape_v=v1[1:2], scale_v=v2[1:2])
  r_lba2 <- rlba_frechet(samples, A=A[2], b=b[2], t0 = t0[2], shape_v=v1[3:4], scale_v=v2[3:4], st0 = st0[1])

  t1 <- tryCatch.W.E(ks.test(r_lba1$rt[r_lba1$response==1], normalised_n1CDF, A=A[1]+0.4, b=b[1]+0.8, t0 = t0[1], shape_v=v1[1:2], scale_v=v2[1:2], distribution = "frechet"))
  expect_less_than(t1$value$p.value, p_min)
  conditional_save_t(t1, "frechet")
  
  t2 <- tryCatch.W.E(ks.test(r_lba2$rt[r_lba2$response==1], normalised_n1CDF, A=A[2], b=b[2], t0 = t0[2], shape_v=v1[3:4], scale_v=v2[3:4], st0 = 0, distribution = "frechet"))
  expect_less_than(t2$value$p.value, p_min)
  conditional_save_t(t2, "frechet")
  
  #browser()
  t3 <- tryCatch.W.E(ks.test(r_lba2$rt[r_lba2$response==1], normalised_n1CDF, A=A[2], b=b[2], t0 = t0[2], shape_v=v1[3:4], scale_v=v2[3:4], st0 = st0[1]+0.2, distribution = "frechet"))
  expect_less_than(t3$value$p.value, p_min)
  conditional_save_t(t3, "frechet")  
  
  t4 <- tryCatch.W.E(ks.test(r_lba1$rt[r_lba1$response==1], normalised_n1CDF, A=A[1], b=b[1], t0 = t0[1], shape_v=v1[1:2], scale_v=v2[1:2], distribution = "frechet"))
  expect_more_than(t4$value$p.value, p_max)
  conditional_save_t(t4, "frechet")
  
  t5 <- tryCatch.W.E(ks.test(r_lba2$rt[r_lba2$response==1], normalised_n1CDF, A=A[2], b=b[2], t0 = t0[2], shape_v=v1[3:4], scale_v=v2[3:4], st0 = st0[1], distribution = "frechet"))
  #t5 <- tryCatch.W.E(ks.test(r_lba2$rt[r_lba2$response==1]-t0[2], normalised_n1CDF, A=A[2], b=b[2], t0 = 0, shape_v=v1[3:4], scale_v=v2[3:4], st0 = st0[1], distribution = "frechet"))
  expect_more_than(t5$value$p.value, p_max)
  conditional_save_t(t5, "frechet")
  
  #if (any(sapply(list(t1, t2, t3, t4, t5), function(x) !is.null(x$warning)))) browser()
})

test_that("lnorm: n1CDF corresponds to random derivates", {
  testthat::skip_on_cran()
  normalised_n1CDF = function(rt,...) n1CDF(rt,...)/n1CDF(rt=Inf,...) 
  samples <- 1e3
  p_min <- 0.0001
  p_max <- 0.05
  A <- runif(2, 0.3, 0.9)
  b <- A+runif(2, 0, 0.5)
  t0 <- runif(2, 0.1, 0.7)
  v1 <- runif(4, 0.5, 1.5)
  v2 <- runif(4, 0.1, 0.5)
  st0 <- runif(1, 0.1, 0.5)
  r_lba1 <- rlba_lnorm(samples, A=A[1], b=b[1], t0 = t0[1], meanlog_v=v1[1:2], sdlog_v=v2[1:2])
  r_lba2 <- rlba_lnorm(samples, A=A[2], b=b[2], t0 = t0[2], meanlog_v=v1[3:4], sdlog_v=v2[3:4], st0 = st0[1])
  
  t1 <- tryCatch.W.E(ks.test(r_lba1$rt[r_lba1$response==1], normalised_n1CDF, A=A[1], b=b[1]+0.1, t0 = t0[1], meanlog_v=v1[1:2], sdlog_v=v2[1:2], distribution = "lnorm"))
  expect_less_than(t1$value$p.value, p_min)
  conditional_save_t(t1, "lnorm")
  
  t2 <- tryCatch.W.E(ks.test(r_lba2$rt[r_lba2$response==1], normalised_n1CDF, A=A[2], b=b[2], t0 = t0[2], meanlog_v=v1[3:4], sdlog_v=v2[3:4], st0 = 0, distribution = "lnorm"))
  expect_less_than(t2$value$p.value, p_min)
  conditional_save_t(t2, "lnorm")
  
  t3 <- tryCatch.W.E(ks.test(r_lba2$rt[r_lba2$response==1], normalised_n1CDF, A=A[2], b=b[2], t0 = t0[2], meanlog_v=v1[3:4], sdlog_v=v2[3:4], st0 = st0[1]+0.2, distribution = "lnorm"))
  expect_less_than(t3$value$p.value, p_min)
  conditional_save_t(t3, "lnorm")
  
  #t4 <- tryCatch.W.E(ks.test(r_lba1$rt[r_lba1$response==1], normalised_n1CDF, A=A[1], b=b[1], t0 = t0[1], meanlog_v=v1[1:2], sdlog_v=v2[1:2], distribution = "lnorm"))
  t4 <- tryCatch.W.E(ks.test(r_lba1$rt[r_lba1$response==1]-t0[1], normalised_n1CDF, A=A[1], b=b[1], t0 = 0, meanlog_v=v1[1:2], sdlog_v=v2[1:2], distribution = "lnorm"))
  expect_more_than(t4$value$p.value, p_max)
  conditional_save_t(t4, "lnorm")
  
  #t5 <- tryCatch.W.E(ks.test(r_lba2$rt[r_lba2$response==1], normalised_n1CDF, A=A[2], b=b[2], t0 = t0[2], meanlog_v=v1[3:4], sdlog_v=v2[3:4], st0 = st0[1], distribution = "lnorm"))
  t5 <- tryCatch.W.E(ks.test(pmax(r_lba2$rt[r_lba2$response==1]-t0[2],0), normalised_n1CDF, A=A[2], b=b[2], t0 = 0, meanlog_v=v1[3:4], sdlog_v=v2[3:4], st0 = st0[1], distribution = "lnorm"))
  expect_more_than(t5$value$p.value, p_max)
  conditional_save_t(t5, "lnorm")
  
  #if (any(sapply(list(t1, t2, t3, t4, t5), function(x) !is.null(x$warning)))) browser()
})


test_that("Norm: n1CDF corresponds to random derivates with accumulatorwise parameters", {
  testthat::skip_on_cran()
  normalised_n1CDF = function(rt,...) n1CDF(rt,...)/n1CDF(rt=Inf,...) 
  samples <- 1e3
  p_min <- 0.001
  p_max <- 0.05
  A <- runif(4, 0.3, 0.9)
  b <- A+runif(4, 0, 0.5)
  t0 <- runif(4, 0.1, 0.7)
  v1 <- runif(4, 0.5, 1.5)
  v2 <- runif(4, 0.1, 0.5)
  st0 <- runif(1, 0.1, 0.5)
  r_lba1 <- rLBA(samples, A=list(A[1], A[2]), b=list(b[1], b[2]), t0 = list(t0[1], t0[2]), mean_v=v1[1:2], sd_v=v2[1:2])
  r_lba2 <- rLBA(samples, A=list(A[3], A[4]), b=list(b[3], b[4]), t0 = list(t0[3], t0[4]), mean_v=v1[3:4], sd_v=v2[3:4], st0 = st0[1])
  t1 <- tryCatch.W.E(ks.test(r_lba1$rt[r_lba1$response==1], normalised_n1CDF, A=list(A[1], A[2]), b=list(b[1], b[2]), t0 = list(t0[1]+0.1, t0[2]+0.1), mean_v=v1[1:2], sd_v=v2[1:2]))
  expect_less_than(t1$value$p.value, p_min)
  conditional_save_t(t1, "norm")
  
  t2 <- tryCatch.W.E(ks.test(r_lba2$rt[r_lba2$response==1], normalised_n1CDF, A=list(A[3], A[4]), b=list(b[3], b[4]), t0 = list(t0[3], t0[4]), mean_v=v1[3:4], sd_v=v2[3:4], st0 = 0))
  expect_less_than(t2$value$p.value, p_min)
  conditional_save_t(t2, "norm")
  
  t3 <- tryCatch.W.E(ks.test(r_lba2$rt[r_lba2$response==1], normalised_n1CDF, A=list(A[3], A[4]), b=list(b[3], b[4]), t0 = list(t0[3], t0[4]), mean_v=v1[3:4], sd_v=v2[3:4], st0 = st0[1]+0.1))
  expect_less_than(t3$value$p.value, p_min) #+0.003
  conditional_save_t(t3, "norm")
  
  t4 <- tryCatch.W.E(ks.test(r_lba1$rt[r_lba1$response==1], normalised_n1CDF, A=list(A[1], A[2]), b=list(b[1], b[2]), t0 = list(t0[1], t0[2]), mean_v=v1[1:2], sd_v=v2[1:2]))
  expect_more_than(t4$value$p.value, p_max)
  conditional_save_t(t4, "norm")
  
  t5 <- tryCatch.W.E(ks.test(r_lba2$rt[r_lba2$response==1], normalised_n1CDF, A=list(A[3], A[4]), b=list(b[3], b[4]), t0 = list(t0[3], t0[4]), mean_v=v1[3:4], sd_v=v2[3:4], st0 = st0[1]))
  conditional_save_t(t5, "norm")
  
  expect_more_than(t5$value$p.value, p_max)
  
  #if (any(sapply(list(t1, t2, t3, t4, t5), function(x) !is.null(x$warning)))) browser()
})

test_that("rLBA works with trial wise input", {
  set.seed(1)
  rt1 <- rLBA(100, A=rep(c(0.5, 0.4), each = 50), b=rep(c(1, 1.4), each = 50), t0 = rep(c(0.2, 0.4), each = 50), mean_v=list(rep(c(1.5, 0.9), each = 50), rep(c(0.9, 1.5), each = 50)), sd_v=list(rep(c(0.2, 0.4), each = 50), rep(c(0.3, 0.5), each = 50)), st0 = rep(c(0.1, 0.2), each = 50))
  set.seed(1)
  rt2a <- rLBA(50, A=0.5, b=1, t0 = 0.2, mean_v=c(1.5, 0.9), sd_v=list(0.2, 0.3), st0 = 0.1)
  rt2b <- rLBA(50, A=0.4, b=1.4, t0 = 0.4, mean_v=c(0.9, 1.5), sd_v=c(0.4, 0.5), st0 = 0.2)
  expect_identical(rt1, rbind(rt2a, rt2b))
  
})

.Random.seed <<- x