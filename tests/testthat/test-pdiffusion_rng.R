context("pdiffusion functions: RNG is equivalent to pdiffusion")

#x <- .Random.seed
#set.seed(3)

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

test_that("Norm: pdiffusion corresponds to random derivates with specific values", {
  testthat::skip_on_cran()
  #testthat::skip_on_travis()
  normalised_pdiffusion = function(rt,...) pdiffusion(rt,...)/pdiffusion(rt=Inf,...) 
  samples <- 1e3
  p_min <- 0.01
  p_max <- 0.01
  diffusion_pars <- structure(list(par = structure(c(1.32060063610882, 3.27271614698074, 0.338560144920614, 0.34996447540773, 0.201794924457386, 1.05516829794661), .Names = c("a", "v", "t0", "sz", "st0", "sv"))))
  
  r_diff1 <- rdiffusion(samples, a=diffusion_pars$par["a"], v=diffusion_pars$par["v"], t0=diffusion_pars$par["t0"], sz=diffusion_pars$par["sz"], st0=diffusion_pars$par["st0"], sv=diffusion_pars$par["sv"])
  t1 <- tryCatch.W.E(ks.test(r_diff1$rt[r_diff1$response=="upper"], normalised_pdiffusion, a=1.1, v=diffusion_pars$par["v"], t0=diffusion_pars$par["t0"], sz=diffusion_pars$par["sz"], st0=diffusion_pars$par["st0"], sv=diffusion_pars$par["sv"]))
  expect_lt(t1$value$p.value, p_min)
  conditional_save_t(t1, "norm")
  
  t2 <- tryCatch.W.E(ks.test(r_diff1$rt[r_diff1$response=="upper"], normalised_pdiffusion, a=diffusion_pars$par["a"], v=diffusion_pars$par["v"]-0.5, t0=diffusion_pars$par["t0"], sz=diffusion_pars$par["sz"], st0=diffusion_pars$par["st0"], sv=diffusion_pars$par["sv"]))
  expect_lt(t2$value$p.value, p_min)
  conditional_save_t(t2, "norm")
  
  t3 <- tryCatch.W.E(ks.test(r_diff1$rt[r_diff1$response=="upper"], normalised_pdiffusion, a=diffusion_pars$par["a"], v=diffusion_pars$par["v"]-0.5, t0=diffusion_pars$par["t0"], sz=diffusion_pars$par["sz"], st0=diffusion_pars$par["st0"], sv=diffusion_pars$par["sv"]))
  expect_lt(t3$value$p.value, p_min)
  conditional_save_t(t3, "norm")
  
  t4 <- tryCatch.W.E(ks.test(r_diff1$rt[r_diff1$response=="upper"], normalised_pdiffusion, a=diffusion_pars$par["a"], v=diffusion_pars$par["v"], t0=diffusion_pars$par["t0"]-0.1, sz=diffusion_pars$par["sz"], st0=diffusion_pars$par["st0"], sv=diffusion_pars$par["sv"]))
  expect_lt(t4$value$p.value, p_min)
  conditional_save_t(t4, "norm")
  
  t5 <- tryCatch.W.E(ks.test(r_diff1$rt[r_diff1$response=="upper"], normalised_pdiffusion, a=diffusion_pars$par["a"], v=diffusion_pars$par["v"], t0=diffusion_pars$par["t0"], sz=diffusion_pars$par["sz"], st0=diffusion_pars$par["st0"], sv=diffusion_pars$par["sv"]))
  conditional_save_t(t5, "norm")
  
  expect_gt(t5$value$p.value, p_max)
  
  #if (any(sapply(list(t1, t2, t3, t4, t5), function(x) !is.null(x$warning)))) browser()
})



test_that("Norm: pdiffusion corresponds to random derivates with random values", {
  testthat::skip_on_cran()
  #testthat::skip_on_travis()
  normalised_pdiffusion = function(rt,...) pdiffusion(rt,...)/pdiffusion(rt=Inf,...) 
  samples <- 5e3
  p_min <- 0.01
  p_max <- 0.01
  diffusion_pars <- list(
    a = runif(1, 0.5, 1.5),
    v = runif(1, 2, 3.5),
    t0 = runif(1, 0.2, 0.4),
    sz = runif(1, 0.1, 0.2),
    st0 = runif(1, 0.1, 0.2),
    sv = runif(1, 0.5, 1.5)
  )
    
  r_diff1 <- rdiffusion(samples, a=diffusion_pars[["a"]], v=diffusion_pars[["v"]], t0=diffusion_pars[["t0"]], sz=diffusion_pars[["sz"]], st0=diffusion_pars[["st0"]], sv=diffusion_pars[["sv"]])
  t1 <- tryCatch.W.E(ks.test(r_diff1$rt[r_diff1$response=="upper"], normalised_pdiffusion, a=diffusion_pars[["a"]]-0.3, v=diffusion_pars[["v"]], t0=diffusion_pars[["t0"]], sz=diffusion_pars[["sz"]], st0=diffusion_pars[["st0"]], sv=diffusion_pars[["sv"]]))
  expect_lt(t1$value$p.value, p_min)
  conditional_save_t(t1, "norm")
  
  t2 <- tryCatch.W.E(ks.test(r_diff1$rt[r_diff1$response=="upper"], normalised_pdiffusion, a=diffusion_pars[["a"]], v=pmin(0, diffusion_pars[["v"]]-1), t0=diffusion_pars[["t0"]], sz=diffusion_pars[["sz"]], st0=diffusion_pars[["st0"]], sv=diffusion_pars[["sv"]]))
  expect_lt(t2$value$p.value, p_min)
  conditional_save_t(t2, "norm")
  
  t3 <- tryCatch.W.E(ks.test(r_diff1$rt[r_diff1$response=="upper"], normalised_pdiffusion, a=diffusion_pars[["a"]], v=diffusion_pars[["v"]], t0=diffusion_pars[["t0"]], sz=diffusion_pars[["sz"]]+0.25, st0=0, sv=diffusion_pars[["sv"]]-0.5))
  expect_lt(t3$value$p.value, p_min)
  conditional_save_t(t3, "norm")
  
  t4 <- tryCatch.W.E(ks.test(r_diff1$rt[r_diff1$response=="upper"], normalised_pdiffusion, a=diffusion_pars[["a"]], v=diffusion_pars[["v"]], t0=diffusion_pars[["t0"]]-0.08, sz=diffusion_pars[["sz"]], st0=diffusion_pars[["st0"]], sv=diffusion_pars[["sv"]]))
  expect_lt(t4$value$p.value, p_max)
  conditional_save_t(t4, "norm")
  
  t5 <- tryCatch.W.E(ks.test(r_diff1$rt[r_diff1$response=="upper"], normalised_pdiffusion, a=diffusion_pars[["a"]], v=diffusion_pars[["v"]], t0=diffusion_pars[["t0"]], sz=diffusion_pars[["sz"]], st0=diffusion_pars[["st0"]], sv=diffusion_pars[["sv"]]))
  conditional_save_t(t5, "norm")
  expect_gt(t5$value$p.value, p_max)
  
  #if (any(sapply(list(t1, t2, t3, t4, t5), function(x) !is.null(x$warning)))) browser()
})


test_that("Norm: pdiffusion corresponds to random derivates with random values 2", {
  testthat::skip_on_cran()
  #testthat::skip_on_travis()
  normalised_pdiffusion = function(rt,...) pdiffusion(rt,...)/pdiffusion(rt=Inf,...) 
  samples <- 5e3
  p_min <- 0.01
  p_max <- 0.01
  diffusion_pars <- list(
    a = runif(1, 0.5, 1.5),
    v = runif(1, 2, 3.5),
    t0 = runif(1, 0.2, 0.4),
    sz = runif(1, 0, 0.4),
    st0 = runif(1, 0, 0.2),
    sv = runif(1, 0, 1)
  )
    
  r_diff1 <- rdiffusion(samples, a=diffusion_pars[["a"]], v=diffusion_pars[["v"]], t0=diffusion_pars[["t0"]], sz=diffusion_pars[["sz"]], st0=diffusion_pars[["st0"]], sv=diffusion_pars[["sv"]])
  t1 <- tryCatch.W.E(ks.test(r_diff1$rt[r_diff1$response=="upper"], normalised_pdiffusion, a=diffusion_pars[["a"]]-0.3, v=diffusion_pars[["v"]], t0=diffusion_pars[["t0"]], sz=diffusion_pars[["sz"]], st0=diffusion_pars[["st0"]], sv=diffusion_pars[["sv"]]))
  expect_lt(t1$value$p.value, p_min)
  conditional_save_t(t1, "norm")
  
  t2 <- tryCatch.W.E(ks.test(r_diff1$rt[r_diff1$response=="upper"], normalised_pdiffusion, a=diffusion_pars[["a"]], v=diffusion_pars[["v"]]+4, t0=diffusion_pars[["t0"]], sz=diffusion_pars[["sz"]], st0=diffusion_pars[["st0"]], sv=diffusion_pars[["sv"]]))
  expect_lt(t2$value$p.value, p_min)
  conditional_save_t(t2, "norm")
  
  t3 <- tryCatch.W.E(ks.test(r_diff1$rt[r_diff1$response=="upper"], normalised_pdiffusion, a=diffusion_pars[["a"]], v=diffusion_pars[["v"]], t0=diffusion_pars[["t0"]], sz=diffusion_pars[["sz"]]+0.25, st0=0, sv=0))
  expect_lt(t3$value$p.value, p_min)
  conditional_save_t(t3, "norm")
  
  t4 <- tryCatch.W.E(ks.test(r_diff1$rt[r_diff1$response=="upper"], normalised_pdiffusion, a=diffusion_pars[["a"]], v=diffusion_pars[["v"]], t0=diffusion_pars[["t0"]]-0.08, sz=diffusion_pars[["sz"]], st0=diffusion_pars[["st0"]], sv=diffusion_pars[["sv"]]))
  expect_lt(t4$value$p.value, p_max)
  conditional_save_t(t4, "norm")
  
  t5 <- tryCatch.W.E(ks.test(r_diff1$rt[r_diff1$response=="upper"], normalised_pdiffusion, a=diffusion_pars[["a"]], v=diffusion_pars[["v"]], t0=diffusion_pars[["t0"]], sz=diffusion_pars[["sz"]], st0=diffusion_pars[["st0"]], sv=diffusion_pars[["sv"]]))
  conditional_save_t(t5, "norm")
  expect_gt(t5$value$p.value, p_max)
  
  #if (any(sapply(list(t1, t2, t3, t4, t5), function(x) !is.null(x$warning)))) browser()
})

