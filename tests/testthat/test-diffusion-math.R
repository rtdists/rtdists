
context("Diffusion Model: Compare with RWiener")

test_that("ddiffusion is equal to dwiener", {
  if (require(RWiener)) {
    for (a in seq(0.5, 2.0, length.out = 6)) {
      for (v in seq(0.5, 2.0, length.out = 6)) {
        for (t0 in seq(0.05, 0.5, length.out = 6)) {
          for (z in seq(0.4, 0.6, length.out = 6)) {
            expect_equivalent(
              ddiffusion(seq(0, 3, length.out = 15), a=a, v=v, t0=t0, z = z*a)
              ,
              dwiener(seq(0, 3, length.out = 15), resp = rep("upper", 15), alpha=a, delta=v, tau = t0, beta = z)
            )
            expect_equivalent(
              ddiffusion(seq(0, 3, length.out = 16), c("upper", "lower"), a=a, v=v, t0=t0, z = z*a)
              ,
              dwiener(seq(0, 3, length.out = 16), resp = rep(c("upper", "lower"), 8), alpha=a, delta=v, tau = t0, beta = z)
            )
          }
        }
      }
    }
  }
  
})

test_that("pdiffusion is equal to pwiener", {
  testthat::skip_on_cran()
  if (require(RWiener)) {
    for (a in seq(0.5, 2.0, length.out = 10)) {
      for (v in seq(0.5, 2.0, length.out = 10)) {
        for (t0 in seq(0.05, 0.5, length.out = 10)) {
          for (z in seq(0.4, 0.6, length.out = 7)) {
            expect_equal(
              pdiffusion(seq(0, 3, length.out = 15), a=a, v=v, t0=t0, z = z*a)
              ,
              pwiener(seq(0, 3, length.out = 15), resp = rep("upper", 15), alpha=a, delta=v, tau = t0, beta = z)
            , tolerance = 0.01)
          }
        }
      }
    }
  }
})


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


test_that("Norm: pdiffusion corresponds to random derivates", {
  testthat::skip_on_cran()
  #testthat::skip_on_travis()
  normalised_pdiffusion <- function(rt,...) pdiffusion(rt,...)/pdiffusion(rt=10, ...) 
  normalised_pwiener <- function(q,...) pwiener(q,  resp = rep("upper", length(q)), ...)/pwiener(q=10, resp = "upper", ...)
  samples <- 1e4
  p_min <- 0.001
  p_max <- 0.01
  a <- runif(1, 0.3, 0.9)
  t0 <- runif(1, 0.1, 0.5)
  v <- runif(1, 0.5, 2.5)
  z <- runif(1, 0.5, 0.6)
  r_diffusion <- rdiffusion(samples, a=a, t0=t0, v=v, z=z*a)
  t1 <- tryCatch.W.E(ks.test(r_diffusion$rt[r_diffusion$response=="upper"], normalised_pdiffusion, a=a*2, t0=t0, v=v*2, z=z*a))
  expect_lt(t1$value$p.value, p_min)
  
  t2 <- tryCatch.W.E(ks.test(r_diffusion$rt[r_diffusion$response=="upper"], normalised_pdiffusion, a=a, t0=t0, v=v, z=z*a))
  expect_gt(t2$value$p.value, p_max)
  
  t3 <- tryCatch.W.E(ks.test(r_diffusion$rt[r_diffusion$response=="upper"], normalised_pwiener, alpha=a, delta=v, tau = t0, beta = z))
  expect_gt(t3$value$p.value, p_max)
  
})

test_that("Norm: pdiffusion corresponds to random derivates (with variabilities)", {
  testthat::skip_on_cran()
  #testthat::skip_on_travis()
  normalised_pdiffusion <- function(rt,...) pdiffusion(rt,...)/pdiffusion(rt=10, ...) 
  samples <- 1e4
  p_min <- 0.001
  p_max <- 0.01
  a <- runif(1, 0.3, 0.9)
  t0 <- runif(1, 0.1, 0.5)
  v <- runif(1, 0.5, 2.5)
  sv <- runif(1, 0.1, 0.5)
  sz <- runif(1, 0.05, 0.2)
  z <- runif(1, 0.5, 0.6)
  r_diffusion <- rdiffusion(samples, a=a, t0=t0, v=v, z=z*a, sz=sz, sv = sv)
  t1 <- tryCatch.W.E(ks.test(r_diffusion$rt[r_diffusion$response=="upper"], normalised_pdiffusion, a=a, t0=t0, v=v, z=z*a, sv=1, sz = 0.6*a))
  expect_lt(t1$value$p.value, p_min)
  
  t2 <- tryCatch.W.E(ks.test(r_diffusion$rt[r_diffusion$response=="upper"], normalised_pdiffusion, a=a, t0=t0, v=v, z=z*a,sv=sv, sz=sz))
  expect_gt(t2$value$p.value, p_max)
  
})