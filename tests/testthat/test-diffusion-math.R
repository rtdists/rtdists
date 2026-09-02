
context("Diffusion Model: Compare with RWiener")

test_that("ddiffusion is equal to dwiener", {
  if (require(RWiener)) {
    for (a in seq(0.5, 2.0, length.out = 6)) {
      for (v in seq(0.5, 2.0, length.out = 6)) {
        for (t0 in seq(0.05, 0.5, length.out = 6)) {
          for (z in seq(0.01, 0.99, length.out = 9)) {
            expect_equal(
              ddiffusion(seq(0, 3, length.out = 15), a=a, v=v, t0=t0, z = z*a)
              ,
              dwiener(seq(0, 3, length.out = 15), resp = rep("upper", 15), alpha=a, delta=v, tau = t0, beta = z)
            )
            expect_equal(
              ddiffusion(seq(0, 3, length.out = 16), c("upper", "lower"), a=a, v=v, t0=t0, z = z*a)
              ,
              dwiener(seq(0, 3, length.out = 16), resp = rep(c("upper", "lower"), 8), alpha=a, delta=v, tau = t0, beta = z)
              , tolerance = 0.00001
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
    for (a in seq(0.5, 2.0, length.out = 6)) {
      for (v in seq(0.5, 2.0, length.out = 6)) {
        for (t0 in seq(0.05, 0.5, length.out = 6)) {
          for (z in seq(0.01, 0.99, length.out = 10)) {
            expect_equal(
              pdiffusion(seq(0, 3, length.out = 15), a=a, v=v, t0=t0, z = z*a)
              ,
              pwiener(seq(0, 3, length.out = 15), resp = rep("upper", 15), alpha=a, delta=v, tau = t0, beta = z)
            , tolerance = 0.01)
            expect_equivalent(
              pdiffusion(seq(0, 3, length.out = 16), c("upper", "lower"), a=a, v=v, t0=t0, z = z*a)
              ,
              pwiener(seq(0, 3, length.out = 16), resp = rep(c("upper", "lower"), 8), alpha=a, delta=v, tau = t0, beta = z)
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

## ---- WienR comparison tests ----
## Parameter mapping: rtdists -> WienR
##   z (absolute) -> w = z/a (relative)
##   sz (absolute) -> sw = sz/a (relative)
##   t0, st0, sv, a, v map directly
##   rtdists s defaults to 1; WienR has no s parameter (fixed to 1)
## requires precision > 3 to work for all examples, see also: https://github.com/rtdists/rtdists/issues/28

test_that("ddiffusion with variability parameters is equal to WienerPDF", {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  if (require(WienR)) {
    gridsize <- 3
    precision <- 6
    rt <- seq(0.3, 3, length.out = 20)
    for (a in seq(0.5, 2.0, length.out = gridsize)) {
      for (v in seq(0.5, 2.0, length.out = gridsize)) {
        for (t0 in seq(0.05, 0.5, length.out = gridsize)) {
          for (z in seq(0.3, 0.7, length.out = gridsize)) {
            z_abs <- z * a
            for (sv in seq(0, 0.5, length.out = gridsize)) {
              for (sz in seq(0, 0.25, length.out = gridsize)) {
                # skip if starting point range extends beyond [0, a]
                if (z_abs - sz/2 <= 0 || z_abs + sz/2 >= a) next
                for (st0 in seq(0, 0.2, length.out = gridsize)) {
                  expect_equal(
                    ddiffusion(rt, a = a, v = v, t0 = t0, z = z_abs, sv = sv, sz = sz, st0 = st0, precision = precision)
                    ,
                    WienerPDF(t = rt, response = rep("upper", length(rt)),
                              a = a, v = v, w = z, t0 = t0, sv = sv, sw = sz / a, st0 = st0)$value
                  , tolerance = 0.001)

                  # all variability, lower boundary
                  expect_equal(
                    ddiffusion(rt, "lower", a = a, v = v, t0 = t0, z = z_abs, sv = sv, sz = sz, st0 = st0, precision = precision)
                    ,
                    WienerPDF(t = rt, response = rep("lower", length(rt)),
                              a = a, v = v, w = z, t0 = t0, sv = sv, sw = sz / a, st0 = st0)$value
                  , tolerance = 0.001)
                }
              }
            }
          }
        }
      }
    }
  }
})

test_that("pdiffusion with variability parameters is equal to WienerPDF", {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  if (require(WienR)) {
    gridsize <- 3
    precision <- 3
    rt <- seq(0.3, 3, length.out = 20)
    for (a in seq(0.5, 2.0, length.out = gridsize)) {
      for (v in seq(0.5, 2.0, length.out = gridsize)) {
        for (t0 in seq(0.05, 0.5, length.out = gridsize)) {
          for (z in seq(0.3, 0.7, length.out = gridsize)) {
            z_abs <- z * a
            for (sv in seq(0, 0.5, length.out = gridsize)) {
              for (sz in seq(0, 0.25, length.out = gridsize)) {
                # skip if starting point range extends beyond [0, a]
                if (z_abs - sz/2 <= 0 || z_abs + sz/2 >= a) next
                for (st0 in seq(0, 0.2, length.out = gridsize)) {
                  expect_equal(
                    pdiffusion(rt, a = a, v = v, t0 = t0, z = z_abs, sv = sv, sz = sz, st0 = st0, precision = precision)
                    ,
                    WienerCDF(t = rt, response = rep("upper", length(rt)),
                              a = a, v = v, w = z, t0 = t0, sv = sv, sw = sz / a, st0 = st0)$value
                  , tolerance = 1e-3)

                  # all variability, lower boundary
                  expect_equal(
                    pdiffusion(rt, "lower", a = a, v = v, t0 = t0, z = z_abs, sv = sv, sz = sz, st0 = st0, precision = precision)
                    ,
                    WienerCDF(t = rt, response = rep("lower", length(rt)),
                              a = a, v = v, w = z, t0 = t0, sv = sv, sw = sz / a, st0 = st0)$value
                  , tolerance = 1e-3)

                }
              }
            }
          }
        }
      }
    }
  }
})
