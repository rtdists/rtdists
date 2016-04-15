
context("Diffusion Model: quantile bugs")

test_that("qdiffusion does not fail for certain values", {
  qex <- c(a = 1.97092532512193, t0 = 0.22083875893466, sv = 0.462630056494015, v = -2.59881372383245, resp_prop = 0.99877722499984, z = 0.5, sz= 0.1)
  
  #pdiffusion(0.28165, boundary = "lower", a=qex["a"], v=qex["v"], t0 = qex["t0"], sv = qex["sv"], sz = 0.1, z = 0.5)
  
  expect_true(!is.na(
    qdiffusion(p = 0.0998, boundary = "lower", a=qex["a"], v=qex["v"], t0 = qex["t0"], sv = qex["sv"], sz = 0.1, z = 0.5)
  ))

})
