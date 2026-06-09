
test_that("rdiffusion does not return rt > t0", {
  set.seed(678)
  df1 <- rdiffusion(1e6, a=2,st0=0, sv=1, v = -1, t0 = 0.4)
  expect_false(any(df1$rt <= 0.4))
})
