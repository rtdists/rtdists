
test_that("rdiffusion does not return rt > t0", {
  set.seed(18433)
  df1 <- rdiffusion(5, a=2,st0=0, sv=0, v = 1, t0 = 0.4)
  expect_false(any(df1$rt <= 0.4))
})
