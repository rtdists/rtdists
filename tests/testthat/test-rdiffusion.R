
test_that("rdiffusion does not return rt > t0", {
  set.seed(18433)
  df1 <- rdiffusion(5, a=2,st0=0, sv=0, v = 1, t0 = 0.4)
  expect_false(any(df1$rt <= 0.4))
})

test_that("length(n) > 1 is taken as the number of observations (#6)", {
  expect_equal(nrow(rdiffusion(rep(500, 2), a = 1, v = 2, t0 = 0.5)), 2)
  expect_equal(nrow(rLBA(1:3, A = 0.5, b = 1, t0 = 0.5, mean_v = c(2.4, 1.6), sd_v = c(1, 1.2), silent = TRUE)), 3)
  expect_equal(nrow(rRDM(1:4, A = 0.5, b = 1, t0 = 0.5, v = c(2, 1), silent = TRUE)), 4)
})
