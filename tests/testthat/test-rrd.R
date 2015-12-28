
context("diffusion parameter input (via rdiffusion)")

test_that("check individual parameters:", {
  expect_that(rdiffusion(10, a=1, z=0.5, v=2, t0=0.5, d=0, sz = 0, sv = 0, st0 = 0), is_a("data.frame"))
  expect_that(suppressWarnings(rdiffusion(10, a=1, z=0.5, v=2, t0=0.5, d=0, sz = 0, sv = NULL, st0 = 0)), throws_error("Not enough parameters"))
  expect_that(rdiffusion(10, a=1, z=0.5, v=2, t0=0.5, d=0, sz = 0, sv = Inf, st0 = 0), throws_error())
  expect_that(suppressWarnings(rdiffusion(10, a=1, z=NA, v=2, t0=0.5,  d=0, sz = 0, sv = 0, st0 = 0)), throws_error())  
})

# test_that("check parameters:", {
#   p1 <- c(1, 0.5, 2, 0.5, 0, 0, 0, 0)
#   expect_that(rdiffusion(10, parameters = p1), is_a("data.frame"))
#   expect_that(rdiffusion(10, parameters = p1[1:7]), throws_error())
#   names(p1) <- c("a", "z", "v","t0","d", "sz","sv","st0")
#   expect_that(rdiffusion(10, parameters = p1), is_a("data.frame"))
#   names(p1) <- c(c("a","v","t0","z"), sample(c("sz","sv","st0", "d")))
#   expect_that(rdiffusion(10, parameters = p1), is_a("data.frame"))
#   names(p1)[3] <- "xx"
#   expect_that(rdiffusion(10, parameters = p1), throws_error())
#   names(p1) <- NULL
#   p1[1] <- NA
#   expect_that(rdiffusion(10, parameters = p1), throws_error())
#   p1[1] <- Inf
#   expect_that(rdiffusion(10, parameters = p1), throws_error())
# })


context("rdiffusion: random number generation for the diffusion model")
test_that("rdiffusion works", {
  rrds1 <- rdiffusion(10, a=1, z=0.5, v=2, t0=0.5, d=0, sz = 0, sv = 0, st0 = 0)
  rrds2 <- rdiffusion(10, a=1, z=0.5, v=2, t0=0.5, d=0, sz = 0, sv = 0, st0 = 0)
  expect_that(rrds1, is_a("data.frame"))
  expect_that(isTRUE(all.equal(rrds1, rrds2)), is_false())
  
  set.seed(1)
  rrds1 <- rdiffusion(10, a=1, z=0.5, v=2, t0=0.5, d=0, sz = 0, sv = 0, st0 = 0)
  set.seed(1)
  rrds2 <- rdiffusion(10, a=1, z=0.5, v=2, t0=0.5, d=0, sz = 0, sv = 0, st0 = 0)
  expect_that(rrds1, equals(rrds2))
  set.seed(NULL)
})
