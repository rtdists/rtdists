context("Diffusion pdiffusion and ddiffusion functions.")
tolerance = 0.0001)

test_that("diffusion functions work with numeric and factor boundaries", {
  n_test <- 20
  rts <- rdiffusion(n_test, a=1, z=0.5, v=2, t0=0.5, d=0, sz = 0, sv = 0, st0 = 0)
  expect_is(ddiffusion(rts$rt, response = rts$response, a=1, z=0.5, v=2, t0=0.5, d=0, sz = 0, sv = 0, st0 = 0), "numeric")
  expect_is(pdiffusion(sort(rts$rt), response = rts$response, a=1, z=0.5, v=2, t0=0.5, d=0, sz = 0, sv = 0, st0 = 0), "numeric")
  expect_is(ddiffusion(rts$rt, response = sample(1:2, 20, replace = TRUE), a=1, z=0.5, v=2, t0=0.5, d=0, sz = 0, sv = 0, st0 = 0), "numeric")
  expect_is(pdiffusion(sort(rts$rt), response = sample(1:2, 20, replace = TRUE), a=1, z=0.5, v=2, t0=0.5, d=0, sz = 0, sv = 0, st0 = 0), "numeric")
  expect_error(ddiffusion(rts$rt, rep_len(1:3, length.out=20), a=1, z=0.5, v=2, t0=0.5, d=0, sz = 0, sv = 0, st0 = 0), "response")
  expect_error(pdiffusion(sort(rts$rt), rep_len(1:3, length.out=20), a=1, z=0.5, v=2, t0=0.5, d=0, sz = 0, sv = 0, st0 = 0), "response")
})

test_that("diffusion functions are identical with all input options", {
  rt1 <- rdiffusion(500, a=1, v=2, t0=0.5)
  # get density for random RTs:
  ref <- sum(log(ddiffusion(rt1$rt, rt1$response, a=1, v=2, t0=0.5)))  # response is factor
  expect_identical(sum(log(ddiffusion(rt1$rt, as.numeric(rt1$response), a=1, v=2, t0=0.5))),
                   ref)
  expect_identical(sum(log(ddiffusion(rt1$rt, as.character(rt1$response), a=1, v=2, t0=0.5))),
                   ref)
  expect_identical(sum(log(ddiffusion(rt1, a=1, v=2, t0=0.5))), ref)

  rt2 <- rt1[order(rt1$rt),]

  ref2 <- pdiffusion(rt2$rt, rt2$response, a=1, v=2, t0=0.5)
  expect_identical(pdiffusion(rt2$rt, as.numeric(rt2$response), a=1, v=2, t0=0.5), ref2)
  expect_identical(pdiffusion(rt2$rt, as.character(rt2$response), a=1, v=2, t0=0.5), ref2)
  expect_identical(pdiffusion(rt2, a=1, v=2, t0=0.5), ref2)

#   rt3 <- data.frame(p = rep(seq(0.1, 0.9, 0.2), 2),
#                     response = rep(c("upper", "lower"), each = 5))

  rt3 <- data.frame(p = rep(c(0.05, 0.1), 2),
                    response = rep(c("upper", "lower"), each = 2))
  ref3 <- qdiffusion(rt3$p, rt3$response, a=1, v=2, t0=0.5)
  expect_identical(qdiffusion(rt3$p, as.numeric(rt3$response), a=1, v=2, t0=0.5), ref3)
  expect_identical(qdiffusion(rt3$p, as.character(rt3$response), a=1, v=2, t0=0.5), ref3)
  expect_identical(qdiffusion(rt3, a=1, v=2, t0=0.5), ref3)

})


test_that("qdiffusion is equivalent to manual calculation",{
  p11_fit <- structure(list(par = structure(c(1.32060063610882, 3.27271614698074, 0.338560144920614, 0.34996447540773, 0.201794924457386, 1.05516829794661), .Names = c("a", "v", "t0", "sz", "st0", "sv"))))
  q <- c(0.1, 0.3, 0.5, 0.7, 0.9)

#   i_pdiffusion <- function(x, args, value, response) {
#     abs(value - do.call(pdiffusion, args = c(rt = x, args, response = response)))
#   }
  #pred_dir <- sapply(q*prop_correct, function(y) optimize(i_pdiffusion, c(0, 3), args = as.list(p11_fit$par), value = y, response = "upper")[[1]])

  expect_equal(qdiffusion(q, response = "upper", a=p11_fit$par["a"], v=p11_fit$par["v"], t0=p11_fit$par["t0"], sz=p11_fit$par["sz"]*p11_fit$par["a"], st0=p11_fit$par["st0"], sv=p11_fit$par["sv"], scale_p = TRUE),c(0.474993255765253, 0.548947327845059, 0.607841745594437, 0.681887193854516, 0.844859938530477), tolerance=0.0001)

  expect_equal(suppressWarnings(qdiffusion(q, response = "lower", a=p11_fit$par["a"], v=p11_fit$par["v"], t0=p11_fit$par["t0"], sz=p11_fit$par["sz"]*p11_fit$par["a"], st0=p11_fit$par["st0"], sv=p11_fit$par["sv"])),as.numeric(rep(NA, 5)))
})

test_that("s works as expected", {
  set.seed(1)
  x <- rdiffusion(n = 100, a = 1, v = 2, t0 = 0.3, z = 0.5, s = 1)
  set.seed(1)
  y <- rdiffusion(n = 100, a = 0.1, v = 0.2, t0 = 0.3, z = 0.05, s = 0.1)
  expect_identical(x, y)
  set.seed(1)
  z <- rdiffusion(n = 100, a = 0.1, v = 0.2, t0 = 0.3, s = 0.1)
  expect_identical(x, z)
  expect_identical(
    ddiffusion(x[x$response == "upper", "rt"], a = 1, v = 2, t0 = 0.3, z = 0.5, s=1),
    ddiffusion(x[x$response == "upper", "rt"], a = 0.1, v = 0.2, t0 = 0.3, z = 0.05, s=0.1)
    )
  expect_identical(
    ddiffusion(x[x$response == "upper", "rt"], a = 1, v = 2, t0 = 0.3, z = 0.5, s=1),
    ddiffusion(x[x$response == "upper", "rt"], a = 0.1, v = 0.2, t0 = 0.3, s=0.1)
    )
  expect_identical(
    pdiffusion(sort(x[x$response == "upper", "rt"]), a = 1, v = 2, t0 = 0.3, z = 0.5, s=1),
    pdiffusion(sort(x[x$response == "upper", "rt"]), a = 0.1, v = 0.2, t0 = 0.3, z = 0.05, s=0.1)
  )
  expect_identical(
    pdiffusion(sort(x[x$response == "upper", "rt"]), a = 1, v = 2, t0 = 0.3, z = 0.5, s=1),
    pdiffusion(sort(x[x$response == "upper", "rt"]), a = 0.1, v = 0.2, t0 = 0.3, s=0.1)
  )
  expect_identical(
    qdiffusion(0.6, a = 1, v = 2, t0 = 0.3, z = 0.5, s=1),
    qdiffusion(0.6, a = 0.1, v = 0.2, t0 = 0.3, z = 0.05, s=0.1)
  )
  expect_identical(
    qdiffusion(0.6, a = 1, v = 2, t0 = 0.3, z = 0.5, s=1),
    qdiffusion(0.6, a = 0.1, v = 0.2, t0 = 0.3, s=0.1)
  )
})

test_that("scale_p works as expected", {
  (max_p <- pdiffusion(20, a=1, v=2, t0=0.5, st0=0.2, sz = 0.1, sv = 0.5, response="u"))
  # [1] 0.8705141
  # to get predicted quantiles, scale required quantiles by maximally predicted response rate:
  qs <- c(.1, .3, .5, .7, .9)
  expect_equal(qdiffusion(qs*max_p, a=1, v=2, t0=0.5, st0=0.2, sz = 0.1, sv = 0.5, response="u"),
                    qdiffusion(qs, a=1, v=2, t0=0.5, st0=0.2, sz = 0.1, sv = 0.5, response="u", scale_p = TRUE))

})

test_that("rdiffusion recovers Table 1 from Wagenmakers et al. (2007)", {
  set.seed(4)
  n <- 1e4 # number of samples
  # take parameter valeus from Table 2 and set s to 0.1
  george <- rdiffusion(n, a = 0.12, v = 0.25, t0 = 0.3, s = 0.1)
  rich   <- rdiffusion(n, a = 0.12, v = 0.25, t0 = 0.25, s = 0.1)
  amy    <- rdiffusion(n, a = 0.08, v = 0.25, t0 = 0.3, s = 0.1)
  mark   <- rdiffusion(n, a = 0.08, v = 0.25, t0 = 0.25, s = 0.1)

  george$id <- "george"
  rich$id <- "rich"
  amy$id <- "amy"
  mark$id <- "mark"

  wag <- rbind(george, rich, amy, mark)
  wag$id <- factor(wag$id, levels = c("george", "rich", "amy", "mark"))

  expect_equal(aggregate(rt ~ id, wag, mean)$rt, c(0.517, 0.467, 0.422, 0.372), tolerance = 0.003)

  expect_equal(aggregate(as.numeric(response)-1 ~ id, wag, mean)[,2], c(0.953, 0.953, 0.881, 0.881), tolerance = 0.01)

  expect_equal(aggregate(rt ~ id, wag, var)$rt, c(0.024, 0.024, 0.009, 0.009), tolerance = 0.01)
})


test_that("pdiffusion recovers proportions of Table 1 from Wagenmakers et al. (2007)", {

  props <- pdiffusion(rep(Inf, 4), a = rep(c(0.12, 0.08), each = 2), v = 0.25, t0 = c(0.3, 0.25), s = 0.1)
  expect_equal(props, c(0.953, 0.953, 0.881, 0.881), tolerance = 0.001)

  props <- pdiffusion(rep(Inf, 4), a = rep(c(0.12, 0.08), each = 2), v = 0.25, t0 = c(0.3, 0.25), z = rep(c(0.06, 0.04), each = 2), s = 0.1)
  expect_equal(props, c(0.953, 0.953, 0.881, 0.881), tolerance = 0.001)
})
