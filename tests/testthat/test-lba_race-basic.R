
context("LBA race functions= f*(1-F)")

test_that("dLBA with norm and posdrift=FALSE works as expected", {
  rt <- 2
  vs <- c(0.8, 0.4)
  sd <- 0.2
  
  f <- dlba_norm(rt, A=0.5, b=1, t0 = 0.5, mean_v=vs[1], sd_v=sd, posdrift = FALSE)
  F <- plba_norm(rt, A=0.5, b=1, t0 = 0.5, mean_v=vs[2], sd_v=sd, posdrift = FALSE)
  
  expect_equivalent(n1PDF(rt, A=0.5, b=1, t0 = 0.5, mean_v=vs, sd_v=sd, 
                          args.dist = list(posdrift = FALSE)), f*(1-F))
})

test_that("dLBA with lnorm works as expected", {
  rt <- 2
  vs <- c(0.8, 0.4)
  sd <- 1
  
  f <- dlba_lnorm(rt, A=0.5, b=1, t0 = 0.5, meanlog_v = vs[1], sdlog_v = sd)
  
  F <- plba_lnorm(rt, A=0.5, b=1, t0 = 0.5, meanlog_v = vs[2], sdlog_v = sd)
  
  expect_equivalent(f*(1-F), 
                    n1PDF(rt, A=0.5, b=1, t0 = 0.5,distribution = "lnorm", 
                          meanlog_v = vs, sdlog_v = sd))
})

test_that("dLBA with gamma works as expected", {
  rt <- 2
  vs <- c(1, 1.2)
  sd <- 1
  
  f <- dlba_gamma(rt, A=0.5, b=1, t0 = 0.5, shape_v = vs[1], scale_v = sd)
  
  F <- plba_gamma(rt, A=0.5, b=1, t0 = 0.5, shape_v = vs[2], scale_v = sd)
  
  expect_equivalent(f*(1-F), 
                    n1PDF(rt, A=0.5, b=1, t0 = 0.5,distribution = "gamma", 
                          shape_v = vs, scale_v = sd))
})

test_that("dLBA with frechet works as expected", {
  rt <- 2
  vs <- c(1, 1.2)
  sd <- 1
  
  f <- dlba_frechet(rt, A=0.5, b=1, t0 = 0.5, shape_v = vs[1], scale_v = sd)
  
  F <- plba_frechet(rt, A=0.5, b=1, t0 = 0.5, shape_v = vs[2], scale_v = sd)
  
  expect_equivalent(f*(1-F), 
                    n1PDF(rt, A=0.5, b=1, t0 = 0.5,distribution = "frechet", 
                          shape_v = vs, scale_v = sd))
})