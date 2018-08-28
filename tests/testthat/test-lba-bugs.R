
context("dLBA: Known Bugs")

test_that("dLBA: List and trialwise input for A and b", {
  samples <- 2
  A <- runif(4, 0.3, 0.9)
  b <- A+runif(4, 0, 0.5)
  t0 <- runif(2, 0.1, 0.7)
  v1 <- runif(4, 0.5, 1.5)
  v2 <- runif(4, 0.1, 0.5)
  st0 <- runif(1, 0.1, 0.5)
  r_lba <- rLBA(samples, A=A[1], b=b[1], t0 = t0[1], mean_v=v1[1:2], 
                sd_v=v2[1:2])
  
  p1 <- dLBA(rt = r_lba$rt, response = c(1, 2), 
             A=list(A[1:2],A[3:4]), 
             b=b[1], t0 = t0[1], mean_v=v1[1:2], 
             sd_v=v2[1:2], silent = TRUE)
  p2 <- dLBA(rt = r_lba$rt[1], response = 1, 
             A=list(A[1],A[3]), 
             b=b[1], t0 = t0[1], mean_v=v1[1:2], 
             sd_v=v2[1:2], silent = TRUE)
  p3 <- dLBA(rt = r_lba$rt[2], response = 2, 
             A=list(A[2],A[4]), 
             b=b[1], t0 = t0[1], mean_v=v1[1:2], 
             sd_v=v2[1:2], silent = TRUE)
  expect_identical(p1, c(p2, p3))
  
  p2n1 <- n1PDF(rt = r_lba$rt[1], 
             A=list(A[1],A[3]), 
             b=b[1], t0 = t0[1], mean_v=v1[1:2], 
             sd_v=v2[1:2], silent = TRUE)
  expect_identical(p2, p2n1)
  
  p3n1 <- n1PDF(rt = r_lba$rt[2], 
             A=list(A[4],A[2]), 
             b=b[1], t0 = t0[1], mean_v=v1[2:1], 
             sd_v=v2[2:1], silent = TRUE)
  expect_identical(p3, p3n1)
  
  pb1 <- dLBA(r_lba$rt, c(1, 2), 
              A=A[1:2], 
              b=list(b[1:2],b[3:4]), 
              t0 = t0[1], mean_v=v1[1:2], 
              sd_v=v2[1:2], silent = TRUE)
  pb2 <- dLBA(r_lba$rt[1], 1, A=A[1], 
              b=list(b[1],b[3]), 
              t0 = t0[1], mean_v=v1[1:2], 
              sd_v=v2[1:2], silent = TRUE)
  pb3 <- dLBA(r_lba$rt[2], 2, A=A[2], 
              b = list(b[2],b[4]), 
              t0 = t0[1], 1, mean_v=v1[1:2], 
              sd_v=v2[1:2], silent = TRUE)
  expect_identical(pb1, c(pb2, pb3))
  
  pb2n1 <- n1PDF(rt = r_lba$rt[1], 
                A=A[1], 
                b=list(b[1],b[3]),
                t0 = t0[1], mean_v=v1[1:2], 
                sd_v=v2[1:2], silent = TRUE)
  expect_identical(pb2, pb2n1)
  
  pb3n1 <- n1PDF(rt = r_lba$rt[2], 
                A=A[2], 
                b=list(b[4],b[2]),
                t0 = t0[1], mean_v=v1[2:1], 
                sd_v=v2[2:1], silent = TRUE)
  expect_identical(pb3, pb3n1)
})

context("n1PDF: Known Bugs")

test_that("n1PDF and n1CDF pass arguments correctly", {
  skip_if_not_installed("glba")
  data(bh08, package = "glba")
  bh08 <- bh08[bh08$rt>.180&bh08$rt<2,]
  ny <- dim(bh08)[1]
  
  set.seed(3)
  sddr <- rep(0.2,ny)
  sp <- rep(rnorm(1,.3,.02),ny)
  bound <- rep(rnorm(1,.1,.02),ny)
  nond  <- rep(rnorm(1,.2,.02),ny)
  drift1 <- rep(rnorm(1,.75,.05),ny)
  drift2 <- 1-drift1
  
  parsMat <- matrix(c(sddr,sp,bound,nond,drift1,drift2),ncol=6,nrow=ny)
  o1 <- n1PDF(bh08$rt,A=sp[1],b=bound[1]+sp[1], t0=nond[1], 
              mean_v=c(drift1[1],drift2[1]), sd_v=sddr[1],dist="norm")
  o2 <- n1PDF(bh08$rt,A=sp[1],b=bound[1]+sp[1], t0=nond[1], 
              mean_v=c(drift1[1],drift2[1]), sd_v=sddr[1],dist="norm", 
              args.dist = list(posdrift = TRUE))
  expect_identical(o1, o2)
  o3 <- n1PDF(bh08$rt,A=sp[1],b=bound[1]+sp[1], t0=nond[1], 
              mean_v=c(drift1[1],drift2[1]), sd_v=sddr[1],dist="norm", 
              args.dist = list(posdrift = FALSE))
  expect_false(all(o1 == o3))
  o4 <- n1PDF(bh08$rt,A=sp[1],b=bound[1]+sp[1], t0=nond[1], 
              mean_v=c(drift1[1],drift2[1]), sd_v=sddr[1],dist="norm", 
              args.dist = list(posdrift = FALSE, robust = TRUE))
  expect_false(all(o1 == o4))
  
  c1 <- n1CDF(bh08$rt,A=sp[1],b=bound[1]+sp[1], t0=nond[1], 
              mean_v=c(drift1[1],drift2[1]), sd_v=sddr[1],dist="norm")
  c2 <- n1CDF(bh08$rt,A=sp[1],b=bound[1]+sp[1], t0=nond[1], 
              mean_v=c(drift1[1],drift2[1]), sd_v=sddr[1],dist="norm", 
              args.dist = list(posdrift = TRUE))
  expect_identical(c1, c2)
  c3 <- n1CDF(bh08$rt,A=sp[1],b=bound[1]+sp[1], t0=nond[1], 
              mean_v=c(drift1[1],drift2[1]), sd_v=sddr[1],dist="norm", 
              args.dist = list(posdrift = FALSE))
  expect_false(all(c1 == c3))
  c4 <- n1CDF(bh08$rt,A=sp[1],b=bound[1]+sp[1], t0=nond[1], 
              mean_v=c(drift1[1],drift2[1]), sd_v=sddr[1],dist="norm", 
              args.dist = list(posdrift = FALSE, robust = TRUE))
  expect_false(all(c1 == c4))
  set.seed(NULL)
})

test_that("named parameter vectors do not cause havoc", {
  xx <- rLBA(10, A=0.5, b=1, t0 = 0.5, mean_v=1.2, sd_v=0.2)
  expect_is(n1PDF(xx$rt, A=0.5, b=1, t0 = 0.5, 
                  mean_v=c(1.2, 1.0), sd_v=0.2, 
                  st0 = c(xx = 0.1), silent =TRUE), 
            "numeric")
  expect_is(n1PDF(xx$rt, A=0.5, b=1, t0 = c(aa=0.5), 
                  mean_v=c(1.2, 1.0), sd_v=c(xx=0.2), 
                  silent =TRUE), 
            "numeric")
  expect_is(n1PDF(xx$rt, A=c(xx=0.5), b=c(A = 1), t0 = 0.5, 
                  mean_v=c(1.2, 1.0), sd_v=c(xx=0.2), 
                  silent =TRUE), 
            "numeric")
  expect_is(n1PDF(xx$rt, A=0.5, b=1, t0 = c(aa=0.5), 
                  mean_v=c(1.2, 1.0), sd_v=0.2, st0 = 0.1, 
                  silent =TRUE), 
            "numeric")
  
  expect_is(n1CDF(xx$rt, A=0.5, b=1, t0 = 0.5, 
                  mean_v=c(1.2, 1.0), sd_v=0.2, 
                  st0 = c(xx = 0.1), silent =TRUE), 
            "numeric")
  expect_is(n1CDF(xx$rt, A=0.5, b=1, t0 = c(aa=0.5), 
                  mean_v=c(1.2, 1.0), sd_v=c(xx=0.2), silent =TRUE), 
            "numeric")
  expect_is(n1CDF(xx$rt, A=c(xx=0.5), b=c(A = 1), t0 = 0.5, 
                  mean_v=c(1.2, 1.0), sd_v=c(xx=0.2), silent =TRUE), 
            "numeric")
  expect_is(n1CDF(xx$rt, A=0.5, b=1, t0 = c(aa=0.5), 
                  mean_v=c(1.2, 1.0), sd_v=0.2, st0 = 0.1, silent =TRUE), 
            "numeric")
})


test_that("PDFs and CDFs do not return NaN for A = 0", {
  expect_true(all(is.finite(dlba_norm(rt = c(0, 0.0000001, 0.5),  
                                      A=0, b=1, t0 = 0, 
                                      mean_v=1.2, sd_v=0.2))))
  expect_true(all(is.finite(dlba_gamma(rt = c(0, 0.0000001, 0.5), 
                                       A=0, b=1, t0 = 0, 
                                       shape_v=1.2, rate_v=0.2))))
  expect_true(all(is.finite(dlba_frechet(rt = c(0, 0.0000001, 0.5),  
                                         A=0, b=1, t0 = 0, 
                                         shape_v=1.2, scale_v=0.2))))
  expect_true(all(is.finite(dlba_lnorm(rt = c(0, 0.0000001, 0.5),  
                                       A=0, b=1, t0 = 0, 
                                       meanlog_v = 1.2, sdlog_v = 0.2))))
  
  expect_true(all(is.finite(plba_norm(rt = c(0, 0.0000001, 0.5),  
                                      A=0, b=1, t0 = 0, 
                                      mean_v=1.2, sd_v=0.2))))
  expect_true(all(is.finite(plba_gamma(rt = c(0, 0.0000001, 0.5),  
                                       A=0, b=1, t0 = 0, 
                                       shape_v=1.2, rate_v=0.2))))
  expect_true(all(is.finite(plba_frechet(rt = c(0, 0.0000001, 0.5),  
                                         A=0, b=1, t0 = 0, 
                                         shape_v=1.2, scale_v=0.2))))
  expect_true(all(is.finite(plba_lnorm(rt = c(0, 0.0000001, 0.5),  
                                       A=0, b=1, t0 = 0, 
                                       meanlog_v = 1.2, sdlog_v = 0.2))))
  
})


test_that("LBA-norm: PDF and CDF work with various parameter values", {
  testthat::skip_on_cran()
  rts <- c(0, 0.0000001, 0.5, 1.5, 2)
  seq_parameters <- seq(0, 1, length.out = 5)
  
  for (A in seq_parameters) {
    for (b in seq_parameters) {
      for (t0 in seq_parameters) {
        for (d1 in seq_parameters) {
          for (d2 in c(0.0001, seq_parameters[-1])) {
            expect_true(all(is.finite(
              dlba_norm(
                rt = rts,
                A = A,
                b = (A + b),
                t0 = t0,
                mean_v = d1,
                sd_v = d2
              )
            )),
            info = paste0("A=", A, ", b=", b, ", t0=", t0, 
                          ", mean_v=", d1, ", sd_v=", d2))
            expect_true(all(is.finite(
              dlba_norm(
                rt = rts,
                A = A,
                b = (A + b),
                t0 = t0,
                mean_v = d1,
                sd_v = d2,
                posdrift = FALSE
              )
            )),
            info = paste0("A=", A, ", b=", b, ", t0=", t0, 
                          ", mean_v=", d1, ", sd_v=", d2))
            expect_true(all(is.finite(
              dlba_norm(
                rt = rts,
                A = A,
                b = (A + b),
                t0 = t0,
                mean_v = d1,
                sd_v = d2,
                robust = TRUE
              )
            )),
            info = paste0("A=", A, ", b=", b, ", t0=", t0, 
                          ", mean_v=", d1, ", sd_v=", d2))
            expect_true(all(is.finite(
              dlba_norm(
                rt = rts,
                A = A,
                b = (A + b),
                t0 = t0,
                mean_v = d1,
                sd_v = d2,
                robust = TRUE,
                posdrift = FALSE
              )
            )),
            info = paste0("A=", A, ", b=", b, ", t0=", t0, 
                          ", mean_v=", d1, ", sd_v=", d2))
            
            expect_true(all(is.finite(
              plba_norm(
                rt = rts,
                A = A,
                b = (A + b),
                t0 = t0,
                mean_v = d1,
                sd_v = d2
              )
            )),
            info = paste0("A=", A, ", b=", b, ", t0=", t0, 
                          ", mean_v=", d1, ", sd_v=", d2))
            expect_true(all(is.finite(
              plba_norm(
                rt = rts,
                A = A,
                b = (A + b),
                t0 = t0,
                mean_v = d1,
                sd_v = d2,
                posdrift = FALSE
              )
            )),
            info = paste0("A=", A, ", b=", b, ", t0=", t0, 
                          ", mean_v=", d1, ", sd_v=", d2))
            expect_true(all(is.finite(
              plba_norm(
                rt = rts,
                A = A,
                b = (A + b),
                t0 = t0,
                mean_v = d1,
                sd_v = d2,
                robust = TRUE
              )
            )),
            info = paste0("A=", A, ", b=", b, ", t0=", t0, 
                          ", mean_v=", d1, ", sd_v=", d2))
            expect_true(all(is.finite(
              plba_norm(
                rt = rts,
                A = A,
                b = (A + b),
                t0 = t0,
                mean_v = d1,
                sd_v = d2,
                robust = TRUE,
                posdrift = FALSE
              )
            )),
            info = paste0("A=", A, ", b=", b, ", t0=", t0, 
                          ", mean_v=", d1, ", sd_v=", d2))
          }
        }
      }
    }
  }
  
})

test_that("LBA-gamma: PDF and CDF work with various parameter values", {
  testthat::skip_on_cran()
  rts <- c(0, 0.0000001, 0.5, 1.5, 2)
  seq_parameters <- seq(0, 1, length.out = 5)
  
  for (A in seq_parameters) {
    for (b in seq_parameters) {
      for (t0 in seq_parameters) {
        for (d1 in seq_parameters) {
          for (d2 in c(0.0001, seq_parameters[-1])) {
            suppressWarnings(expect_true(
              all(is.finite(
                dlba_gamma(
                  rt = rts,
                  A = A,
                  b = (A + b),
                  t0 = t0,
                  shape_v = d1,
                  scale_v = d2
                )
              )),
              info = paste0("A=", A, ", b=", b, ", t0=", t0, 
                            ", shape_v=", d1, ", scale_v=", d2)
            ))
            
            suppressWarnings(expect_true(
              all(is.finite(
                plba_gamma(
                  rt = rts,
                  A = A,
                  b = (A + b),
                  t0 = t0,
                  shape_v = d1,
                  scale_v = d2
                )
              )),
              info = paste0("A=", A, ", b=", b, ", t0=", t0, 
                            ", shape_v=", d1, ", scale_v=", d2)
            ))
          }
        }
      }
    }
  }
  
})

test_that("LBA-frechet: PDF and CDF work with various parameter values", {
  testthat::skip_on_cran()
  rts <- c(0, 0.0000001, 0.5, 1.5, 2)
  seq_parameters <- seq(0, 1, length.out = 5)
  
  for (A in seq_parameters) {
    for (b in seq_parameters) {
      for (t0 in seq_parameters) {
        for (d1 in c(0.0001, seq_parameters[-1])) {
          for (d2 in c(0.0001, seq_parameters[-1])) {
            expect_true(all(is.finite(
              dlba_frechet(
                rt = rts,
                A = A,
                b = (A + b),
                t0 = t0,
                shape_v = d1,
                scale_v = d2
              )
            )),
            info = paste0("A=", A, ", b=", b, ", t0=", t0, 
                          ", shape_v=", d1, ", scale_v=", d2))
            
            expect_true(all(is.finite(
              plba_frechet(
                rt = rts,
                A = A,
                b = (A + b),
                t0 = t0,
                shape_v = d1,
                scale_v = d2
              )
            )),
            info = paste0("A=", A, ", b=", b, ", t0=", t0, 
                          ", shape_v=", d1, ", scale_v=", d2))
          }
        }
      }
    }
  }
  
})

test_that("LBA-lnorm: PDF and CDF work with various parameter values", {
  testthat::skip_on_cran()
  rts <- c(0, 0.0000001, 0.5, 1.5, 2)
  seq_parameters <- seq(0, 1, length.out = 5)
  
  for (A in seq_parameters) {
    for (b in seq_parameters) {
      for (t0 in seq_parameters) {
        for (d1 in c(0.0001, seq_parameters[-1])) {
          for (d2 in c(0.0001, seq_parameters[-1])) {
            expect_true(
              all(is.finite(
                dlba_lnorm(
                  rt = rts,
                  A = A,
                  b = (A + b),
                  t0 = t0,
                  meanlog_v = d1,
                  sdlog_v = d2
                )
              )),
              info = paste0("A=", A, ", b=", b, ", t0=", t0, 
                            ", meanlog_v=", d1, ", sdlog_v=", d2)
            )
            expect_true(
              all(is.finite(
                dlba_lnorm(
                  rt = rts,
                  A = A,
                  b = (A + b),
                  t0 = t0,
                  meanlog_v = d1,
                  sdlog_v = d2
                )
              ), robust = TRUE),
              info = paste0("A=", A, ", b=", b, ", t0=", t0, 
                            ", meanlog_v=", d1, ", sdlog_v=", d2)
            )
            
            expect_true(
              all(is.finite(
                plba_lnorm(
                  rt = rts,
                  A = A,
                  b = (A + b),
                  t0 = t0,
                  meanlog_v = d1,
                  sdlog_v = d2
                )
              )),
              info = paste0("A=", A, ", b=", b, ", t0=", t0, 
                            ", meanlog_v=", d1, ", sdlog_v=", d2)
            )
            expect_true(
              all(is.finite(
                plba_lnorm(
                  rt = rts,
                  A = A,
                  b = (A + b),
                  t0 = t0,
                  meanlog_v = d1,
                  sdlog_v = d2,
                  robust = TRUE
                )
              )),
              info = paste0("A=", A, ", b=", b, ", t0=", t0, 
                            ", meanlog_v=", d1, ", sdlog_v=", d2)
            )
          }
        }
      }
    }
  }
  
})

context("glba and rtdists are in agreement")

test_that("glba and rtdists agree", {
  obj <-
    function(rt,pars,loglink,weights) {
      # vectorized loglike function
      # rt: a vector with response times
      # pars: matrix with 4+nrcat parameters on each row to model each rt
      # the drift pars are ordered: the drift for the given response first, the others 
      # after that (order in the remaining drifts does not make a difference)
      for(i in 1:4) if(loglink[i]) pars[,i]=exp(pars[,i])
      ndrift <- dim(pars)[2]-4
      if(ndrift<2) stop("nr of drift pars should at least be two")
      ll <- numeric(length(rt))
      
      ll <- glba:::n1PDF(t=rt-pars[,4], x0max=pars[,2],
                         chi=pars[,2]+pars[,3], sdI=pars[,1], # sdI=0.15, # Scott: I fit chi-x0max.
                         drift=pars[,5:(4+ndrift)])	
      
      # 	return(logl=-sum(log(pmax(weights*ll,1e-10)))) # this has weird effects due to the contaminant model ...
      return(logl=log(weights*ll))
    }
  
  skip_if_not_installed("glba")
  data(bh08, package = "glba")
  # remove extreme RTs
  bh08 <- bh08[bh08$rt>.180&bh08$rt<2,]
  
  ny <- dim(bh08)[1]
  
  set.seed(3)
  sddr <- rep(0.2,ny)
  sp <- rep(rnorm(1,.3,.02),ny)
  bound <- rep(rnorm(1,.1,.02),ny)
  nond  <- rep(rnorm(1,.2,.02),ny)
  drift1 <- rep(rnorm(1,.75,.05),ny)
  drift2 <- 1-drift1
  
  parsMat <- matrix(c(sddr,sp,bound,nond,drift1,drift2),ncol=6,nrow=ny)
  
  head(parsMat)
  
  ll1 <- obj(bh08$rt,parsMat,loglink = c(FALSE,FALSE,FALSE,FALSE),rep(1,ny))
  
  
  ll2 <- log(n1PDF(bh08$rt,A=sp[1],b=bound[1]+sp[1], t0=nond[1], 
                   mean_v=c(drift1[1],drift2[1]), sd_v=sddr[1],
                   dist="norm", args.dist = list(posdrift = FALSE)))
  
  expect_identical(ll1, ll2)
})


test_that("n1PDF works with named lists", {
  rt1 <- rLBA(500, A=0.5, b=1, t0 = 0.5, 
              mean_v=list(a=2.4, b=1.6), sd_v=list(v=1,A=1.2))
  expect_is(sum(log(n1PDF(rt1$rt, A=list(r1=0.5,r2=.5), b=1, t0 = 0.5, 
                          mean_v=list(b=seq(2.0, 2.4, length.out = 500), c=1.6), 
                          sd_v=c(xx=1,hans=1.2)))), 
            "numeric")
  expect_is(sum(log(n1PDF(rt1$rt, A=.5, b=list(r1=0.5,r2=.5), t0 = 0.5, 
                          mean_v=list(b=seq(2.0, 2.4, length.out = 500), c=1.6), 
                          sd_v=c(xx=1,hans=1.2)))), 
            "numeric")
})

test_that("lba_lnorm work with A = 0", {
  A <- 0
  b <- 1  #Can compare to log-normal if b=1
  t0 <- 0
  meanlog_v=0
  sdlog_v=.5
  
  ##########
  set.seed(1)
  check<-rlba_lnorm(1000, A=0, b=1, t0 = 0, meanlog_v=0, sdlog_v=.5)
  rt<-check[,"rt"]
  
  expect_equal(
    sum(log(dlba_lnorm(rt=rt,A=A,b=b,t0=0,meanlog_v = 0,sdlog_v=.5))),
    sum(log(dlnorm(rt,0,.5)))
  )
  
  set.seed(2)
  check<-rlba_lnorm(1000, A=0, b=1, t0 = 0, meanlog_v=.5, sdlog_v=.5)
  rt<-check[,"rt"]
  
  expect_equal(
    sum(log(dlba_lnorm(rt=rt,A=A,b=b,t0=0,meanlog_v = .5,sdlog_v=.5))),
    sum(log(dlnorm(rt,-.5,.5)))  
  )
  
  #x<-  plba_lnorm(rt=rt,A=A,b=b,t0=0,meanlog_v = .5,sdlog_v=.5)- plnorm(rt,-.5,.5)
  expect_equal(
    plba_lnorm(rt=rt,A=A,b=b,t0=0,meanlog_v = .5,sdlog_v=.5),
    plnorm(rt,-.5,.5)
    )
  #CDF is working too
  
  
  #what about b=.5
  set.seed(3)
  check2<-rlba_lnorm(1000, A=0, b=.5, t0 = 0, meanlog_v=0, sdlog_v=.5)
  rt<-check[,"rt"]
  
  expect_equal(
  sum(log(dlba_lnorm(rt=rt,A=A,b=.5,t0=0,meanlog_v = 0,sdlog_v=.5))),
  # [1] -Inf
  sum(log(dlnorm(rt/.5,0,.5)/.5))
  # [1] -261.9191
  )
  
  expect_equal(
    plba_lnorm(rt=rt,A=A,b=.5,t0=0,meanlog_v = 0,sdlog_v=.5),
    plnorm(rt/.5,0,.5)
  )
  
})

test_that("lba_gamma works with A=0", {
  check_gamma <- rlba_gamma(10, A=0.5, b=1, t0 = 0.5, 
                            shape_v=c(1.2, 1), scale_v=c(0.2,0.3))
  rt<-check_gamma[,"rt"]
  expect_equal(
    sum(log(dlba_gamma(rt=rt,A=0.00001, b=1, t0 = 0.5, 
                       shape_v=1.2, scale_v=0.2))),
    sum(log(dlba_gamma(rt=rt,A=0, b=1, t0 = 0.5, 
                       shape_v=1.2, scale_v=0.2)))  
  , tolerance = 0.00001)
  
  expect_equal(
    plba_gamma(rt=rt,A=0.00001, b=1, t0 = 0.5, 
               shape_v=1.2, scale_v=0.2),
    plba_gamma(rt=rt,A=0, b=1, t0 = 0.5, 
               shape_v=1.2, scale_v=0.2)  
  , tolerance = 0.00001)
  
  A <- runif(1, 0.3, 0.9)
  b <- A+runif(1, 0, 0.5)
  t0 <- runif(1, 0.1, 0.7)
  v1 <- runif(2, 0.5, 1.5)
  v2 <- runif(2, 0.1, 0.5)
  
  expect_equal(
    sum(log(dlba_gamma(rt=rt,A=0.00001, b=b, t0 = t0, shape_v=v1, scale_v=v2))),
    sum(log(dlba_gamma(rt=rt,A=0, b=b, t0 = t0, shape_v=v1, scale_v=v2)))  
    , tolerance = 0.00001)
  
  expect_equal(
    sum(log(dlba_gamma(rt=rt,A=0.00001, b=b, t0 = t0, shape_v=v1, scale_v=v2))),
    sum(log(dlba_gamma(rt=rt,A=0, b=b, t0 = t0, shape_v=v1, scale_v=v2)))  
  , tolerance = 0.00001)
  
})

test_that("args.dist is passed through correctly for dLBA, pLBA, qLBA", {
  # see: https://github.com/rtdists/rtdists/issues/7
  d1 <- dLBA(100,1, 10, 100, 0,  
             mean_v=c(3,1), sd_v=c(1,1),
             args.dist = list(posdrift = FALSE))  
  d2 <- dlba_norm(100, 10, 100, 0, 3, 1, posdrift = F, robust = FALSE) * 
    (1-plba_norm(100, 10, 100, 0, 1, 1, posdrift = F, robust = FALSE))
  d3 <- n1PDF(100, 10, 100, 0,  
              mean_v=c(3,1), sd_v=c(1,1), 
              args.dist = list(posdrift = FALSE))
  
  expect_identical(d1, d2)
  expect_identical(d1, d3)
  
  p1 <- pLBA(100,1, 10, 100, 0,  
             mean_v=c(3,1), sd_v=c(1,1),
             args.dist = list(posdrift = FALSE))
  p2 <- n1CDF(100, 10, 100, 0,  
              mean_v=c(3,1), sd_v=c(1,1),
              args.dist = list(posdrift = FALSE))  
  
  expect_identical(p1, p2)
  
  q1 <- qLBA(0.5, 1, 10, 100, 0,  
             mean_v=c(3,1), sd_v=c(1,1), scale_p = TRUE, interval = c(0, 100))
  q2 <- qLBA(0.5, 1, 10, 100, 0,  
             mean_v=c(3,1), sd_v=c(1,1), scale_p = TRUE, interval = c(0, 100), 
             args.dist = list(posdrift = FALSE))
  expect_true(q1 != q2)
  
})
