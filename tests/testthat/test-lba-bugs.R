
context("n1PDF: Known Bugs")

test_that("n1PDF and n1CDF pass arguments correctly", {
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
  o1 <- n1PDF(bh08$rt,A=sp[1],b=bound[1]+sp[1], t0=nond[1], mean_v=c(drift1[1],drift2[1]), sd_v=sddr[1],dist="norm")
  o2 <- n1PDF(bh08$rt,A=sp[1],b=bound[1]+sp[1], t0=nond[1], mean_v=c(drift1[1],drift2[1]), sd_v=sddr[1],dist="norm", args.dist = list(posdrift = TRUE))
  expect_identical(o1, o2)
  o3 <- n1PDF(bh08$rt,A=sp[1],b=bound[1]+sp[1], t0=nond[1], mean_v=c(drift1[1],drift2[1]), sd_v=sddr[1],dist="norm", args.dist = list(posdrift = FALSE))
  expect_false(all(o1 == o3))
  o4 <- n1PDF(bh08$rt,A=sp[1],b=bound[1]+sp[1], t0=nond[1], mean_v=c(drift1[1],drift2[1]), sd_v=sddr[1],dist="norm", args.dist = list(posdrift = FALSE, robust = TRUE))
  expect_false(all(o1 == o4))
  
  c1 <- n1CDF(bh08$rt,A=sp[1],b=bound[1]+sp[1], t0=nond[1], mean_v=c(drift1[1],drift2[1]), sd_v=sddr[1],dist="norm")
  c2 <- n1CDF(bh08$rt,A=sp[1],b=bound[1]+sp[1], t0=nond[1], mean_v=c(drift1[1],drift2[1]), sd_v=sddr[1],dist="norm", args.dist = list(posdrift = TRUE))
  expect_identical(c1, c2)
  c3 <- n1CDF(bh08$rt,A=sp[1],b=bound[1]+sp[1], t0=nond[1], mean_v=c(drift1[1],drift2[1]), sd_v=sddr[1],dist="norm", args.dist = list(posdrift = FALSE))
  expect_false(all(c1 == c3))
  c4 <- n1CDF(bh08$rt,A=sp[1],b=bound[1]+sp[1], t0=nond[1], mean_v=c(drift1[1],drift2[1]), sd_v=sddr[1],dist="norm", args.dist = list(posdrift = FALSE, robust = TRUE))
  expect_false(all(c1 == c4))
  set.seed(NULL)
})

test_that("named parameter vectors do not cause havoc", {
  xx <- rlba_norm(10, A=0.5, b=1, t0 = 0.5, mean_v=1.2, sd_v=0.2)
  expect_is(n1PDF(xx$rt, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1.0), sd_v=0.2, st0 = c(xx = 0.1), silent =TRUE), "numeric")
  expect_is(n1PDF(xx$rt, A=0.5, b=1, t0 = c(aa=0.5), mean_v=c(1.2, 1.0), sd_v=c(xx=0.2), silent =TRUE), "numeric")
  expect_is(n1PDF(xx$rt, A=c(xx=0.5), b=c(A = 1), t0 = 0.5, mean_v=c(1.2, 1.0), sd_v=c(xx=0.2), silent =TRUE), "numeric")
  expect_is(n1PDF(xx$rt, A=0.5, b=1, t0 = c(aa=0.5), mean_v=c(1.2, 1.0), sd_v=0.2, st0 = 0.1, silent =TRUE), "numeric")
  
  expect_is(n1CDF(xx$rt, A=0.5, b=1, t0 = 0.5, mean_v=c(1.2, 1.0), sd_v=0.2, st0 = c(xx = 0.1), silent =TRUE), "numeric")
  expect_is(n1CDF(xx$rt, A=0.5, b=1, t0 = c(aa=0.5), mean_v=c(1.2, 1.0), sd_v=c(xx=0.2), silent =TRUE), "numeric")
  expect_is(n1CDF(xx$rt, A=c(xx=0.5), b=c(A = 1), t0 = 0.5, mean_v=c(1.2, 1.0), sd_v=c(xx=0.2), silent =TRUE), "numeric")
  expect_is(n1CDF(xx$rt, A=0.5, b=1, t0 = c(aa=0.5), mean_v=c(1.2, 1.0), sd_v=0.2, st0 = 0.1, silent =TRUE), "numeric")
})


test_that("PDFs and CDFs do not return NaN for A = 0", {
  expect_true(all(is.finite(dlba_norm(rt = c(0, 0.0000001, 0.5),  A=0, b=1, t0 = 0, mean_v=1.2, sd_v=0.2))))
  expect_true(all(is.finite(dlba_gamma(rt = c(0, 0.0000001, 0.5),  A=0, b=1, t0 = 0, shape_v=1.2, rate_v=0.2))))
  expect_true(all(is.finite(dlba_frechet(rt = c(0, 0.0000001, 0.5),  A=0, b=1, t0 = 0, shape_v=1.2, scale_v=0.2))))
  expect_true(all(is.finite(dlba_lnorm(rt = c(0, 0.0000001, 0.5),  A=0, b=1, t0 = 0, meanlog_v = 1.2, sdlog_v = 0.2))))
  
  expect_true(all(is.finite(plba_norm(rt = c(0, 0.0000001, 0.5),  A=0, b=1, t0 = 0, mean_v=1.2, sd_v=0.2))))
  expect_true(all(is.finite(plba_gamma(rt = c(0, 0.0000001, 0.5),  A=0, b=1, t0 = 0, shape_v=1.2, rate_v=0.2))))
  expect_true(all(is.finite(plba_frechet(rt = c(0, 0.0000001, 0.5),  A=0, b=1, t0 = 0, shape_v=1.2, scale_v=0.2))))
  expect_true(all(is.finite(plba_lnorm(rt = c(0, 0.0000001, 0.5),  A=0, b=1, t0 = 0, meanlog_v = 1.2, sdlog_v = 0.2))))
  
})


test_that("LBA-norm: PDF and CDF work with various parameter values", {
  
  rts <- c(0, 0.0000001, 0.5, 1.5, 2)
  seq_parameters <- seq(0, 1, length.out = 5)
  
  for (A in seq_parameters) {
    for (b in seq_parameters) {
      for (t0 in seq_parameters) {
        for (d1 in seq_parameters) {
          for (d2 in c(0.0001, seq_parameters[-1])) {
            expect_true(all(is.finite(dlba_norm(rt = rts,  A=A, b=(A+b), t0 = t0, mean_v=d1, sd_v=d2))), info = paste0("A=", A, ", b=", b, ", t0=", t0, ", mean_v=", d1, ", sd_v=", d2))
            expect_true(all(is.finite(dlba_norm(rt = rts,  A=A, b=(A+b), t0 = t0, mean_v=d1, sd_v=d2, posdrift = FALSE))), info = paste0("A=", A, ", b=", b, ", t0=", t0, ", mean_v=", d1, ", sd_v=", d2))
            expect_true(all(is.finite(dlba_norm(rt = rts,  A=A, b=(A+b), t0 = t0, mean_v=d1, sd_v=d2, robust = TRUE))), info = paste0("A=", A, ", b=", b, ", t0=", t0, ", mean_v=", d1, ", sd_v=", d2))
            expect_true(all(is.finite(dlba_norm(rt = rts,  A=A, b=(A+b), t0 = t0, mean_v=d1, sd_v=d2, robust = TRUE, posdrift = FALSE))), info = paste0("A=", A, ", b=", b, ", t0=", t0, ", mean_v=", d1, ", sd_v=", d2))
            
            expect_true(all(is.finite(plba_norm(rt = rts,  A=A, b=(A+b), t0 = t0, mean_v=d1, sd_v=d2))), info = paste0("A=", A, ", b=", b, ", t0=", t0, ", mean_v=", d1, ", sd_v=", d2))
            expect_true(all(is.finite(plba_norm(rt = rts,  A=A, b=(A+b), t0 = t0, mean_v=d1, sd_v=d2, posdrift = FALSE))), info = paste0("A=", A, ", b=", b, ", t0=", t0, ", mean_v=", d1, ", sd_v=", d2))
            expect_true(all(is.finite(plba_norm(rt = rts,  A=A, b=(A+b), t0 = t0, mean_v=d1, sd_v=d2, robust = TRUE))), info = paste0("A=", A, ", b=", b, ", t0=", t0, ", mean_v=", d1, ", sd_v=", d2))
            expect_true(all(is.finite(plba_norm(rt = rts,  A=A, b=(A+b), t0 = t0, mean_v=d1, sd_v=d2, robust = TRUE, posdrift = FALSE))), info = paste0("A=", A, ", b=", b, ", t0=", t0, ", mean_v=", d1, ", sd_v=", d2))
          }
        }
      }
    }
  }
  
})

test_that("LBA-gamma: PDF and CDF work with various parameter values", {
  
  rts <- c(0, 0.0000001, 0.5, 1.5, 2)
  seq_parameters <- seq(0, 1, length.out = 5)
  
  for (A in seq_parameters) {
    for (b in seq_parameters) {
      for (t0 in seq_parameters) {
        for (d1 in seq_parameters) {
          for (d2 in c(0.0001, seq_parameters[-1])) {
            suppressWarnings(expect_true(all(is.finite(dlba_gamma(rt = rts,  A=A, b=(A+b), t0 = t0, shape_v=d1, scale_v=d2))), info = paste0("A=", A, ", b=", b, ", t0=", t0, ", shape_v=", d1, ", scale_v=", d2)))
            
            suppressWarnings(expect_true(all(is.finite(plba_gamma(rt = rts,  A=A, b=(A+b), t0 = t0, shape_v=d1, scale_v=d2))), info = paste0("A=", A, ", b=", b, ", t0=", t0, ", shape_v=", d1, ", scale_v=", d2)))
          }
        }
      }
    }
  }
  
})

test_that("LBA-frechet: PDF and CDF work with various parameter values", {
  
  rts <- c(0, 0.0000001, 0.5, 1.5, 2)
  seq_parameters <- seq(0, 1, length.out = 5)
  
  for (A in seq_parameters) {
    for (b in seq_parameters) {
      for (t0 in seq_parameters) {
        for (d1 in c(0.0001, seq_parameters[-1])) {
          for (d2 in c(0.0001, seq_parameters[-1])) {
            expect_true(all(is.finite(dlba_frechet(rt = rts,  A=A, b=(A+b), t0 = t0, shape_v=d1, scale_v=d2))), info = paste0("A=", A, ", b=", b, ", t0=", t0, ", shape_v=", d1, ", scale_v=", d2))
            
            expect_true(all(is.finite(plba_frechet(rt = rts,  A=A, b=(A+b), t0 = t0, shape_v=d1, scale_v=d2))), info = paste0("A=", A, ", b=", b, ", t0=", t0, ", shape_v=", d1, ", scale_v=", d2))
          }
        }
      }
    }
  }
  
})

test_that("LBA-lnorm: PDF and CDF work with various parameter values", {
  
  rts <- c(0, 0.0000001, 0.5, 1.5, 2)
  seq_parameters <- seq(0, 1, length.out = 5)
  
  for (A in seq_parameters) {
    for (b in seq_parameters) {
      for (t0 in seq_parameters) {
        for (d1 in c(0.0001, seq_parameters[-1])) {
          for (d2 in c(0.0001, seq_parameters[-1])) {
            expect_true(all(is.finite(dlba_lnorm(rt = rts,  A=A, b=(A+b), t0 = t0, meanlog_v=d1, sdlog_v=d2))), info = paste0("A=", A, ", b=", b, ", t0=", t0, ", meanlog_v=", d1, ", sdlog_v=", d2))
            expect_true(all(is.finite(dlba_lnorm(rt = rts,  A=A, b=(A+b), t0 = t0, meanlog_v=d1, sdlog_v=d2)), robust = TRUE), info = paste0("A=", A, ", b=", b, ", t0=", t0, ", meanlog_v=", d1, ", sdlog_v=", d2))
            
            expect_true(all(is.finite(plba_lnorm(rt = rts,  A=A, b=(A+b), t0 = t0, meanlog_v=d1, sdlog_v=d2))), info = paste0("A=", A, ", b=", b, ", t0=", t0, ", meanlog_v=", d1, ", sdlog_v=", d2))
            expect_true(all(is.finite(plba_lnorm(rt = rts,  A=A, b=(A+b), t0 = t0, meanlog_v=d1, sdlog_v=d2, robust = TRUE))), info = paste0("A=", A, ", b=", b, ", t0=", t0, ", meanlog_v=", d1, ", sdlog_v=", d2))
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
  
  
  ll2 <- log(n1PDF(bh08$rt,A=sp[1],b=bound[1]+sp[1], t0=nond[1], mean_v=c(drift1[1],drift2[1]), sd_v=sddr[1],dist="norm", args.dist = list(posdrift = FALSE)))
  
  expect_identical(ll1, ll2)
})



