
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



