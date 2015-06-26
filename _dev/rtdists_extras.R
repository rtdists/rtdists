# LNR n-choice
# LNR go-nogo
# LBA go-nogo
# LNR stop signal 

####################### LNR n-choice

rlnr <- function (n, meanlog, sdlog, t0, st0 = 0) 
# Race among length(meanlog) accumulators, t0 can be
# a) a scalar, b) a vector of length number of accumulators or
# c) a matrix with 1 row per accumulator, when start times differ on each trial
# st0, range of non-decison time variability, must be a scalar, as the same
# variability is assumed in a common encoding/production stage

{
    n_acc <- length(meanlog)
    dt <- matrix(rlnorm(n = n*n_acc, meanlog = meanlog, sdlog = sdlog), 
        nrow = n_acc)
    dt <- dt + t0
    winner <- apply(dt,2,which.min)    
    if (st0[1]==0) data.frame(RT=dt[cbind(winner,1:n)],R=winner) else
      data.frame(RT=dt[cbind(winner,1:n)]+runif(n,0,st0[1]),R=factor(winner))
}

n1PDFfixedt0.lnr=function(dt,meanlog,sdlog) 
# Generates defective PDF for responses on node= 1
# dt (decison time) is a matrix with length(meanlog) rows, one row for
# each accumulator to allow for different start times
{
  
  dt[1,] <- dlnorm(dt[1,],meanlog[1],sdlog[1])
  for (i in 2:length(meanlog))
    dt[1,] <- dt[1,]*plnorm(dt[i,],meanlog[i],sdlog[i],lower.tail=FALSE)
  dt[1,]
}

n1PDF.lnr <- function(dt,meanlog,sdlog,t0,st0=0)
# dt (decision time) is a vector, meanlog and sdlog have same length = number of 
# accumulators, t0 is the lower bound of non-decision time, it can be:
# a) a scalar, b) a vector of length number of accumulators or
# c) a matrix with 1 row per accumulator, when start times differ on each trial
# st0, range of non-decison time variability, must be a scalar, as the same
# variability is assumed in a common encoding/production stage
{
  
  # NOTE: No need to flag negative dt in dt-t0 as plnorm/dlnorm return 0
             
  if ( st0 < 1e-3 ) # smaller values can cause numerical integration problems
    return(n1PDFfixedt0.lnr(meanlog=meanlog,sdlog=sdlog,dt=matrix(
      pmax(rep(dt,each=length(meanlog))-t0,0),nrow=length(meanlog)))) else
  {    
    
    integrate.f <- function(dt,meanlog,sdlog,t0,st0) 
      n1PDFfixedt0.lnr(meanlog=meanlog,sdlog=sdlog,dt=matrix(
        pmax(rep(dt,each=length(meanlog))-t0,0),nrow=length(meanlog)))/st0
    
    outs <- numeric(length(dt))
    for (i in 1:length(outs)) {
      tmp <- try(integrate(f=integrate.f,
        lower=dt[i]-st0[1],upper=dt[i],
        meanlog=meanlog,sdlog=sdlog,t0=t0,st0=st0[1])$value,silent=TRUE)
      if (is.numeric(tmp)) 
        outs[i] <- tmp else
        outs[i] <- 0
    }
    return(outs)
  }
}


# # Check
# n=1e6
# meanlog=c(.5,1,1); sdlog=c(1,1,1)
# t0=c(.2,1,1)
# sim <- rlnr(n=n,meanlog,sdlog,t0,st0=0)
# dns <- plot.cell.density(sim,C=1,xlim=c(0,7),save.density=TRUE)
# # n1PDFfixedt0 check
# dt=matrix(rep(dns$correct$x,each=length(meanlog))-t0,nrow=length(meanlog))
# d <- n1PDFfixedt0.lnr(dt,meanlog,sdlog)
# # n1PDF check
# d <- n1PDF.lnr(dns$correct$x,meanlog,sdlog,t0,st0=0)
# # red=black?
# plot(dns$correct$x,dns$correct$y,type="l")
# lines(dns$correct$x,d,col="red")
# # integration check
# st0=1
# sim <- rlnr(n=1e6,meanlog,sdlog,t0,st0)
# dns <- plot.cell.density(sim,C=1,c(0,7),save.density=TRUE)
# d <- n1PDF.lnr(dt=dns$correct$x,meanlog,sdlog,t0,st0)
# plot(dns$correct$x,dns$correct$y,type="l")
# lines(dns$correct$x,d,col="red")

####################### LNR go-nogo

rlnrgng <- function (n, meanlog, sdlog, t0, st0 = 0) 
# Race among length(meanlog) accumulators, first of which is a no-go 
# accumulator. For trials with winning first accumulator RT set to NA.   
{
   out <- rlnr(n,meanlog,sdlog,t0,st0)
   out[out$R==1,"RT"] <- NA
   out$R <- factor(as.numeric(out$R))
   data.frame(out)
}
 

n1PDFfixedt0.lnrgng=function(dt,meanlog,sdlog) 
# Same as n1PDFfixedt0.lnr except dt=NA done by integration
{
  
  stopfn <- function(t,meanlog,sdlog) 
  {
    n1PDFfixedt0.lnr(
      matrix(rep(t,each=length(meanlog)),nrow=length(meanlog)),
      meanlog,sdlog
    )
  }
  
  n.trials <- dim(dt)[2]
  out <- numeric(n.trials)
  is.stop <- is.na(dt[1,])
  dt[1,!is.stop] <- n1PDFfixedt0.lnr(dt[,!is.stop,drop=F],meanlog,sdlog)
  tmp <- try(integrate(f=stopfn,lower=0,upper=Inf,
      meanlog=meanlog,sdlog=sdlog)$value,silent=TRUE)
  if (!is.numeric(tmp)) tmp <- 0
  dt[1,is.na(dt[1,])] <- tmp 
  dt[1,]
}

n1PDF.lnrgng <- function(dt,meanlog,sdlog,t0,st0=0)
# Same as n1PDF.lnr except NAs dealt with by integration
{
    
  if ( st0 < 1e-3 ) # smaller values can cause numerical integration problems
    return(n1PDFfixedt0.lnrgng(meanlog=meanlog,sdlog=sdlog,dt=matrix(
      pmax(rep(dt,each=length(meanlog))-t0,0),nrow=length(meanlog)))
    ) else
  {    
    
    integrate.f <- function(dt,meanlog,sdlog,t0,st0) 
      n1PDFfixedt0.lnrgng(meanlog=meanlog,sdlog=sdlog,dt=matrix(pmax(
        rep(dt,each=length(meanlog))-t0,0),nrow=length(meanlog)))/st0
    
    outs <- numeric(length(dt))
    for (i in 1:length(outs)) {
      tmp <- try(integrate(f=integrate.f,
        lower=dt[i]-st0[1],upper=dt[i],
        meanlog=meanlog,sdlog=sdlog,t0=t0,st0=st0[1])$value,silent=TRUE)
      if (is.numeric(tmp)) 
        outs[i] <- tmp else
        outs[i] <- 0
    }
    return(outs)
  }
}


# # Check 
# n=1e6
# meanlog=c(.5,.75); sdlog=c(1,1); t0=.2
# sim <- rlnrgng(n=n,meanlog,sdlog,t0,st0=0)
# dns <- plot.cell.density(sim,xlim=c(0,7),save.density=TRUE)
# d <- n1PDF.lnrgng(dns$'2'$x,meanlog[c(2,1)],sdlog[c(2,1)],t0,st0=0)
# lines(dns$'2'$x,d,col="red")
# 
# # p(Stop check)
# mean(is.na(sim$RT))
# n1PDFfixedt0.lnrgng(dt=matrix(rep(NA,2),ncol=1),meanlog,sdlog)

####################### LBA go-nogo

make.r <- function (drifts, b, A, n_v, t0, st0 = 0, n) 
{
    drifts[drifts < 0] <- 0
    starts <- matrix(runif(min = 0, max = A, n = n * n_v), nrow = n_v)
    ttf <- t0 + (b - starts)/drifts
    resp <- apply(ttf, 2, which.min)
    rt <- ttf[cbind(resp,1:n)] + runif(min = 0, max = st0[1], n = n)
    bad <- !is.finite(rt)
    if (any(bad)) {
        warning(paste(sum(bad), "infinite RTs removed and less than", 
            n, "rts returned"))
        resp <- resp[!bad]
        rt <- rt[!bad]
    }
    data.frame(rt = rt, response = resp)
}

rlba.norm <- function (n,A,b,t0,mean_v,sd_v,st0=0,posdrift=TRUE) 
{
    if (any(b < A)) 
        stop(error_message_b_smaller_A)
    n_v <- length(mean_v)
    if (posdrift) 
        drifts <- matrix(rtnorm(n = n * n_v, mean = mean_v, sd = sd_v, 
            lower = 0), nrow = n_v)
    else drifts <- matrix(rnorm(n = n * n_v, mean = mean_v, sd = sd_v), 
        nrow = n_v)
    make.r(drifts = drifts, b = b, A = A, n_v = n_v, t0 = t0, st0 = st0, n = n)
}

rlba.normgng <- function (n, A, b, t0, mean_v, sd_v, st0 = 0, posdrift = TRUE) 
# Race among length(mean_v) accumulators, first of which is a no-go 
# accumulator. For trials with winning first accumulator RT set to NA. 
{
   out <- rlba.norm(n,A,b,t0,mean_v,sd_v,st0,posdrift)
   out[out$response==1,"rt"] <- NA
   out$response <- factor(as.numeric(out$response))
   data.frame(out)
}


dlba.norm.core <- function (t,A,b,mean_v,sd_v,posdrift=TRUE,robust=FALSE,nn) 
# like dlba_norm_core but t0 dealt with outside (removed from t), and
# default added to nn
{
  
  pnormP <- function (x, mean = 0, sd = 1, lower.tail = TRUE) 
    ifelse(abs(x) < 7, pnorm(x, mean = mean, sd = sd, lower.tail = lower.tail), 
      ifelse(x < 0, 0, 1))
  
  dnormP <- function (x, mean = 0, sd = 1) 
    ifelse(abs(x) < 7, dnorm(x, mean = mean, sd = sd), 0)

    if (robust) {
        pnorm1 <- pnormP
        dnorm1 <- dnormP
    }
    else {
        pnorm1 <- pnorm
        dnorm1 <- dnorm
    }
    if (posdrift) 
        denom <- pmax(pnorm1(mean_v/sd_v), 1e-10)
    else denom <- rep(1, nn)
    if (any(A < 1e-10)) {
        A_small <- A < 1e-10
        out_A <- pmax(0, ((b[A_small]/t[A_small]^2) * 
                            dnorm1(b[A_small]/t[A_small], 
            mean_v[A_small], sd = sd_v[A_small]))/denom[A_small])
        zs <- t[!A_small] * sd_v[!A_small]
        zu <- t[!A_small] * mean_v[!A_small]
        chiminuszu <- b[!A_small] - zu
        chizu <- chiminuszu/zs
        chizumax <- (chiminuszu - A[!A_small])/zs
        out_o <- pmax(0, (mean_v[!A_small] * (pnorm1(chizu) - 
            pnorm1(chizumax)) + sd_v[!A_small] * (dnorm1(chizumax) - 
            dnorm1(chizu)))/(A[!A_small] * denom[!A_small]))
        out <- numeric(nn)
        out[!A_small] <- out_o
        out[A_small] <- out_A
        return(out)
    }
    else {
        zs <- t * sd_v
        zu <- t * mean_v
        chiminuszu <- b - zu
        chizu <- chiminuszu/zs
        chizumax <- (chiminuszu - A)/zs
        return(pmax(0, (mean_v * (pnorm1(chizu) - pnorm1(chizumax)) + 
            sd_v * (dnorm1(chizumax) - dnorm1(chizu)))/(A * denom)))
    }
}

plba.norm.core <- function (t,A,b,mean_v,sd_v,posdrift=TRUE,robust=FALSE,nn) 
# like plba_norm_core but t0 dealt with outside (removed from t), and
# default added to nn
{
  
    pnormP <- function (x, mean = 0, sd = 1, lower.tail = TRUE) 
      ifelse(abs(x) < 7, pnorm(x, mean = mean, sd = sd, lower.tail=lower.tail), 
        ifelse(x < 0, 0, 1))
  
    dnormP <- function (x, mean = 0, sd = 1) 
      ifelse(abs(x) < 7, dnorm(x, mean = mean, sd = sd), 0)

    if (robust) {
        pnorm1 <- pnormP
        dnorm1 <- dnormP
    }
    else {
        pnorm1 <- pnorm
        dnorm1 <- dnorm
    }
    if (posdrift) 
        denom <- pmax(pnorm1(mean_v/sd_v), 1e-10)
    else denom <- 1
    if (any(A < 1e-10)) {
        A_small <- A < 1e-10
        out_A <- pmin(1, pmax(0, (pnorm1(b[A_small]/t[A_small], 
            mean = mean_v[A_small], sd = sd_v[A_small], 
            lower.tail = FALSE))/denom[A_small]))
        zs <- t[!A_small] * sd_v[!A_small]
        zu <- t[!A_small] * mean_v[!A_small]
        chiminuszu <- b[!A_small] - zu
        xx <- chiminuszu - A[!A_small]
        chizu <- chiminuszu/zs
        chizumax <- xx/zs
        tmp1 <- zs * (dnorm1(chizumax) - dnorm1(chizu))
        tmp2 <- xx * pnorm1(chizumax) - chiminuszu * pnorm1(chizu)
        out_o <- pmin(pmax(0, (1 + (tmp1 + tmp2)/A[!A_small])/denom[!A_small]), 
            1)
        out <- numeric(length(mean_v))
        out[!A_small] <- out_o
        out[A_small] <- out_A
        return(out)
    }
    else {
        zs <- t * sd_v
        zu <- t * mean_v
        chiminuszu <- b - zu
        xx <- chiminuszu - A
        chizu <- chiminuszu/zs
        chizumax <- xx/zs
        tmp1 <- zs * (dnorm1(chizumax) - dnorm1(chizu))
        tmp2 <- xx * pnorm1(chizumax) - chiminuszu * pnorm1(chizu)
        return(pmin(pmax(0, (1 + (tmp1 + tmp2)/A)/denom), 1))
    }
}

dlba.norm <- function (t, A, b, mean_v, sd_v, posdrift = TRUE, robust = FALSE) 
{
    nn <- length(t)
    A <- rep(A, length.out = nn)
    b <- rep(b, length.out = nn)
    mean_v <- rep(mean_v, length.out = nn)
    sd_v <- rep(sd_v, length.out = nn)
    if (any(b < A)) # b cannot be smaller than A! 
        return(rep(0,nn))
    dlba.norm.core(t = t, A = A, b = b, mean_v = mean_v, 
        sd_v = sd_v, posdrift = posdrift, robust = robust, nn = nn)
}

plba.norm <- function (t, A, b, mean_v, sd_v, posdrift = TRUE, robust = FALSE) 
{
    nn <- length(t)
    A <- rep(A, length.out = nn)
    b <- rep(b, length.out = nn)
    mean_v <- rep(mean_v, length.out = nn)
    sd_v <- rep(sd_v, length.out = nn)
    if (any(b < A)) # b cannot be smaller than A! 
        return(rep(0,nn))
    plba.norm.core(t = t, A = A, b = b, mean_v = mean_v, 
        sd_v = sd_v, posdrift = posdrift, robust = robust, nn = nn)
}

n1PDFfixedt0.norm=function(dt,A,b,mean_v,sd_v, 
                           posdrift=TRUE,robust = FALSE) 
# Generates defective PDF for responses on node= 1
# dt (decison time) is a matrix with length(mean_v) rows, one row for
# each accumulator to allow for different start times
{
  
  dt[1,] <- dlba.norm(dt[1,],A=A[1],b=b[1],mean_v=mean_v[1],sd_v=sd_v[1],
                      posdrift=posdrift,robust=robust)
  for (i in 2:length(mean_v))
    dt[1,] <- dt[1,]*(1-plba.norm(dt[i,],
      A=A[i],b=b[i],mean_v=mean_v[i],sd_v=sd_v[i],
      posdrift=posdrift,robust=robust))
  dt[1,]
}

 
n1PDFfixedt0.normgng=function(dt,A,b,mean_v,sd_v,st0=0,posdrift=TRUE) 
# Same as n1PDFfixedt0 except dt=NA done by integration
{
  
  stopfn <- function(t,A,b,mean_v,sd_v,st0=0,posdrift=TRUE) 
  {
    n1PDFfixedt0.norm(
      matrix(rep(t,each=length(mean_v)),nrow=length(mean_v)),
      A=A,b=b,mean_v=mean_v,sd_v=sd_v,posdrift=posdrift
    )
  }
  
  n.trials <- dim(dt)[2]
  out <- numeric(n.trials)
  is.stop <- is.na(dt[1,])
  dt[1,!is.stop] <- n1PDFfixedt0.norm(dt[,!is.stop,drop=F],
    A=A,b=b,mean_v=mean_v,sd_v=sd_v,posdrift=posdrift)
  tmp <- try(integrate(f=stopfn,lower=0,upper=Inf,
      A=A,b=b,mean_v=mean_v,sd_v=sd_v,posdrift=posdrift)$value,silent=TRUE)
  if (!is.numeric(tmp)) tmp <- 0
  dt[1,is.na(dt[1,])] <- tmp 
  dt[1,]
}

n1PDF.normgng <- function(dt,A,b,t0,mean_v,sd_v,st0=0,posdrift=TRUE)
# Same as n1PDF except NAs dealt with by integration
{
    
  if ( st0 < 1e-3 ) # smaller values can cause numerical integration problems
    return(n1PDFfixedt0.normgng(
      A=A,b=b,mean_v=mean_v,sd_v=sd_v,posdrift=posdrift,
      dt=matrix(pmax(rep(dt,each=length(mean_v))-t0,0),nrow=length(mean_v))
    )) else
  {    
    
    integrate.f <- function(A,b,t0,mean_v,sd_v,st0,posdrift) 
      n1PDFfixedt0.normgng(
        A=A,b=b,mean_v=mean_v,sd_v=sd_v,posdrift=posdrift,
        dt=matrix(pmax(rep(dt,each=length(mean_v))-t0,0),
                  nrow=length(mean_v)))/st0
    
    outs <- numeric(length(dt))
    for (i in 1:length(outs)) {
      tmp <- try(integrate(f=integrate.f,
        lower=dt[i]-st0[1],upper=dt[i],A=A,b=b,t0=t0,mean_v=mean_v,sd_v=sd_v,
        st0=st0[1],posdrift=posdrift)$value,silent=TRUE)
      if (is.numeric(tmp)) 
        outs[i] <- tmp else
        outs[i] <- 0
    }
    return(outs)
  }
}


# # Check 
# n=1e6
# A=c(1,1);b=c(2,2);t0=.2;mean_v=c(1,0);sd_v=c(1,1);st0=0;posdrift=TRUE
# sim <- rlba.normgng(n,A,b,t0,mean_v,sd_v,st0,posdrift)
# names(sim) <- c("RT","R")
# dns <- plot.cell.density(sim,xlim=c(0,7),save.density=TRUE)
# d <- n1PDF.normgng(dt=dns$'2'$x,A=A[2:1],b=b[2:1],mean_v=mean_v[2:1],
#   sd_v=sd_v[2:1],t0=t0,st0=st0,posdrift=posdrift)
# # d <- n1PDF.normgng(dt=dns$'2'$x,A=A,b=b,mean_v=mean_v,
# #   sd_v=sd_v,t0=t0,st0=st0,posdrift=posdrift)
# 
# lines(dns$'2'$x,d,col="red")
# 
# # p(Stop check)
# mean(is.na(sim$RT))
# n1PDFfixedt0.normgng(dt=matrix(rep(NA,2),ncol=1),A=A,b=b,
#                      mean_v=mean_v,sd_v=sd_v)



####################### LNR stop signal

rlnrss <- function (n, meanlog, sdlog, t0, st0 = 0, SSD=Inf, tf=0) 
# Race among length(meanlog) accumulators, first of which is a stop accumulator.
# Acts the same as rlnr(meanlog=meanlog[-1]) except NA returned for RT when
# winner = 1. Optional SSD arguement can be used to adjust t0 for first 
# accumulator to t0+SSD. SSD can be a scalar or vector length n. 
# For trials with winning first accumulator RT and R set to NA. 
# Adds SSD column to output. tf = trigger failure probability: sets SSD <> Inf
# to Inf with probability tf.
{
  if ( length(SSD)==1 ) SSD <- rep(SSD,n)
  if ( any(is.na(SSD)) || length(SSD) != n )
    stop("SSD cannot have NAs and must be a scalar or same length as n")
  n_acc <- length(meanlog) 
  if ( length(t0)==1 ) t0 <- rep(t0,n_acc)
  if ( !is.matrix(t0) ) 
    t0 <- matrix(rep(t0,times=n),nrow=n_acc)
  t0[1,] <- t0[1,] + SSD
  out <- rlnr(n,meanlog,sdlog,t0,st0)
  if ( tf>0 ) {
    is.tf <- logical(length(SSD))
    is.tf[is.finite(SSD)][as.logical(rbinom(sum(is.finite(SSD)),1,tf))] <- TRUE  
    if ( any(is.tf) ) { 
      out[is.tf,] <- 
        rlnr(sum(is.tf),meanlog[-1],sdlog[-1],t0[-1,is.tf,drop=FALSE],st0)
      out[is.tf,"R"] <- out[is.tf,"R"]+1
    }
  }
  out[out$R==1,"RT"] <- NA
  cbind.data.frame(out,SSD=SSD)
}
 

n1PDFfixedt0.lnrss=function(dt,meanlog,sdlog,SSD,Si,tf) 
# Same as n1PDFfixedt0.lnr except SSD is subtracted from dt[Si,]
# (i.e., stop accumulator RT) and dt=NA done by integration.
# tf=trigger failure, where L = tf*L(N-1)+(1-tf)[L(N)+p(S)],
# L(N-1) = choice race likelihood (no stop accumulator), 
# L(N) = full N unit race likelihood given they did respond, 
# p(S) probability of stop winning
{
  
  stopfn <- function(t,meanlog,sdlog,SSD,Si) 
  {
    dt <- matrix(rep(t+SSD,each=length(meanlog)),nrow=length(meanlog))
    dt[Si,] <- dt[Si,]-SSD
    i <- c(Si,c(1:length(meanlog))[-Si])
    n1PDFfixedt0.lnr(dt[i,],meanlog[i],sdlog[i])
  }
  
  is.stop <- is.na(dt[1,])  
  dt[Si,] <- dt[Si,] - SSD
  if ( any(!is.stop) ) 
  {
    if ( tf > 0 ) 
    {
      dt[1,!is.stop] <- 
        tf*n1PDFfixedt0.lnr(dt[-Si,!is.stop,drop=FALSE],meanlog[-Si],sdlog[-Si]) +
        (1-tf)*n1PDFfixedt0.lnr(dt[,!is.stop,drop=FALSE],meanlog,sdlog)
    } else 
      dt[1,!is.stop] <- n1PDFfixedt0.lnr(dt[,!is.stop,drop=FALSE],meanlog,sdlog)
  }
  if ( any(is.stop) ) for ( i in unique(SSD[is.stop]) ) {
    tmp <- try(integrate(f=stopfn,lower=0,upper=Inf,
      meanlog=meanlog,sdlog=sdlog,SSD=i,Si=Si)$value,silent=TRUE)
    if (!is.numeric(tmp)) tmp <- 0
    dt[1,is.stop & (SSD==i)] <- tmp*(1-tf) 
  }
  dt[1,]
}

n1PDF.lnrss <- function(dt,meanlog,sdlog,t0,st0,SSD,Si,tf)
# Same as n1PDF.lnr except SSD is either a scalar or vector of length(dt)
# stop accumulator must have name "NR"
{
  
  # NOTE: No need to flag negative dt in dt-t0 as plnorm/dlnorm return 0
  
  if ( length(SSD)==1 ) SSD <- rep(SSD,length(dt))
  if (length(SSD) != length(dt))
    stop("SSD must be a scalar or same length as dt")
  if ( st0 < 1e-3 ) # smaller values can cause numerical integration problems
    return(n1PDFfixedt0.lnrss(meanlog=meanlog,sdlog=sdlog,dt=matrix(
      rep(dt,each=length(meanlog))-t0,nrow=length(meanlog)),
      SSD=SSD,Si=Si,tf=tf[1])

    ) else
  {    
    
    integrate.f <- function(dt,meanlog,sdlog,t0,st0,SSD,Si) 
      n1PDFfixedt0.lnrss(meanlog=meanlog,sdlog=sdlog,dt=matrix(
        rep(dt,each=length(meanlog))-t0,nrow=length(meanlog)),
        SSD=SSD,Si=Si,tf=tf[1])/st0
    
    outs <- numeric(length(dt))
    for (i in 1:length(outs)) {
      tmp <- try(integrate(f=integrate.f,
        lower=dt[i]-st0[1],upper=dt[i],
        meanlog=meanlog,sdlog=sdlog,t0=t0,st0=st0[1],
        SSD=SSD[i],Si=Si,tf=tf[1])$value,silent=TRUE)
      if (is.numeric(tmp)) 
        outs[i] <- tmp else
        outs[i] <- 0
    }
    return(outs)
  }
}


# # No fail check, two different SSDs
# n=1e5
# meanlog=c(.75,.75,1); sdlog=c(.5,1,1)
# t0=.2; SSD = rep(c(1,10)/10,each=n/2)
# st0=0; tf=0; 
# sim <- rlnrss(n=n,meanlog,sdlog,t0,st0=st0,SSD=SSD,tf=tf)
# mean(is.na(sim$RT))
# 
# par(mfrow=c(1,2))
# dns1 <- plot.cell.density(sim[sim$SSD==.1,],xlim=c(0,7),save.density=TRUE)
# dns2 <- plot.cell.density(sim[sim$SSD!=.1,],xlim=c(0,7),save.density=TRUE)
# 
# dat <- sim; dat <- dat[!is.na(dat$RT),]; dat$R <- factor(as.character(dat$R))
# round(tapply(dat$RT,dat[,c("R","SSD")],mean),2)
# 
# x <- dns1$'2'$x
# SSD <- c(rep(.1,length(x)),rep(1,length(x)))
# r1 <- c(2,1,3)
# d.r1 <- n1PDF.lnrss(c(x,x),meanlog[r1],sdlog[r1],t0,st0=st0,SSD=SSD,Si=2,tf=tf)
# r2 <- c(3,1,2)
# d.r2 <- n1PDF.lnrss(c(x,x),meanlog[r2],sdlog[r2],t0,st0=st0,SSD=SSD,Si=2,tf=tf)
# 
# par(mfrow=c(2,2))
# # red=black?
# plot(x,dns1$'2'$y,type="l",main="1, SSD=.1",ylab="Density",xlab="RT",
#      ylim=c(0,max(dns1$'2'$y)))
# lines(x,d.r1[1:length(x)],col="red")
# # red=black?
# plot(x,dns1$'3'$y,type="l",main="2: SSD=.1",ylab="Density",xlab="RT",
#      ylim=c(0,max(dns1$'3'$y)))
# lines(x,d.r2[1:length(x)],col="red")
# # red=black?
# plot(x,dns2$'2'$y,type="l",main="1, SSD=1",ylab="Density",xlab="RT",
#      ylim=c(0,max(dns2$'2'$y)))
# lines(x,d.r1[(length(x)+1):(2*length(x))],col="red")
# # red=black?
# plot(x,dns2$'3'$y,type="l",main="2: SSD=1",ylab="Density",xlab="RT",
#      ylim=c(0,max(dns2$'3'$y)))
# lines(x,d.r2[(length(x)+1):(2*length(x))],col="red")
#  
# # p(Stop check)
# tapply(is.na(sim$RT),sim$SSD,mean)
# n1PDFfixedt0.lnrss(matrix(rep(NA,3),ncol=1),meanlog,sdlog,.1,Si=1,tf=tf)/(1-tf)
# n1PDFfixedt0.lnrss(matrix(rep(NA,3),ncol=1),meanlog,sdlog,1,Si=1,tf=tf)/(1-tf)


# # failure check, one SSD
# stopfn <- function(t,meanlog,sdlog,SSD,Si) 
#   {
#     dt <- matrix(rep(t+SSD,each=length(meanlog)),nrow=length(meanlog))
#     dt[Si,] <- dt[Si,]-SSD
#     i <- c(Si,c(1:length(meanlog))[-Si])
#     n1PDFfixedt0.lnr(dt[i,],meanlog[i],sdlog[i])
#   }
# 
# go3fn <- function(t,meanlog,sdlog,SSD,Si) 
#   {
#     dt <- matrix(rep(t,each=length(meanlog)),nrow=length(meanlog))
#     dt[Si,] <- dt[Si,]-SSD
#     n1PDFfixedt0.lnr(dt,meanlog,sdlog)
#   }
# 
# go2fn <- function(t,meanlog,sdlog,SSD,Si) 
#   {
#     dt <- matrix(rep(t,each=length(meanlog)),nrow=length(meanlog))
#     n1PDFfixedt0.lnr(dt,meanlog,sdlog)
#   }
# 
# 
# 
# meanlog=c(-1,0,0.5); sdlog=c(1,1,1); n=1e5
# t0=.2; ssd <- .5; SSD = rep(ssd,n)
# st0=0; tf=.5; 
# sim <- rlnrss(n=n,meanlog,sdlog,t0,st0=st0,SSD=SSD,tf=tf)
# 
# # p(Stop check)
# pS <- n1PDFfixedt0.lnrss(matrix(rep(NA,3),ncol=1),meanlog,sdlog,ssd,Si=1,tf=tf)
# print(pS)
# mean(is.na(sim$RT))
# pS=integrate(f=stopfn,lower=0,upper=Inf,meanlog=meanlog,sdlog=sdlog,SSD=ssd,Si=1)$value*(1-tf)
# 
# lower=-Inf # adds to 1 with this but shouldnt lower=0????
# mean(sim$R==2)
# p2=integrate(f=go3fn,lower=0,upper=Inf,meanlog=meanlog[c(2,3,1)],sdlog=sdlog[c(2,3,1)],
#           SSD=ssd,Si=3)$value*(1-tf)+
# integrate(f=go2fn,lower=0,upper=Inf,meanlog=meanlog[c(2,3)],sdlog=sdlog[c(2,3)])$value*tf
# p2
# 
# mean(sim$R==3)
# p3=integrate(f=go3fn,lower=0,upper=Inf,meanlog=meanlog[c(3,2,1)],sdlog=sdlog[c(3,2,1)],
#           SSD=ssd,Si=3)$value*(1-tf) +
# integrate(f=go2fn,lower=0,upper=Inf,meanlog=meanlog[c(3,2)],sdlog=sdlog[c(3,2)])$value*tf
# p3
# 
# p2+p3+pS
# 
# 
# 
# 
# dns1 <- plot.cell.density(sim,xlim=c(0,7),save.density=TRUE)
# 
# dat <- sim; dat <- dat[!is.na(dat$RT),]; dat$R <- factor(as.character(dat$R))
# round(tapply(dat$RT,dat[,c("R","SSD")],mean),2)
# 
# x <- dns1$'2'$x
# SSD <- c(rep(ssd,length(x)))
# r1 <- c(2,1,3)
# d.r1 <- n1PDF.lnrss(x,meanlog[r1],sdlog[r1],t0,st0=st0,SSD=SSD,Si=2,tf=tf)
# r2 <- c(3,1,2)
# d.r2 <- n1PDF.lnrss(x,meanlog[r2],sdlog[r2],t0,st0=st0,SSD=SSD,Si=2,tf=tf)
# 
# par(mfrow=c(1,2))
# # red=black?
# plot(x,dns1$'2'$y,type="l",main="1",ylab="Density",xlab="RT",
#      ylim=c(0,max(dns1$'2'$y)))
# lines(x,d.r1,col="red")
# # red=black?
# plot(x,dns1$'3'$y,type="l",main="2",ylab="Density",xlab="RT",
#      ylim=c(0,max(dns1$'3'$y)))
# lines(x,d.r2,col="red")



