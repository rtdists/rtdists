#' The linear Ballistic Accumulator (LBA)
#' 
#' Density, distribution function, and random generation for the LBA model with 4 parameters: \code{A} (upper value of starting point), \code{b} (response threshold), \code{v} (driftrate), and \code{sv} (inter-trial-variability of drift). In addition, all functions are available with different distribution functions underlying the drift rate variability \code{sv}: Normal (\code{norm}), Gamma (\code{gamma}), Frechet (\code{frechet}), and lognormal (\code{lnorm})
#' 
#' @param t a vector of RTs.
#' @param n desired number of observations.
#' 
#' @param A start point interval or evidence in accumulator before beginning of decision process. Start point varies from trial to trial in the interval [0, \code{A}] (uniform distribution). Average amount of evidence before evidence accumulation across trials is \code{A}/2.
#' @param b response threshold. (\code{b} - \code{A}/2) is a measure of "response caution".
#' @param v drift rate. Rate at which evidence is accumulated. See details for random number generation.
#' @param t0 non-decision time or response time constant (in seconds). Average duration of all non-decisional processes (encoding and response execution).
#' @param sv variability of drift rate. Distribution depends on function. See details for random number generation.
#' @param st0 variability of non-decision time. Uniformly distributed around \code{t0} +/- \code{st0}/2.
#' 
#' @param truncdrifts logical. Should truncated normal be used for \code{rlba_norm} prohibiting drift rates < 0. Default is \code{TRUE}.
#' 
#' 
#' @details For random number generation \code{v} and 
#' 
#' @importFrom evd rfrechet dfrechet pfrechet
#' @importFrom msm rtnorm
#' 
#' @name LBA
#' 
#' @example examples/examples.lba.R
#' 
NULL

gamma_inc <- function(a,x)  pgamma(x,a,lower.tail=FALSE)*gamma(a)


make_r <- function(drifts, n,b,A,n_v,t0,st0=0) {
  drifts <- drifts[1:n,]
  drifts[drifts<0] <- 0
  starts <- matrix(runif(min=0,max=A,n=n*n_v),ncol=n_v,byrow=TRUE)
  ttf <- t((b-t(starts)))/drifts
  rt <- apply(ttf,1,min)+t0+runif(min=-st0/2,max=+st0/2,n=n)
  resp <- apply(ttf,1,which.min)
  data.frame(rt=rt,response=resp)
}

rem_t0 <- function(t, t0) pmax(t - t0, 0)

#1================2=====================3==============4=================5
# Andrew Terry
# File: generalisedlba-math.r
#1================2=====================3==============4=================5

####### Normal:

#' @rdname LBA
#' @export dlba_norm
dlba_norm <- function(t,A,b,v,t0,sv) {
  t <- rem_t0(t, t0)
  if (A<1e-10) return( (b/t^2)*dnorm(b/t,mean=v,sd=sv)) 
  zs <- t*sv
  zu <- t*v
  chiminuszu <- b-zu
  chizu <- chiminuszu/zs
  chizumax <- (chiminuszu-A)/zs
  return((v*(pnorm(chizu)-pnorm(chizumax)) + sv*(dnorm(chizumax)-dnorm(chizu)))/A)
}

#' @rdname LBA
#' @export plba_norm
plba_norm <- function(t,A,b,v,t0,sv) {
  t <- rem_t0(t, t0)
  if (A<1e-10) return(pnorm(b/t,mean=v,sd=sv,lower.tail=F))
  zs <- t*sv
  zu <- t*v 
  chiminuszu <- b-zu
  xx <- chiminuszu-A
  chizu <- chiminuszu/zs
  chizumax <- xx/zs
  tmp1 <- zs*(dnorm(chizumax)-dnorm(chizu))
  tmp2 <- xx*pnorm(chizumax)-chiminuszu*pnorm(chizu)
  return(1+(tmp1+tmp2)/A)  
}

#' @rdname LBA
#' @export rlba_norm
rlba_norm <- function(n,A,b,v,t0,sv,st0=0,truncdrifts=TRUE){
  if (truncdrifts) drifts <- matrix(rtnorm(n=n*length(v), mean=v, sd=sv, lower=0),ncol=length(v),byrow=TRUE)  
  else drifts <- matrix(rnorm(n=n*length(v), mean=v, sd=sv),ncol=length(v),byrow=TRUE)
  make_r(drifts=drifts, n=n, b=b,A=A, n_v=length(v), t0=t0, st0=st0)
}


####### Gamma:

#' @rdname LBA
#' @export dlba_gamma
dlba_gamma <- function(t,A,b,v,t0,sv) {
  t <- rem_t0(t, t0)
  alpha <- v
  beta <- sv
  min <- (b-A)/t
  max <- b/t
  Gmax <- pgamma(max, alpha, rate=beta)
  Gmin <- pgamma(min, alpha, rate=beta)
  Gmax2 <- pgamma(max, (alpha+1), rate=beta)
  Gmin2 <- pgamma(min, (alpha+1), rate=beta)
  zgamma <- ( ((Gmax2-Gmin2)*gamma(alpha+1))/((Gmax-Gmin)*beta*gamma(alpha)) )
  
  diffG <- function(t,point,alpha, beta) {
    (-point/(t^2))*dgamma(point/t,alpha,rate = beta)
  } #NB:point refers to the constants b OR b-A.
  u <- (Gmax2-Gmin2)
  v <- (Gmax-Gmin)
  udash <- (diffG(t, b, alpha+1, beta)- diffG(t, (b-A), alpha+1, beta))
  vdash <- (diffG(t, b, alpha, beta)- diffG(t, (b-A), alpha, beta))
  const <- gamma(alpha+1)/(beta*gamma(alpha))
  diffzgamma <- ((udash*v - vdash*u)/(v^2))*const #quotient rule
  term1 <- (Gmax - Gmin)*(zgamma + (t*diffzgamma))
  term2 <- diffG(t,b,alpha,beta)*((zgamma*t)-b)
  term3 <- diffG(t,(b-A),alpha,beta)*(b-A-(zgamma*t))
  out.value <- ((term1+term2+term3)/A)
  out.value[!is.finite(out.value)] <- 0 # Set NaN or -Inf or Inf to pdf=0
  return(out.value)
}


#' @rdname LBA
#' @export plba_gamma  
plba_gamma <- function(t,A,b,v,t0,sv) {
  t <- rem_t0(t, t0)
  alpha <- v
  beta <- sv
  min <- (b-A)/t
  max <- b/t
  Gmax <- pgamma(max, alpha, rate=beta)
  Gmin <- pgamma(min, alpha, rate=beta)
  Gmax2 <- pgamma(max, (alpha+1), rate=beta)
  Gmin2 <- pgamma(min, (alpha+1), rate=beta)
  zgamma <- ((Gmax2-Gmin2)*gamma(alpha+1))/((Gmax-Gmin)*beta*gamma(alpha)) 
  
  term1 <- ((t*zgamma) - b)/A
  term2 <- (b-A-(t*zgamma))/A
  pmax <- pgamma(max, alpha, rate = beta)
  pmin <- pgamma(min, alpha, rate = beta)
  out.value <- (1 + pmax*term1 + pmin*term2)
  out.value[t==Inf] <- 1 # term1=Inf and term2=-Inf cancel in this case
  out.value[!is.finite(out.value)] <- 0 # Set NaN or -Inf to CDF=0
  return(out.value)
}

#' @rdname LBA
#' @export rlba_gamma
rlba_gamma <- function(n,A,b,v,t0,sv,st0=0){
  alpha <- v
  beta <- sv
  # beta is rate and alpha is shape
  drifts <- matrix(rgamma(n=n*length(alpha),shape = alpha,rate = beta),ncol=length(v),byrow=TRUE)
  
  make_r(drifts=drifts, n=n, b=b,A=A, n_v = length(v), t0=t0, st0=st0)
}


####### Frechet:

#' @rdname LBA
#' @export dlba_frechet
dlba_frechet <- function(t,A,b,v,t0,sv) {
  t <- rem_t0(t, t0)
  if (any(c(b,b-A,v,sv)<0)) return(rep(0,length(t))) #protection for pfrechet()
  mew <- 1/v
  alpha <- sv
  t <- pmax(t,0)
  min <- (b-A)/t
  max <- b/t
  Gmax <- pfrechet(max, loc=0, scale=1/mew, shape=alpha)
  Gmin <- pfrechet(min, loc=0, scale=1/mew, shape=alpha)
  D <- Gmax - Gmin
  gam <- gamma_inc(1-(1/alpha), (mew*max)^(-alpha))-gamma_inc(1-(1/alpha), (mew*min)^(-alpha))
  zfrechet <- gam/(mew*D)
  diffG1 <- ((-b/(t^2))*dfrechet(b/t, loc=0, scale=1/mew, shape=alpha))
  diffG2 <- ((-(b-A)/(t^2))*dfrechet((b-A)/t, loc=0, scale=1/mew, shape=alpha))    
  diffD <- diffG1 - diffG2    
  diffgam <- (-alpha*(((mew*b)^(-alpha+1))/(t^(-alpha+2)))*exp(-(mew*b/t)^(-alpha))) - (-alpha*(((mew*(b-A))^(-alpha+1))/(t^(-alpha+2)))*exp(-(mew*(b-A)/t)^(-alpha)))
  diffzfrechet <- (mew^(-1))*(((-D^(-2))*diffD)*gam + (diffgam*(D^(-1))))
  term1 <- (Gmax - Gmin)*(zfrechet + (t*diffzfrechet))
  term2 <- diffG1*((zfrechet*t)-b)
  term3 <- diffG2*(b-A-(zfrechet*t))
  out.value <- ((term1+term2+term3)/A)
  out.value[!is.finite(out.value)] <- 0 # Set NaN or -Inf or Inf to pdf=0
  return(out.value)    
}

#' @rdname LBA
#' @export plba_frechet
plba_frechet <- function(t,A,b,v,t0,sv) {
  if (any(c(b,b-A,v,sv)<0)) return(rep(0,length(t)))#Protection for the pfrechet()
  t <- pmax(t,0)
  mew <- 1/v
  alpha <- sv
  min <- (b-A)/t
  max <- b/t
  pmax <- pfrechet(max, loc=0, scale=1/mew, shape=alpha)
  pmin <- pfrechet(min, loc=0, scale=1/mew, shape=alpha)
  zfrechet <- (gamma_inc(1-(1/alpha),(mew*max)^(-alpha))-gamma_inc(1-(1/alpha),(mew*min)^(-alpha)))/(mew*(pmax-pmin))    
  term1 <- ((t*zfrechet) - b)/A
  term2 <- (b-A-(t*zfrechet))/A 
  out.value <- (1 + pmax*term1 + pmin*term2)
  out.value[t==Inf] <- 1 # term1=Inf and term2=-Inf cancel in this case
  out.value[!is.finite(out.value)] <- 0 # Set NaN or -Inf to CDF=0
  return(out.value)
}


#' @rdname LBA
#' @export rlba_frechet
rlba_frechet <- function(n,A,b,v,t0,sv,st0=0){
  mew <- v
  alpha <- sv
  # mew is rate and alpha is shape
  drifts <- matrix(rfrechet(n=n*length(mew), loc=0, scale=mew, shape=alpha),ncol=length(v),byrow=TRUE)
  
  make_r(drifts=drifts, n=n, b=b,A=A, n_v=length(v), t0=t0, st0=st0)
}


####### Log-Normal:

#' @rdname LBA
#' @export dlba_lnorm
dlba_lnorm <- function(t,A,b,v,t0,sv) {
  t <- rem_t0(t, t0)
  mean <- v
  sd <- sv
  min <- (b-A)/t
  max <- b/t
  
  zlognorm <- (exp(mean+(sd^2)/2)*(pnorm((log(max)-mean-(sd^2))/sd)-pnorm((log(min)-mean-(sd^2))/sd))) / (pnorm((log(max)-mean)/sd)-pnorm((log(min)-mean)/sd))
  Gmax <- plnorm(max,meanlog=mean,sdlog=sd) 
  Gmin <- plnorm(min,meanlog=mean,sdlog=sd)
  
  u <- (pnorm((log(max)-mean-(sd)^2)/sd)-pnorm((log(min)-mean-(sd)^2)/sd))
  v <- (pnorm((log(max)-mean)/sd)-pnorm((log(min)-mean)/sd))
  
  udash <- (((-1/(sd*t))*dnorm((log(b/t)-mean-(sd)^2)/sd)) - ((-1/(sd*t))*dnorm((log((b-A)/t)-mean-(sd)^2)/sd)))
  vdash <- (((-1/(sd*t))*dnorm((log(b/t)-mean)/sd)) - ((-1/(sd*t))*dnorm((log((b-A)/t)-mean)/sd)))
  const <- exp(mean+((sd)^2)/2)
  
  diffzlognorm <- ((udash*v - vdash*u)/(v^2))*const #quotient rule
  term1 <- (Gmax - Gmin)*(zlognorm + (t*diffzlognorm))
  term2 <- ((-b/(t^2))*dlnorm(b/t,meanlog=mean,sdlog=sd))*((zlognorm*t)-b)
  term3 <- (b-A-(zlognorm*t))*((-(b-A)/(t^2))*dlnorm((b-A)/t,meanlog=mean,sdlog=sd))
  out.value <- ((term1+term2+term3)/A)
  out.value[!is.finite(out.value)] <- 0 # Set NaN or -Inf or Inf to pdf=0
  return(out.value)  
}

#' @rdname LBA
#' @export plba_lnorm
plba_lnorm <- function(t,A,b,v,t0,sv) {
  t <- rem_t0(t, t0)
  mean <- v
  sd <- sv
  min <- (b-A)/t
  max <- b/t
  zlognorm <- (exp(mean+(sd^2)/2)*(pnorm((log(max)-mean-(sd^2))/sd)-pnorm((log(min)-mean-(sd^2))/sd))) / (pnorm((log(max)-mean)/sd)-pnorm((log(min)-mean)/sd))
  term1 <- ((t*zlognorm) - b)/A
  term2 <- (b-A-(t*zlognorm))/A 
  pmax <- plnorm(max, meanlog=mean, sdlog=sd) 
  pmin <- plnorm(min, meanlog=mean, sdlog=sd)
  out.value <- (1 + pmax*term1 + pmin*term2)
  out.value[t==Inf] <- 1 # term1=Inf and term2=-Inf cancel in this case
  out.value[!is.finite(out.value)] <- 0 # Set NaN or -Inf to CDF=0
  return(out.value)
}


#' @rdname LBA
#' @export rlba_lnorm
rlba_lnorm <- function(n,A,b,v,t0,sv,st0=0){
  mean <- v
  sd <- sv
  drifts=matrix(rlnorm(n=n*length(mean),meanlog = mean,sdlog=sd),ncol=length(v),byrow=TRUE)
  make_r(drifts=drifts, n=n, b=b, A=A, n_v=length(v), t0=t0, st0=st0)
}

