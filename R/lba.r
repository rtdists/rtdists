#' The linear Ballistic Accumulator (LBA)
#' 
#' Density, distribution function, and random generation for the LBA model with different distribution functions underlying the drift rate: Normal (\code{norm}), Gamma (\code{gamma}), Frechet (\code{frechet}), and lognormal (\code{lnorm})
#' 
#' @importFrom evd rfrechet
#' @importFrom msm dtnorm
#' 
#' @name LBA
#' 
NULL

gamma_inc <- function(a,x)  pgamma(x,a,lower=FALSE)*gamma(a)


make_r <- function(drifts, n,b,A,vs,s,t0,st0=0,truncdrifts=TRUE) {
  drifts=drifts[1:n,]
  drifts[drifts<0]=0
  starts=matrix(runif(min=0,max=A,n=n*length(vs)),ncol=length(vs),byrow=TRUE)
  ttf=t((b-t(starts)))/drifts
  rt=apply(ttf,1,min)+t0+runif(min=-st0/2,max=+st0/2,n=n)
  resp=apply(ttf,1,which.min)
  list(rt=rt,resp=resp)
}


#1================2=====================3==============4=================5
# Andrew Terry
# File: generalisedlba-math.r
#1================2=====================3==============4=================5

####### Normal:

#' @rdname LBA
#' @export dlba_norm
dlba_norm <- function(z,x0max,chi,driftrate,sddrift) {
  if (x0max<1e-10) return( (chi/z^2)*dnorm(chi/z,mean=driftrate,sd=sddrift)) 
  zs=z*sddrift ; zu=z*driftrate ; chiminuszu=chi-zu
  chizu=chiminuszu/zs ; chizumax=(chiminuszu-x0max)/zs
  return ((driftrate*(pnorm(chizu)-pnorm(chizumax)) + 
             sddrift*(dnorm(chizumax)-dnorm(chizu)))/x0max)
}

#' @rdname LBA
#' @export plba_norm
plba_norm <- function(z,x0max,chi,driftrate,sddrift) {
  if (x0max<1e-10) return(pnorm(chi/z,mean=driftrate,sd=sddrift,lower.tail=F))
  zs=z*sddrift ; zu=z*driftrate ; chiminuszu=chi-zu ; xx=chiminuszu-x0max;
  chizu=chiminuszu/zs ; chizumax=xx/zs;
  tmp1=zs*(dnorm(chizumax)-dnorm(chizu));
  tmp2=xx*pnorm(chizumax)-chiminuszu*pnorm(chizu);
  return(1+(tmp1+tmp2)/x0max)  
}

#' @rdname LBA
#' @export rlba_norm
rlba_norm <- function(n,b,A,vs,s,t0,st0=0,truncdrifts=TRUE){
  mew=vs ; alpha=s;
  drifts=matrix(rtnorm(n=n*length(vs), mean=vs, sd=s, lower=0),ncol=length(vs),byrow=TRUE)  
  make_r(drifts=drifts, n=n, b=b,A=A, vs=vs, s=s, t0=t0, st0=st0, truncdrifts=truncdrifts)
}


####### Gamma:

#' @rdname LBA
#' @export dlba_gamma
dlba_gamma <- function(z,x0max,chi,driftrate,sddrift) {
  alpha=driftrate ; beta=sddrift;
  min = (chi-x0max)/z ; max = chi/z;
  Gmax=pgamma(max, alpha, rate=beta);Gmin=pgamma(min, alpha, rate=beta)
  Gmax2=pgamma(max, (alpha+1), rate=beta); Gmin2=pgamma(min, (alpha+1), rate=beta)
  zgamma=( ((Gmax2-Gmin2)*gamma(alpha+1))/((Gmax-Gmin)*beta*gamma(alpha)) )
  
  Gmax=pgamma(max, alpha, rate=beta); Gmin=pgamma(min,alpha,rate=beta)
  Gmax2=pgamma(max,(alpha+1),rate=beta); Gmin2=pgamma(min,(alpha+1),rate=beta)
  diffG=function(z,point,alpha, beta) {
    (-point/(z^2))*dgamma(point/z,alpha,rate = beta)
  } #NB:point refers to the constants b OR b-A.
  u = (Gmax2-Gmin2)
  v = (Gmax-Gmin)
  udash = (diffG(z, chi, alpha+1, beta)- diffG(z, (chi-x0max), alpha+1, beta))
  vdash = (diffG(z, chi, alpha, beta)- diffG(z, (chi-x0max), alpha, beta))
  const = gamma(alpha+1)/(beta*gamma(alpha))
  diffzgamma = ((udash*v - vdash*u)/(v^2))*const #quotient rule
  term1 = (Gmax - Gmin)*(zgamma + (z*diffzgamma)); term2=diffG(z,chi,alpha,beta)*((zgamma*z)-chi)
  term3 = diffG(z,(chi-x0max),alpha,beta)*(chi-x0max-(zgamma*z))
  out.value=((term1+term2+term3)/x0max)
  out.value[!is.finite(out.value)] <- 0 # Set NaN or -Inf or Inf to pdf=0
  return(out.value)
}


#' @rdname LBA
#' @export plba_gamma  
plba_gamma <- function(z,x0max,chi,driftrate,sddrift) {
  alpha=driftrate ; beta=sddrift;
  min = (chi-x0max)/z ; max = chi/z;
  Gmax=pgamma(max, alpha, rate=beta);Gmin=pgamma(min, alpha, rate=beta)
  Gmax2=pgamma(max, (alpha+1), rate=beta); Gmin2=pgamma(min, (alpha+1), rate=beta)
  zgamma= ((Gmax2-Gmin2)*gamma(alpha+1))/((Gmax-Gmin)*beta*gamma(alpha)) 
  
  term1 = ((z*zgamma) - chi)/x0max; term2 = (chi-x0max-(z*zgamma))/x0max; 
  pmax = pgamma(max, alpha, rate = beta); pmin = pgamma(min, alpha, rate = beta)
  out.value=(1 + pmax*term1 + pmin*term2)
  out.value[z==Inf] <- 1 # term1=Inf and term2=-Inf cancel in this case
  out.value[!is.finite(out.value)] <- 0 # Set NaN or -Inf to CDF=0
  return(out.value)
}

#' @rdname LBA
#' @export rlba_gamma
rlba_gamma <- function(n,b,A,vs,s,t0,st0=0,truncdrifts=TRUE){
  alpha=vs ; beta=s
  # beta is rate and alpha is shape
  drifts=matrix(rgamma(n=n*length(alpha),shape = alpha,rate = beta),ncol=length(vs),byrow=TRUE)
  
  make_r(drifts=drifts, n=n, b=b,A=A, vs=vs, s=s, t0=t0, st0=st0, truncdrifts=truncdrifts)
}


####### Frechet:

#' @rdname LBA
#' @export dlba_frechet
dlba_frechet <- function(z,x0max,chi,driftrate,sddrift) {
  if (any(c(chi,chi-x0max,driftrate,sddrift)<0)) return(rep(0,length(z))) #protection for pfrechet()
  mew=1/driftrate;  alpha=sddrift;
  z=pmax(z,0)
  min = (chi-x0max)/z; max = chi/z;
  Gmax  = pfrechet(max, loc=0, scale=1/mew, shape=alpha)
  Gmin  = pfrechet(min, loc=0, scale=1/mew, shape=alpha)
  D = Gmax - Gmin
  gam = gamma_inc(1-(1/alpha), (mew*max)^(-alpha))-gamma_inc(1-(1/alpha), (mew*min)^(-alpha))
  zfrechet <- gam/(mew*D)
  diffG1 = ((-chi/(z^2))*dfrechet(chi/z, loc=0, scale=1/mew, shape=alpha))
  diffG2 = ((-(chi-x0max)/(z^2))*dfrechet((chi-x0max)/z, loc=0, scale=1/mew, shape=alpha))    
  diffD = diffG1 - diffG2    
  diffgam = (-alpha*(((mew*chi)^(-alpha+1))/(z^(-alpha+2)))*exp(-(mew*chi/z)^(-alpha))) - 
    (-alpha*(((mew*(chi-x0max))^(-alpha+1))/(z^(-alpha+2)))*exp(-(mew*(chi-x0max)/z)^(-alpha)))
  diffzfrechet = (mew^(-1))*(((-D^(-2))*diffD)*gam + (diffgam*(D^(-1))))
  term1 = (Gmax - Gmin)*(zfrechet + (z*diffzfrechet))
  term2=diffG1*((zfrechet*z)-chi)
  term3 = diffG2*(chi-x0max-(zfrechet*z))
  out.value=((term1+term2+term3)/x0max)
  out.value[!is.finite(out.value)] <- 0 # Set NaN or -Inf or Inf to pdf=0
  return(out.value)    
}

#' @rdname LBA
#' @export plba_frechet
plba_frechet <- function(z,x0max,chi,driftrate,sddrift) {                                                    
  if (any(c(chi,chi-x0max,driftrate,sddrift)<0)) return(rep(0,length(z)))#Protection for the pfrechet()
  z=pmax(z,0)
  mew=1/driftrate; alpha=sddrift;
  min = (chi-x0max)/z; max = chi/z;
  pmax = pfrechet(max, loc=0, scale=1/mew, shape=alpha)
  pmin = pfrechet(min, loc=0, scale=1/mew, shape=alpha)
  zfrechet = (gamma_inc(1-(1/alpha),(mew*max)^(-alpha))-gamma_inc(1-(1/alpha),(mew*min)^(-alpha)))/
    (mew*(pmax-pmin))    
  term1 = ((z*zfrechet) - chi)/x0max ;term2 = (chi-x0max-(z*zfrechet))/x0max ; 
  out.value=(1 + pmax*term1 + pmin*term2)
  out.value[z==Inf] <- 1 # term1=Inf and term2=-Inf cancel in this case
  out.value[!is.finite(out.value)] <- 0 # Set NaN or -Inf to CDF=0
  return(out.value)
}


#' @rdname LBA
#' @export rlba_frechet
rlba_frechet <- function(n,b,A,vs,s,t0,st0=0,truncdrifts=TRUE){
  mew=vs ; alpha=s;
  # mew is rate and alpha is shape
  drifts=matrix(rfrechet(n=n*length(mew), loc=0, scale=mew, shape=alpha),ncol=length(vs),byrow=TRUE)
  
  make_r(drifts=drifts, n=n, b=b,A=A, vs=vs, s=s, t0=t0, st0=st0, truncdrifts=truncdrifts)
}


####### Log-Normal:

#' @rdname LBA
#' @export dlba_lnorm
dlba_lnorm <- function(z,x0max,chi,driftrate,sddrift) {
  mean=driftrate ; sd=sddrift;
  min = (chi-x0max)/z ; max = chi/z;
  
  zlognorm = 
    (exp(mean+(sd^2)/2)*(pnorm((log(max)-mean-(sd^2))/sd)-pnorm((log(min)-mean-(sd^2))/sd))) /
    (pnorm((log(max)-mean)/sd)-pnorm((log(min)-mean)/sd))
  Gmax =plnorm(max,meanlog=mean,sdlog=sd) 
  Gmin =plnorm(min,meanlog=mean,sdlog=sd)
  
  u = (pnorm((log(max)-mean-(sd)^2)/sd)-pnorm((log(min)-mean-(sd)^2)/sd))
  v = (pnorm((log(max)-mean)/sd)-pnorm((log(min)-mean)/sd))
  
  udash = (((-1/(sd*z))*dnorm((log(chi/z)-mean-(sd)^2)/sd)) - 
             ((-1/(sd*z))*dnorm((log((chi-x0max)/z)-mean-(sd)^2)/sd)))
  vdash = (((-1/(sd*z))*dnorm((log(chi/z)-mean)/sd)) - 
             ((-1/(sd*z))*dnorm((log((chi-x0max)/z)-mean)/sd)))
  const = exp(mean+((sd)^2)/2)
  
  diffzlognorm = ((udash*v - vdash*u)/(v^2))*const #quotient rule
  term1 = (Gmax - Gmin)*(zlognorm + (z*diffzlognorm))
  term2=((-chi/(z^2))*dlnorm(chi/z,meanlog=mean,sdlog=sd))*((zlognorm*z)-chi)    
  term3 = (chi-x0max-(zlognorm*z))*
    ((-(chi-x0max)/(z^2))*dlnorm((chi-x0max)/z,meanlog=mean,sdlog=sd))
  out.value=((term1+term2+term3)/x0max)
  out.value[!is.finite(out.value)] <- 0 # Set NaN or -Inf or Inf to pdf=0
  return(out.value)  
}

#' @rdname LBA
#' @export plba_lnorm
plba_lnorm <- function(z,x0max,chi,driftrate,sddrift) {
  mean=driftrate ; sd=sddrift
  min = (chi-x0max)/z ; max = chi/z;
  zlognorm = 
    (exp(mean+(sd^2)/2)*(pnorm((log(max)-mean-(sd^2))/sd)-pnorm((log(min)-mean-(sd^2))/sd))) /
    (pnorm((log(max)-mean)/sd)-pnorm((log(min)-mean)/sd))
  term1 = ((z*zlognorm) - chi)/x0max
  term2 = (chi-x0max-(z*zlognorm))/x0max 
  pmax = plnorm(max, meanlog=mean, sdlog=sd) 
  pmin = plnorm(min, meanlog=mean, sdlog=sd)
  out.value=(1 + pmax*term1 + pmin*term2)
  out.value[z==Inf] <- 1 # term1=Inf and term2=-Inf cancel in this case
  out.value[!is.finite(out.value)] <- 0 # Set NaN or -Inf to CDF=0
  return(out.value)
}


#' @rdname LBA
#' @export rlba_lnorm
rlba_lnorm <- function(n,b,A,vs,s,t0,st0=0,truncdrifts=TRUE){
  mean=vs ; sd=s
  drifts=matrix(rlnorm(n=n*length(mean),meanlog = mean,sdlog=sd),ncol=length(vs),byrow=TRUE)
  make_r(drifts=drifts, n=n, b=b, A=A, vs=vs, s=s, t0=t0, st0=st0, truncdrifts=truncdrifts)
}

