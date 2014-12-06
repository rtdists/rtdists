require(msm)

# z = time, A=x0max, b=chi, v=driftrate, sv=sddrift
fptcdf0=function(z,x0max,chi,driftrate,sddrift) {
  if (x0max<1e-10) # LATER solution
    return( pnorm(chi/z,mean=driftrate,sd=sddrift,lower.tail=FALSE) ) 
  zs=z*sddrift ; zu=z*driftrate ; chiminuszu=chi-zu ; xx=chiminuszu-x0max
  chizu=chiminuszu/zs ; chizumax=xx/zs
  tmp1=zs*(dnorm(chizumax)-dnorm(chizu))
  tmp2=xx*pnorm(chizumax)-chiminuszu*pnorm(chizu)
  1+(tmp1+tmp2)/x0max 
}

# z = time, A=x0max, b=chi, v=driftrate, sv=sddrift
fptpdf0=function(z,x0max,chi,driftrate,sddrift) {
  if (x0max<1e-10) # LATER solution
    return( (chi/z^2)*dnorm(chi/z,mean=driftrate,sd=sddrift) ) 
  zs=z*sddrift ; zu=z*driftrate ; chiminuszu=chi-zu
  chizu=chiminuszu/zs ; chizumax=(chiminuszu-x0max)/zs
  (driftrate*(pnorm(chizu)-pnorm(chizumax)) +
             sddrift*(dnorm(chizumax)-dnorm(chizu)))/x0max
}

# protected normal desity and cdf
pnormP <- function(x,mean=0,sd=1,lower.tail=TRUE){
  ifelse(abs(x)<7,pnorm(x,mean=mean,sd=sd,lower.tail=lower.tail),ifelse(x<0,0,1))}  
dnormP <- function(x,mean=0,sd=1){
  ifelse(abs(x)<7,dnorm(x,mean=mean,sd=sd),0)}

# robust version, 3 times slower!
fptcdfR=function(z,x0max,chi,driftrate,sddrift) {
  if (x0max<1e-10) # LATER solution
    return( pmin(1,pmax(0,pnormP(chi/z,mean=driftrate,sd=sddrift,lower.tail=FALSE))) ) 
  zs=z*sddrift ; zu=z*driftrate ; chiminuszu=chi-zu ; xx=chiminuszu-x0max
  chizu=chiminuszu/zs ; chizumax=xx/zs
  tmp1=zs*(dnormP(chizumax)-dnormP(chizu))
  tmp2=xx*pnormP(chizumax)-chiminuszu*pnormP(chizu)
  pmin(pmax(0,1+(tmp1+tmp2)/x0max),1) 
}

# robust version, 3 times slower!
fptpdfR=function(z,x0max,chi,driftrate,sddrift) {
  if (x0max<1e-10) # LATER solution
    return( pmax(0,(chi/z^2)*dnormP(chi/z,mean=driftrate,sd=sddrift)) ) 
  zs=z*sddrift ; zu=z*driftrate ; chiminuszu=chi-zu
  chizu=chiminuszu/zs ; chizumax=(chiminuszu-x0max)/zs
  pmax(0,(driftrate*(pnormP(chizu)-pnormP(chizumax)) +
             sddrift*(dnormP(chizumax)-dnormP(chizu)))/x0max)
}


fptcdf=function(z,x0max,chi,driftrate,sddrift,posdrift=TRUE,robust=FALSE) {
  if ( robust ) {
    if (!posdrift) fptcdfR(z,x0max,chi,driftrate,sddrift) else
      fptcdfR(z,x0max,chi,driftrate,sddrift)/pmax(pnormP(driftrate/sddrift),1e-10)  
  } else {
    if (!posdrift) fptcdf0(z,x0max,chi,driftrate,sddrift) else
      fptcdf0(z,x0max,chi,driftrate,sddrift)/pmax(pnormP(driftrate/sddrift),1e-10)
  }
}

fptpdf=function(z,x0max,chi,driftrate,sddrift,posdrift=TRUE,robust=FALSE) {
  if ( robust ) {
    if (!posdrift) fptpdfR(z,x0max,chi,driftrate,sddrift) else
      fptpdfR(z,x0max,chi,driftrate,sddrift)/pmax(pnormP(driftrate/sddrift),1e-10)  
  } else {
    if (!posdrift) fptpdf0(z,x0max,chi,driftrate,sddrift) else
      fptpdf0(z,x0max,chi,driftrate,sddrift)/pmax(pnormP(driftrate/sddrift),1e-10)
  }
}

# n = number of samples, v = vs, t0 = minimum non-decision time
# st = width of uniform non-decision time
# If all rates are negative returns Inf and resp=1
rlba=function(n,b,A,vs,s,t0,st0=0,posdrift=TRUE){
  if (posdrift) {
    drifts=matrix(rtnorm(mean=vs,sd=s,n=n*length(vs),lower=0),ncol=length(vs),byrow=TRUE)    
  } else {
    drifts=matrix(rnorm(mean=vs,sd=s,n=n*length(vs)),ncol=length(vs),byrow=TRUE)
    drifts[drifts<0]=0
  }
  starts=matrix(runif(min=0,max=A,n=n*length(vs)),ncol=length(vs),byrow=TRUE)
  ttf=t((b-t(starts)))/drifts
  rt=apply(ttf,1,min)+t0+runif(min=0,max=st0,n=n)
  resp=apply(ttf,1,which.min)
  bad <- !is.finite(rt)
  if (any(bad)) {
    warning(paste(sum(bad),"infinite RTs removed"))
    resp <- resp[!bad]
    rt <- rt[!bad]
  }
  list(rt=rt,resp=resp)
}


# t = time, A=x0max, b=chi, v=drift, sv=sdI
.n1PDFfixedt0=function(t,x0max,chi,drift,sdI,posdrift=TRUE,robust=FALSE) {
  # Generates defective PDF for responses on node #1.
  N=length(drift) # Number of responses.
  if (N>2) {
    tmp=array(dim=c(length(t),N-1))
    for (i in 2:N) tmp[,i-1]=fptcdf(z=t,x0max=x0max[i],chi=chi[i],
      driftrate=drift[i],sddrift=sdI[i],posdrift=posdrift,robust=robust)
    G=apply(1-tmp,1,prod)
  } else {
    G=1-fptcdf(z=t,x0max=x0max[2],chi=chi[2],driftrate=drift[2],
      sddrift=sdI[2],posdrift=posdrift,robust=robust)
  }
  G*fptpdf(z=t,x0max=x0max[1],chi=chi[1],driftrate=drift[1],
    sddrift=sdI[1],posdrift=posdrift,robust=robust)
}

.n1PDF=function(t,x0max,chi,drift,sdI,st0=0,posdrift=TRUE,robust=FALSE) {
  N=length(drift) # Number of responses
  if (length(x0max)<N) x0max=rep(x0max,length.out=N)
  if (length(chi)<N) chi=rep(chi,length.out=N)
  if (length(sdI)<N) sdI=rep(sdI,length.out=N)
  if (length(st0)>1) st0=st0[1] # Only ONE non-decision time.
  if (st0==0) return(.n1PDFfixedt0(t,x0max,chi,drift,sdI,posdrift,robust))
  tmpf=function(t,x0max,chi,drift,sdI,st0,posdrift,robust)
    .n1PDFfixedt0(t,x0max,chi,drift,sdI,posdrift,robust)/st0
  outs=numeric(length(t))
  for (i in 1:length(outs)) outs[i]=integrate(f=tmpf,lower=t[i]-st0,upper=t[i],
    x0max=x0max,chi=chi,drift=drift,sdI=sdI,st0=st0,posdrift=posdrift,robust=robust)$value
  outs
}


.n1CDF=function(t,x0max,chi,drift,sdI,st0=0,posdrift=TRUE,robust=FALSE) {  #, browser=FALSE
  # Generates defective CDF for responses on node #1. 
  N=length(drift) # Number of responses
  if (length(x0max)<N) x0max=rep(x0max,length.out=N)
  if (length(chi)<N) chi=rep(chi,length.out=N)
  if (length(sdI)<N) sdI=rep(sdI,length.out=N)
  if (length(st0)>1) stop("Only one value of st0 allowed.")
  if (st0<1e-6) st0=0 # Integral can fail for small st0.
  outs=numeric(length(t)) ; bounds=c(-st0/2,t)
  for (i in 1:length(t)) {
    tmp="error"
    repeat {
      if (bounds[i]>=bounds[i+1]) {outs[i]=0;break}
      #if(i==1 && browser) browser()
      tmp=try(integrate(f=.n1PDF,lower=bounds[i],upper=bounds[i+1],subdivisions=1000,
        x0max=x0max,chi=chi,drift=drift,sdI=sdI,st0=st0,posdrift=posdrift,robust=robust)$value,silent=T)
      if (is.numeric(tmp)) {outs[i]=tmp;break}
      #browser()
      # Try smart lower bound.
      if (bounds[i]<=0) {
  bounds[i]=max(c((chi-0.98*x0max)/(max(mean(drift),drift[1])+2*sdI)[1],0))
	next
      }
      # Try smart upper bound.
      if (bounds[i+1]==Inf) {
	bounds[i+1]=0.02*max(chi)/(mean(drift)-2*mean(sdI))
	next
      }
      stop("Error in n1CDF that I could not catch.")
    }
  }
  cumsum(outs)
}

n1mean=function(x0max,chi,drift,sdI,posdrift=TRUE,robust=FALSE) {
  # Generates mean RT for responses on node #1. 
   pc=n1CDF(Inf,x0max,chi,drift,sdI,posdrift,robust)
   fn=function(t,x0max,chi,drift,sdI,st0=0,pc,posdrift,robust)
     t*n1PDF(t,x0max,chi,drift,sdI,st0,posdrift,robust)/pc
   tmp=integrate(f=fn,lower=0,upper=100*chi,x0max=x0max,chi=chi,pc=pc,
     drift=drift,sdI=sdI,st0=st0,posdrift=posdrift,robust=robust)$value
   list(mean=tmp,p=pc)
}

