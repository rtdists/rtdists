#' LBA race functions
#' 
#' as described elsewhere.
#' 
#' 
#' @name LBA-race
#' 
#' @example examples/examples.lba-race.R
#' 
NULL

# t = time, A=x0max, b=chi, v=drift, sv=sdI
n1PDFfixedt0 <- function(t,A,b, t0, ..., plba = plba_norm, plba.args = list()) {
  # Generates defective PDF for responses on node #1.
  dots <- list(...)
  if (is.null(names(dots))) stop("... arguments need to be named.")
  n_v <- max(vapply(dots, length, 0))  # Number of responses
  if (n_v>2) {
    tmp=array(dim=c(length(t),n_v-1))
    for (i in 2:n_v) tmp[,i-1] <- do.call(plba, args = c(t=list(t), A=A[i], b=b[i], t0 = t0[i], sapply(dots, "[[", i = i, simplify = FALSE), plba.args = plba.args))
    G <- apply(1-tmp,1,prod)
  } else {
    G <- 1-do.call(plba, args = c(t=list(t), A=A[2], b=b[2], t0 = t0[2], sapply(dots, "[[", i = 2, simplify = FALSE), plba.args = plba.args))
  }
  G*do.call(plba, args = c(t=list(t), A=A[1], b=b[1], t0 = t0[1], sapply(dots, "[[", i = 1, simplify = FALSE), plba.args = plba.args))
}

#' @rdname LBA-race
#' @export
n1PDF <- function(t, A, b, t0, ..., st0=0, plba = plba_norm, plba.args = list()) {
  dots <- list(...)
  if (is.null(names(dots))) stop("... arguments need to be named.")
  n_v <- max(vapply(dots, length, 0))  # Number of responses
  
  if (length(A)<n_v) A <- rep(A,length.out=n_v)
  if (length(b)<n_v) b <- rep(b,length.out=n_v)
  if (length(t0)<n_v) t0 <- rep(t0,length.out=n_v)
  for (i in length(dots)) {
    if (length(dots[[i]]) < n_v) dots[[i]] <- rep(dots[[i]],length.out=n_v)
  }
  if (length(st0)>1) {
    warning("st0 set to st0[1]. Only one non-decision time variability permitted.")
    st0 <- st0[1] # Only ONE non-decision time.
  }
  if (st0==0) return(do.call(n1PDFfixedt0, args = c(t=list(t), A=list(A), b=list(b), t0 = list(t0), dots, plba = plba, plba.args = plba.args)))
  else {
    tmpf <- function(t, A, b, t0, ..., plba, plba.args = list()) do.call(n1PDFfixedt0, args = c(t=list(t), A=list(A), t0 = list(t0), b=list(b), dots, plba = plba, plba.args = plba.args))/st0
    outs=numeric(length(t))
    for (i in 1:length(outs))
      outs[i] <- do.call(integrate, args=c(f=tmpf, lower=t[i]-st0 , upper=t[i], A=list(A), b=list(b), t0 = list(t0), dots, plba = plba, plba.args = plba.args))$value
    return(outs)
  }
}

# t = time, A=x0max, b=chi, v=drift, sv=sdI
n1CDF <- function(t,x0max,chi,drift,sdI,st0=0,posdrift=TRUE,robust=FALSE) {
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
      tmp=try(integrate(f=n1PDF,lower=bounds[i],upper=bounds[i+1],subdivisions=1000,
                        x0max=x0max,chi=chi,drift=drift,sdI=sdI,st0=st0,posdrift=posdrift,robust=robust)$value,silent=T)
      if (is.numeric(tmp)) {outs[i]=tmp;break}
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

