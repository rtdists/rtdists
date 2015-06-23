#' The linear Ballistic Accumulator (LBA)
#' 
#' Density, distribution function, and random generation for the LBA model with the following parameters: \code{A} (upper value of starting point), \code{b} (response threshold), \code{t0} (non-decision time), and driftrate (\code{v}). All functions are available with different distributions underlying the drift rate: Normal (\code{norm}), Gamma (\code{gamma}), Frechet (\code{frechet}), and log normal (\code{lnorm}).
#' 
#' @param t a vector of RTs.
#' @param n desired number of observations (scalar integer).
#' 
#' @param A start point interval or evidence in accumulator before beginning of decision process. Start point varies from trial to trial in the interval [0, \code{A}] (uniform distribution). Average amount of evidence before evidence accumulation across trials is \code{A}/2.
#' @param b response threshold. (\code{b} - \code{A}/2) is a measure of "response caution". 
#' @param t0 non-decision time or response time constant (in seconds). Lower bound for the duration of all non-decisional processes (encoding and response execution).
#' @param st0 variability of non-decision time, such that \code{t0} is uniformly distributed between \code{t0} and \code{t0} + \code{st0}. Only available in random number generation functions \code{rlba_}.
#' 
#' @param mean_v,sd_v mean and standard deviation of normal distribution for drift rate (\code{norm}). See \code{\link{Normal}}
#' @param shape_v,rate_v,scale_v shape, rate, and scale of gamma (\code{gamma}) and scale and shape of Frechet (\code{frechet}) distributions for drift rate. See \code{\link{GammaDist}} or \code{\link[evd]{frechet}}. For Gamma, scale = 1/shape and shape = 1/scale.
#' @param meanlog_v,sdlog_v mean and standard deviation of lognormal distribution on the log scale for drift rate (\code{lnorm}). See \code{\link{Lognormal}}.
#' 
#' @param posdrift logical. Should driftrates be forced to be positive? Default is \code{TRUE}. (Uses truncated normal for random generation).
#' @param robust logical. Should robust normal distributions be used for \code{norm} and \code{lnorm}? Can be helpful in rare cases but is approximately three times slower than the non-robust versions. Default is \code{FALSE}.
#' 
#' 
#' @details For random number generation at least one of the distribution parameters (i.e., \code{mean_v}, \code{sd_v}, \code{shape_v}, \code{scale_v}, \code{rate_v}, \code{meanlog_v}, and \code{sdlog_v}) should be of length > 1 to receive RTs from multiple responses. Shorter vectors are recycled as necessary.\cr
#' Note that for random number generation from a normal distribution for the driftrate the number of returned samples may be less than the number of requested samples if \code{posdrifts==FALSE}.
#' 
#' @return All functions starting with a \code{d} return the density (PDF), all functions starting with \code{p} return the distribution function (CDF), and all functions starting with \code{r} return random response times and responses (in a \code{data.frame}).
#' 
#' @note Density (i.e., \code{dlba_}) and distribution (i.e., \code{plba_}) functions are vectorized for all parameters (i.e., in case parameters are not of the same length as \code{t}, parameters are recycled). Somewhat inconsistently, the random number generation functions \code{rlba_} accept only scalar inputs.
#' 
#' @references 
#' 
#' Brown, S. D., & Heathcote, A. (2008). The simplest complete model of choice response time: Linear ballistic accumulation. \emph{Cognitive Psychology}, 57(3), 153-178. doi:10.1016/j.cogpsych.2007.12.002
#' 
#' Donkin, C., Averell, L., Brown, S., & Heathcote, A. (2009). Getting more from accuracy and response time data: Methods for fitting the linear ballistic accumulator. \emph{Behavior Research Methods}, 41(4), 1095-1110. doi:10.3758/BRM.41.4.1095
#' 
#' Heathcote, A., & Love, J. (2012). Linear deterministic accumulator models of simple choice. \emph{Frontiers in Psychology}, 3, 292. doi:10.3389/fpsyg.2012.00292
#' 
#' @importFrom evd rfrechet dfrechet pfrechet
#' @importFrom msm rtnorm
#' @importFrom gsl gamma_inc
#' 
#' @name LBA
#' 
#' @example examples/examples.lba.R
#' 
NULL

# protected normal desity and cdf
pnormP <- function(x,mean=0,sd=1,lower.tail=TRUE) ifelse(abs(x)<7,pnorm(x, mean=mean, sd=sd,lower.tail=lower.tail),ifelse(x<0,0,1))
dnormP <- function(x,mean=0,sd=1) ifelse(abs(x)<7,dnorm(x,mean=mean,sd=sd),0)

make_r <- function(drifts, n,b,A,n_v,t0,st0=0) {
  drifts <- drifts[1:n,]
  drifts[drifts<0] <- 0
  starts <- matrix(runif(min=0,max=A,n=n*n_v),ncol=n_v,byrow=TRUE)
  ttf <- t((b-t(starts)))/drifts
  rt <- apply(ttf,1,min)+t0+runif(min=0,max=st0,n=n)
  resp <- apply(ttf,1,which.min)
  bad <- !is.finite(rt)
  if (any(bad)) {
    warning(paste(sum(bad),"infinite RTs removed and less than", n, "rts returned"))
    resp <- resp[!bad]
    rt <- rt[!bad]
  }
  data.frame(rt=rt,response=resp)
}

rem_t0 <- function(t, t0) pmax(t - t0, 0)

check_single_arg <- function(...) {
  mc <- match.call()
  vars <- all.vars(mc)
  arguments <- list(...)
  for(i in seq_along(arguments)) {
    if (length(arguments[[i]]) != 1) stop(paste(vars[i], "needs to be of length 1!"))
    if (!is.numeric(arguments[[i]]) | !is.finite(arguments[[i]])) stop(paste(vars[i], "needs to be numeric and finite!"))
  }
}

check_vector <- function(...) {
  mc <- match.call()
  vars <- all.vars(mc)
  dots <- list(...)
  for (i in seq_along(dots)) {
    if ((vars[i] == "t") && (any(dots[[i]] < 0))) stop("t needs to contain only positive values.") 
    if (!is.vector(dots[[i]], "numeric")) stop(paste(vars[[i]], "needs to be a numeric vector!"))
    if (length(dots[[i]]) < 1) stop(paste(vars[[i]], "needs to have a length >= 1."))
  }
}

error_message_b_smaller_A <- "b cannot be smaller than A!"

####### Normal:

#' @rdname LBA
#' @export dlba_norm
dlba_norm <- function(t,A,b, t0, mean_v, sd_v, posdrift=TRUE, robust = FALSE) {
  #check_single_arg(A=A, b=b, t0=t0, mean_v=mean_v, sd_v=sd_v)
  check_vector(t, A, b, t0, mean_v, sd_v)
  # bring all arguments to length of t
  nn <- length(t)
  A <- rep(A, length.out = nn)
  b <- rep(b, length.out = nn)
  t0 <- rep(t0, length.out = nn)
  mean_v <- rep(mean_v, length.out = nn)
  sd_v <- rep(sd_v, length.out = nn)
  if (any(b < A)) stop(error_message_b_smaller_A)
  dlba_norm_core(t = t, A = A, b = b, t0 = t0, mean_v = mean_v, sd_v = sd_v, posdrift = posdrift, robust = robust, nn = nn)
}

## this functions expects all arguments to have the samel length (which is nn)
dlba_norm_core <- function(t,A,b, t0, mean_v, sd_v, posdrift=TRUE, robust = FALSE, nn) {
  if (robust) { # robust == TRUE uses robust versions of the normal distributions
    pnorm1 <- pnormP
    dnorm1 <- dnormP
  } else {
    pnorm1 <- pnorm
    dnorm1 <- dnorm
  }
  t <- rem_t0(t, t0) # rmove t0 from t
  if (posdrift) denom <- pmax(pnorm1(mean_v/sd_v),1e-10) else denom <- rep(1, nn)
  
  if (any(A<1e-10)) {
    # for A<1e-10 save results in out_A
    A_small <- A<1e-10
    out_A <- pmax(0, ((b[A_small]/t[A_small]^2)*dnorm1(b[A_small]/t[A_small],mean_v[A_small],sd=sd_v[A_small]))/denom[A_small]) 
    # calculate other results into out_o
    zs <- t[!A_small]*sd_v[!A_small]
    zu <- t[!A_small]*mean_v[!A_small]
    chiminuszu <- b[!A_small]-zu
    chizu <- chiminuszu/zs
    chizumax <- (chiminuszu-A[!A_small])/zs
    out_o <- pmax(0,(mean_v[!A_small]*(pnorm1(chizu)-pnorm1(chizumax)) + sd_v[!A_small]*(dnorm1(chizumax)-dnorm1(chizu)))/(A[!A_small]*denom[!A_small]))
    # combine out_A and out_o
    out <- numeric(nn)
    out[!A_small] <- out_o
    out[A_small] <- out_A
    return(out) 
  } else {
    zs <- t*sd_v
    zu <- t*mean_v
    chiminuszu <- b-zu
    chizu <- chiminuszu/zs
    chizumax <- (chiminuszu-A)/zs
    return(pmax(0,(mean_v*(pnorm1(chizu)-pnorm1(chizumax)) + sd_v*(dnorm1(chizumax)-dnorm1(chizu)))/(A*denom))) 
  }
}

#' @rdname LBA
#' @export plba_norm
plba_norm <- function(t,A,b,t0,mean_v, sd_v,posdrift=TRUE, robust = FALSE) {
  check_vector(t, A, b, t0, mean_v, sd_v)
  nn <- length(t)
  A <- rep(A, length.out = nn)
  b <- rep(b, length.out = nn)
  t0 <- rep(t0, length.out = nn)
  mean_v <- rep(mean_v, length.out = nn)
  sd_v <- rep(sd_v, length.out = nn)
  if (any(b < A)) stop(error_message_b_smaller_A)
  plba_norm_core(t = t, A = A, b = b, t0 = t0, mean_v = mean_v, sd_v = sd_v, posdrift = posdrift, robust = robust, nn = nn)
}

plba_norm_core <- function(t,A,b,t0,mean_v, sd_v,posdrift=TRUE, robust = FALSE, nn) {
  if (robust) { # robust == TRUE uses robust versions of the normal distributions
    pnorm1 <- pnormP
    dnorm1 <- dnormP
  } else {
    pnorm1 <- pnorm
    dnorm1 <- dnorm  
  }

  t <- rem_t0(t, t0)
  if (posdrift) denom <- pmax(pnorm1(mean_v/sd_v),1e-10) else denom <- 1
  
  if (any(A<1e-10)) {
    # for A<1e-10 save results in out_A
    A_small <- A<1e-10
    out_A <- pmin(1, pmax(0, (pnorm1(b[A_small]/t[A_small],mean=mean_v[A_small],sd=sd_v[A_small],lower.tail=FALSE))/denom[A_small]))

    # calculate other results into out_o
    zs <- t[!A_small]*sd_v[!A_small]
    zu <- t[!A_small]*mean_v[!A_small]
    chiminuszu <- b[!A_small]-zu
    xx <- chiminuszu-A[!A_small]
    chizu <- chiminuszu/zs
    chizumax <- xx/zs
    tmp1 <- zs*(dnorm1(chizumax)-dnorm1(chizu))
    tmp2 <- xx*pnorm1(chizumax)-chiminuszu*pnorm1(chizu)
    out_o <- pmin(pmax(0,(1+(tmp1+tmp2)/A[!A_small])/denom[!A_small]), 1)
    
    # combine out_A and out_o
    out <- numeric(nn)
    out[!A_small] <- out_o
    out[A_small] <- out_A
    return(out)
  } else {
    zs <- t*sd_v
    zu <- t*mean_v
    chiminuszu <- b-zu
    xx <- chiminuszu-A
    chizu <- chiminuszu/zs
    chizumax <- xx/zs
    tmp1 <- zs*(dnorm1(chizumax)-dnorm1(chizu))
    tmp2 <- xx*pnorm1(chizumax)-chiminuszu*pnorm1(chizu)
    return(pmin(pmax(0,(1+(tmp1+tmp2)/A)/denom), 1))
  }
}

#' @rdname LBA
#' @export rlba_norm
rlba_norm <- function(n,A,b,t0,mean_v, sd_v, st0=0,posdrift=TRUE) {
  check_single_arg(n, A, b, t0, st0)
  if (b < A) stop(error_message_b_smaller_A)
  n_v <- max(length(mean_v), length(sd_v))
  if (posdrift) drifts <- matrix(rtnorm(n=n*n_v, mean=mean_v, sd=sd_v, lower=0),ncol=n_v,byrow=TRUE)  
  else drifts <- matrix(rnorm(n=n*n_v, mean=mean_v, sd=sd_v),ncol=n_v,byrow=TRUE)
  make_r(drifts=drifts, n=n, b=b,A=A, n_v=n_v, t0=t0, st0=st0)
}


####### Gamma:

#' @rdname LBA
#' @export dlba_gamma
dlba_gamma <- function(t,A,b,t0,shape_v,rate_v, scale_v) {
  
  if (!missing(rate_v) && !missing(scale_v)) stop("specify 'rate_v' or 'scale_v', but not both")
  if (missing(rate_v)) rate_v <- 1/scale_v
  
  check_vector(t, A, b=b, t0, shape_v, rate_v)
  
  nn <- length(t)
  A <- rep(A, length.out = nn)
  b <- rep(b, length.out = nn)
  t0 <- rep(t0, length.out = nn)
  shape_v <- rep(shape_v, length.out = nn)
  rate_v <- rep(rate_v, length.out = nn)
  if (any(b < A)) stop(error_message_b_smaller_A)
  
  dlba_gamma_core(t=t,A=A,b=b,t0=t0, shape_v=shape_v, rate_v=rate_v, nn=nn)
  
}


dlba_gamma_core <- function(t,A,b,t0,shape_v, rate_v, nn) {
  t <- rem_t0(t, t0)
  min <- (b-A)/t
  max <- b/t
  
  Gmax <- pgamma(max, shape_v, rate=rate_v)
  Gmin <- pgamma(min, shape_v, rate=rate_v)
  Gmax2 <- pgamma(max, (shape_v+1), rate=rate_v)
  Gmin2 <- pgamma(min, (shape_v+1), rate=rate_v)
  zgamma <- ( ((Gmax2-Gmin2)*gamma(shape_v+1))/((Gmax-Gmin)*rate_v*gamma(shape_v)) )
  
  diffG <- function(t,point,shape_v, rate_v) {
    (-point/(t^2))*dgamma(point/t,shape_v,rate = rate_v)
  } #NB:point refers to the constants b OR b-A.
  u <- (Gmax2-Gmin2)
  v <- (Gmax-Gmin)
  udash <- (diffG(t, b, shape_v+1, rate_v)- diffG(t, (b-A), shape_v+1, rate_v))
  vdash <- (diffG(t, b, shape_v, rate_v)- diffG(t, (b-A), shape_v, rate_v))
  const <- gamma(shape_v+1)/(rate_v*gamma(shape_v))
  diffzgamma <- ((udash*v - vdash*u)/(v^2))*const #quotient rule
  term1 <- (Gmax - Gmin)*(zgamma + (t*diffzgamma))
  term2 <- diffG(t,b,shape_v,rate_v)*((zgamma*t)-b)
  term3 <- diffG(t,(b-A),shape_v,rate_v)*(b-A-(zgamma*t))
  out.value <- ((term1+term2+term3)/A)
  out.value[!is.finite(out.value)] <- 0 # Set NaN or -Inf or Inf to pdf=0
  return(pmax(0, out.value))
}


#' @rdname LBA
#' @export plba_gamma  
plba_gamma <- function(t,A,b,t0,shape_v, rate_v, scale_v) {
  if (!missing(rate_v) && !missing(scale_v)) stop("specify 'rate_v' or 'scale_v', but not both")
  if (missing(rate_v)) rate_v <- 1/scale_v
  
  check_vector(t, A, b=b, t0, shape_v, rate_v)
  
  nn <- length(t)
  A <- rep(A, length.out = nn)
  b <- rep(b, length.out = nn)
  t0 <- rep(t0, length.out = nn)
  shape_v <- rep(shape_v, length.out = nn)
  rate_v <- rep(rate_v, length.out = nn)
  if (any(b < A)) stop(error_message_b_smaller_A)
  
  plba_gamma_core(t=t,A=A,b=b,t0=t0,shape_v=shape_v, rate_v=rate_v, nn=nn)
}

plba_gamma_core <- function(t,A,b,t0,shape_v, rate_v, nn) {
  
  t <- rem_t0(t, t0)
  min <- (b-A)/t
  max <- b/t
  Gmax <- pgamma(max, shape_v, rate=rate_v)
  Gmin <- pgamma(min, shape_v, rate=rate_v)
  Gmax2 <- pgamma(max, (shape_v+1), rate=rate_v)
  Gmin2 <- pgamma(min, (shape_v+1), rate=rate_v)
  zgamma <- ((Gmax2-Gmin2)*gamma(shape_v+1))/((Gmax-Gmin)*rate_v*gamma(shape_v)) 
  
  term1 <- ((t*zgamma) - b)/A
  term2 <- (b-A-(t*zgamma))/A
  pmax <- pgamma(max, shape_v, rate = rate_v)
  pmin <- pgamma(min, shape_v, rate = rate_v)
  out.value <- (1 + pmax*term1 + pmin*term2)
  out.value[t==Inf] <- 1 # term1=Inf and term2=-Inf cancel in this case
  out.value[!is.finite(out.value)] <- 0 # Set NaN or -Inf to CDF=0
  return(pmin(pmax(0, out.value), 1))
}

#' @rdname LBA
#' @export rlba_gamma
rlba_gamma <- function(n,A,b,t0,shape_v, rate_v, scale_v, st0=0) {
  check_single_arg(n, A, b, t0, st0)
  if (b < A) stop(error_message_b_smaller_A)
  if (!missing(rate_v) && !missing(scale_v)) stop("specify 'rate_v' or 'scale_v', but not both")
  if (missing(rate_v)) rate_v <- 1/scale_v
  n_v <- max(length(shape_v), length(rate_v))  
  drifts <- matrix(rgamma(n=n*n_v,shape = shape_v,rate = rate_v),ncol=n_v,byrow=TRUE)
  make_r(drifts=drifts, n=n, b=b,A=A, n_v = n_v, t0=t0, st0=st0)
}


####### Frechet:

#' @rdname LBA
#' @export dlba_frechet
dlba_frechet <- function(t,A,b,t0,shape_v, scale_v) {

  check_vector(t, A, b, t0, shape_v, scale_v)
  nn <- length(t)
  A <- rep(A, length.out = nn)
  b <- rep(b, length.out = nn)
  t0 <- rep(t0, length.out = nn)
  shape_v <- rep(shape_v, length.out = nn)
  scale_v <- rep(scale_v, length.out = nn)
  if (any(b < A)) stop(error_message_b_smaller_A)
  
  dlba_frechet_core(t=t,A=A,b=b,t0=t0,shape_v=shape_v, scale_v=scale_v, nn=nn)
}

dlba_frechet_core <- function(t,A,b,t0,shape_v, scale_v, nn) {
  
  t <- rem_t0(t, t0)
  
  ps <- cbind(b, b-A, scale_v,shape_v)
  ps_below_zero <- apply(ps, 1, function(x) any(x <= 0))
  
  # t <- pmax(t,0) #not needed, see rem_t0
  t_old <- t
  
  t <- t[!ps_below_zero]
  A <- A[!ps_below_zero]
  b <- b[!ps_below_zero]
  t0 <- t0[!ps_below_zero]
  shape_v <- shape_v[!ps_below_zero]
  scale_v <- scale_v[!ps_below_zero]
  
  min <- (b-A)/t
  max <- b/t
  Gmax <- pfrechet(max, loc=0, scale=scale_v, shape=shape_v)
  Gmin <- pfrechet(min, loc=0, scale=scale_v, shape=shape_v)
  D <- Gmax - Gmin
  gam <- gamma_inc(1-(1/shape_v), (1/scale_v*max)^(-shape_v))-gamma_inc(1-(1/shape_v), (1/scale_v*min)^(-shape_v))
  zfrechet <- gam/(1/scale_v*D)
  diffG1 <- ((-b/(t^2))*dfrechet(b/t, loc=0, scale=scale_v, shape=shape_v))
  diffG2 <- ((-(b-A)/(t^2))*dfrechet((b-A)/t, loc=0, scale=scale_v, shape=shape_v))    
  diffD <- diffG1 - diffG2    
  diffgam <- (-shape_v*(((1/scale_v*b)^(-shape_v+1))/(t^(-shape_v+2)))*exp(-(1/scale_v*b/t)^(-shape_v))) - (-shape_v*(((1/scale_v*(b-A))^(-shape_v+1))/(t^(-shape_v+2)))*exp(-(1/scale_v*(b-A)/t)^(-shape_v)))
  diffzfrechet <- ((1/scale_v)^(-1))*(((-D^(-2))*diffD)*gam + (diffgam*(D^(-1))))
  term1 <- (Gmax - Gmin)*(zfrechet + (t*diffzfrechet))
  term2 <- diffG1*((zfrechet*t)-b)
  term3 <- diffG2*(b-A-(zfrechet*t))
  out.value <- ((term1+term2+term3)/A)
  out.value[!is.finite(out.value)] <- 0 # Set NaN or -Inf or Inf to pdf=0
  
  out <- numeric(nn)
  out[!ps_below_zero] <- out.value
  return(pmax(0, out))
}

#' @rdname LBA
#' @export plba_frechet
plba_frechet <- function(t,A,b,t0,shape_v, scale_v) {
  check_vector(t, A, b, t0, shape_v, scale_v)
  nn <- length(t)
  A <- rep(A, length.out = nn)
  b <- rep(b, length.out = nn)
  t0 <- rep(t0, length.out = nn)
  shape_v <- rep(shape_v, length.out = nn)
  scale_v <- rep(scale_v, length.out = nn)
  if (any(b < A)) stop(error_message_b_smaller_A)
  
  plba_frechet_core(t=t,A=A,b=b,t0=t0,shape_v=shape_v, scale_v=scale_v, nn=nn)
}

plba_frechet_core <- function(t,A,b,t0,shape_v, scale_v, nn) {  
  t <- rem_t0(t, t0)
  
  ps <- cbind(b, b-A, scale_v,shape_v)
  ps_below_zero <- apply(ps, 1, function(x) any(x <= 0))
  
  # t <- pmax(t,0) #not needed, see rem_t0
  t_old <- t
  
  t <- t[!ps_below_zero]
  A <- A[!ps_below_zero]
  b <- b[!ps_below_zero]
  t0 <- t0[!ps_below_zero]
  shape_v <- shape_v[!ps_below_zero]
  scale_v <- scale_v[!ps_below_zero]
  
  # t <- pmax(t,0) #not needed, see rem_t0
  min <- (b-A)/t
  max <- b/t
  pmax <- pfrechet(max, loc=0, scale=scale_v, shape=shape_v)
  pmin <- pfrechet(min, loc=0, scale=scale_v, shape=shape_v)
  zfrechet <- (gamma_inc(1-(1/shape_v),(1/scale_v*max)^(-shape_v))-gamma_inc(1-(1/shape_v),(1/scale_v*min)^(-shape_v)))/(1/scale_v*(pmax-pmin))    
  term1 <- ((t*zfrechet) - b)/A
  term2 <- (b-A-(t*zfrechet))/A 
  out.value <- (1 + pmax*term1 + pmin*term2)
  out.value[t==Inf] <- 1 # term1=Inf and term2=-Inf cancel in this case
  out.value[!is.finite(out.value)] <- 0 # Set NaN or -Inf to CDF=0
  
  out <- numeric(nn)
  out[!ps_below_zero] <- out.value
  return(pmin(pmax(0, out), 1))
}


#' @rdname LBA
#' @export rlba_frechet
rlba_frechet <- function(n,A,b,t0,shape_v, scale_v,st0=0){
  check_single_arg(n, A, b, t0, st0)
  if (b < A) stop(error_message_b_smaller_A)
  n_v <- max(length(shape_v), length(scale_v))
  drifts <- matrix(rfrechet(n=n*n_v, loc=0, scale=scale_v, shape=shape_v),ncol=n_v,byrow=TRUE)
  
  make_r(drifts=drifts, n=n, b=b,A=A, n_v=n_v, t0=t0, st0=st0)
}


####### Log-Normal:

#' @rdname LBA
#' @export dlba_lnorm
dlba_lnorm <- function(t,A,b,t0,meanlog_v, sdlog_v, robust = FALSE) {
  check_vector(t, A, b, t0, meanlog_v, sdlog_v)
  nn <- length(t)
  A <- rep(A, length.out = nn)
  b <- rep(b, length.out = nn)
  t0 <- rep(t0, length.out = nn)
  mean_v <- rep(meanlog_v, length.out = nn)
  sd_v <- rep(sdlog_v, length.out = nn)
  if (any(b < A)) stop(error_message_b_smaller_A)
  
  dlba_lnorm_core(t=t,A=A,b=b,t0=t0,meanlog_v=meanlog_v, sdlog_v=sdlog_v, robust = robust, nn=nn)
}

dlba_lnorm_core <- function(t,A,b,t0,meanlog_v, sdlog_v, robust=FALSE, nn) {
  if (robust) { # robust == TRUE uses robust versions of the normal distributions
    pnorm1 <- pnormP
    dnorm1 <- dnormP
  } else {
    pnorm1 <- pnorm
    dnorm1 <- dnorm  
  }
  
  t <- rem_t0(t, t0)
  
  min <- (b-A)/t
  max <- b/t
  
  zlognorm <- (exp(meanlog_v+(sdlog_v^2)/2)*(pnorm1((log(max)-meanlog_v-(sdlog_v^2))/sdlog_v)-pnorm1((log(min)-meanlog_v-(sdlog_v^2))/sdlog_v))) / (pnorm1((log(max)-meanlog_v)/sdlog_v)-pnorm1((log(min)-meanlog_v)/sdlog_v))
  Gmax <- plnorm(max,meanlog=meanlog_v,sdlog=sdlog_v) 
  Gmin <- plnorm(min,meanlog=meanlog_v,sdlog=sdlog_v)
  
  u <- (pnorm1((log(max)-meanlog_v-(sdlog_v)^2)/sdlog_v)-pnorm1((log(min)-meanlog_v-(sdlog_v)^2)/sdlog_v))
  v <- (pnorm1((log(max)-meanlog_v)/sdlog_v)-pnorm1((log(min)-meanlog_v)/sdlog_v))
  
  udash <- (((-1/(sdlog_v*t))*dnorm1((log(b/t)-meanlog_v-(sdlog_v)^2)/sdlog_v)) - ((-1/(sdlog_v*t))*dnorm1((log((b-A)/t)-meanlog_v-(sdlog_v)^2)/sdlog_v)))
  vdash <- (((-1/(sdlog_v*t))*dnorm1((log(b/t)-meanlog_v)/sdlog_v)) - ((-1/(sdlog_v*t))*dnorm1((log((b-A)/t)-meanlog_v)/sdlog_v)))
  const <- exp(meanlog_v+((sdlog_v)^2)/2)
  
  diffzlognorm <- ((udash*v - vdash*u)/(v^2))*const #quotient rule
  term1 <- (Gmax - Gmin)*(zlognorm + (t*diffzlognorm))
  term2 <- ((-b/(t^2))*dlnorm(b/t,meanlog=meanlog_v,sdlog=sdlog_v))*((zlognorm*t)-b)
  term3 <- (b-A-(zlognorm*t))*((-(b-A)/(t^2))*dlnorm((b-A)/t,meanlog=meanlog_v,sdlog=sdlog_v))
  out.value <- ((term1+term2+term3)/A)
  out.value[!is.finite(out.value)] <- 0 # Set NaN or -Inf or Inf to pdf=0
  return(pmax(0, out.value))
}

#' @rdname LBA
#' @export plba_lnorm
plba_lnorm <- function(t,A,b,t0,meanlog_v, sdlog_v, robust = FALSE) {
  check_vector(t, A, b, t0, meanlog_v, sdlog_v)
  nn <- length(t)
  A <- rep(A, length.out = nn)
  b <- rep(b, length.out = nn)
  t0 <- rep(t0, length.out = nn)
  mean_v <- rep(meanlog_v, length.out = nn)
  sd_v <- rep(sdlog_v, length.out = nn)
  if (any(b < A)) stop(error_message_b_smaller_A)
  
  plba_lnorm_core(t=t,A=A,b=b,t0=t0,meanlog_v=meanlog_v, sdlog_v=sdlog_v, robust=robust, nn=nn)
}

plba_lnorm_core <- function(t,A,b,t0,meanlog_v, sdlog_v, robust = FALSE, nn) {
  if (robust) { # robust == TRUE uses robust versions of the normal distributions
    pnorm1 <- pnormP
  } else {
    pnorm1 <- pnorm 
  }
  
  t <- rem_t0(t, t0)
  min <- (b-A)/t
  max <- b/t
  zlognorm <- (exp(meanlog_v+(sdlog_v^2)/2)*(pnorm1((log(max)-meanlog_v-(sdlog_v^2))/sdlog_v)-pnorm1((log(min)-meanlog_v-(sdlog_v^2))/sdlog_v))) / (pnorm1((log(max)-meanlog_v)/sdlog_v)-pnorm1((log(min)-meanlog_v)/sdlog_v))
  term1 <- ((t*zlognorm) - b)/A
  term2 <- (b-A-(t*zlognorm))/A 
  pmax <- plnorm(max, meanlog=meanlog_v, sdlog=sdlog_v) 
  pmin <- plnorm(min, meanlog=meanlog_v, sdlog=sdlog_v)
  out.value <- (1 + pmax*term1 + pmin*term2)
  out.value[t==Inf] <- 1 # term1=Inf and term2=-Inf cancel in this case
  out.value[!is.finite(out.value)] <- 0 # Set NaN or -Inf to CDF=0
  return(pmin(pmax(0, out.value), 1))
}

#' @rdname LBA
#' @export rlba_lnorm
rlba_lnorm <- function(n,A,b,t0,meanlog_v, sdlog_v, st0=0){
  check_single_arg(n, A, b, t0, st0)
  if (b < A) stop(error_message_b_smaller_A)
  n_v <- max(length(meanlog_v), length(sdlog_v))
  drifts=matrix(rlnorm(n=n*n_v,meanlog = meanlog_v,sdlog=sdlog_v),ncol=n_v,byrow=TRUE)
  make_r(drifts=drifts, n=n, b=b, A=A, n_v=n_v, t0=t0, st0=st0)
}

