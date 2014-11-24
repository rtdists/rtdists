#' LBA race functions
#' 
#' as described elsewhere.
#'
#' @param t a vector of RTs.
#' @param A,b,t0 LBA parameters, see \code{\link{LBA}}.
#' @param ... two named drift rate parameters dependening on \code{distribution} (e.g., \code{mean_v} and \code{sd_v} for \code{distribution=="norm"}). 
#' @param distribution character specifying the distribution of the drift rate.
#' @param args.dist Optional further arguments to the distribution functions (i.e., \code{posdrift} or \code{robust} for \code{distribution=="norm"}).
#' @param st0 variability of \code{t0}.
#' 
#' 
#' @name LBA-race
#' 
#' @example examples/examples.lba-race.R
#' 
NULL

# t = time, A=x0max, b=chi, v=drift, sv=sdI
n1PDFfixedt0 <- function(t,A,b, t0, ..., distribution = c("norm", "gamma", "frechet", "lnorm"), args.dist = list()) {
  # Generates defective PDF for responses on node #1.
  dots <- list(...)
  if (is.null(names(dots))) stop("... arguments need to be named.")
  check_single_arg(t0 = t0)
  distribution <- match.arg(distribution)
  switch(distribution, 
         norm = {
           pdf <- dlba_norm
           cdf <- plba_norm
         },
         gamma = {
           pdf <- dlba_gamma
           cdf <- plba_gamma
         },
         frechet = {
           pdf <- dlba_frechet
           cdf <- plba_frechet
         },
         lnorm = {
           pdf <- dlba_lnorm
           cdf <- plba_lnorm
         }
         )
  n_v <- max(vapply(dots, length, 0))  # Number of responses
  if (n_v>2) {
    tmp=array(dim=c(length(t),n_v-1))
    for (i in 2:n_v) tmp[,i-1] <- do.call(cdf, args = c(t=list(t), A=A[i], b=b[i], t0 = t0, sapply(dots, "[[", i = i, simplify = FALSE), args.dist = args.dist))
    G <- apply(1-tmp,1,prod)
  } else {
    G <- 1-do.call(cdf, args = c(t=list(t), A=A[2], b=b[2], t0 = t0, sapply(dots, "[[", i = 2, simplify = FALSE), args.dist = args.dist))
  }
  G*do.call(pdf, args = c(t=list(t), A=A[1], b=b[1], t0 = t0, sapply(dots, "[[", i = 1, simplify = FALSE), args.dist = args.dist))
}

#' @rdname LBA-race
#' @export
n1PDF <- function(t, A, b, t0, ..., st0=0, distribution = c("norm", "gamma", "frechet", "lnorm"), args.dist = list()) {
  dots <- list(...)
  if (is.null(names(dots))) stop("... arguments need to be named.")
  n_v <- max(vapply(dots, length, 0))  # Number of responses
  if (n_v < 2) stop("There need to be at least two accumulators/drift rates.")
  distribution <- match.arg(distribution)
  check_single_arg(t0 = t0)
  if (length(A)<n_v) A <- rep(A,length.out=n_v)
  if (length(b)<n_v) b <- rep(b,length.out=n_v)
  #if (length(t0)<n_v) t0 <- rep(t0,length.out=n_v)
  for (i in length(dots)) {
    if (length(dots[[i]]) < n_v) dots[[i]] <- rep(dots[[i]],length.out=n_v)
  }
  if (length(st0)>1) {
    warning("st0 set to st0[1]. Only one non-decision time variability permitted.")
    st0 <- st0[1] # Only ONE non-decision time.
  }
  #browser()
  if (st0==0) return(do.call(n1PDFfixedt0, args = c(t=list(t), A=list(A), b=list(b), t0 = list(t0), dots, distribution = distribution, args.dist = args.dist)))
  else {
    tmpf <- function(t, A, b, t0, ..., plba, plba.args = list()) {
      #browser()
      do.call(n1PDFfixedt0, args = c(t=list(pmax(t-t0, 0)), A=list(A), t0 = list(0), b=list(b), dots, distribution = distribution, args.dist = args.dist))/st0
    }
    outs=numeric(length(t))
    #browser()
    for (i in 1:length(outs))
      outs[i] <- do.call(integrate, args=c(f=tmpf, lower=t[i]-t0[1]-st0, upper=t[i]-t0[1], A=list(A), b=list(b), t0=list(0), dots, distribution = distribution, args.dist = args.dist))$value
    return(outs)
  }
}

# t = time, A=x0max, b=chi, v=drift, sv=sdI
#' @rdname LBA-race
#' @export
n1CDF <- function(t,A,b, t0, ..., st0=0, distribution = c("norm", "gamma", "frechet", "lnorm"), args.dist = list()) {
  # Generates defective CDF for responses on node #1. 
  dots <- list(...)
  if (is.null(names(dots))) stop("... arguments need to be named.")
  check_single_arg(t0 = t0)
  n_v <- max(vapply(dots, length, 0))  # Number of responses
  if (n_v < 2) stop("There need to be at least two accumulators/drift rates.")
  distribution <- match.arg(distribution)
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
  if (st0<1e-6) {
    if(!isTRUE(all.equal(0, st0))) warning("st0 set to 0. Integral can fail for small st0.")
    st0=0
    } # 
  outs <- numeric(length(t))
  bounds <- c(0,t)
  for (i in 1:length(t)) {
    tmp <- "error"
    repeat {
      if (bounds[i]>=bounds[i+1]) {
        outs[i]=0
        break
      }
      tmp <- try(do.call(integrate, args=c(f=n1PDF,lower=bounds[i],upper=bounds[i+1],subdivisions=1000, A=list(A), b=list(b), t0 = list(t0[1]), st0 = list(st0), dots, distribution = distribution, args.dist = args.dist))$value,silent=T)
      
      if (is.numeric(tmp)) {
        outs[i]=tmp
        break
      }
      # Try smart lower bound.
      if ((distribution == "norm") && (bounds[i]<=0)) {
        #browser()
        bounds[i] <- max(c((b-0.98*A)/(max(mean(dots$mean_v),dots$mean_v[1])+2*dots$sd_v)[1],0))
        next
      }
      # Try smart upper bound.
      if ((distribution == "norm") && (bounds[i+1]==Inf)) {
        bounds[i+1]=0.02*max(b)/(mean(dots$mean_v)-2*mean(dots$sd_v))
        next
      }
      stop("Error in n1CDF that I could not catch. Please report (with data used) to maintainer email.")
    }
  }
  cumsum(outs)
}

