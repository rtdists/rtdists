#' LBA race functions: Likelihood for first accumulator to win.
#' 
#' n1PDF and n1CDF take RTs, the distribution functions of the \link{LBA}, and corresponding parameter values and put them throughout the race equations and return the likelihood for the first accumulator winning (hence n1) in a set of accumulators.  
#'
#' @param t a vector of RTs.
#' @param A,b LBA parameters, see \code{\link{LBA}}. Can either be a single numeric value or vector (which will be recycled to reach \code{length(t)}) \emph{or} a \code{list} of such vectors in which each list element corresponds to the parameters for this accumulator (i.e., the list needs to be of the same length as there are accumulators).
#' @param t0 \emph{one} scalar \code{t0} parameter (see \code{\link{LBA}}). Multiple \code{t0} parameters are currently not implemented.
#' @param st0 \emph{one} scalar parameter specifying the variability of \code{t0} (which varies uniformly from \code{t0} to \code{t0} + \code{st0}).
#' @param ... two \emph{named} drift rate parameters dependening on \code{distribution} (e.g., \code{mean_v} and \code{sd_v} for \code{distribution=="norm"}). 
#' @param distribution character specifying the distribution of the drift rate. Possible values are \code{c("norm", "gamma", "frechet", "lnorm")}, default is \code{"norm"}.
#' @param args.dist list of optional further arguments to the distribution functions (i.e., \code{posdrift} or \code{robust} for \code{distribution=="norm"}).
#' 
#' 
#' @details For a set of \eqn{N} independent accumulators \eqn{i = 1...N}, the race likelihood for a given accumulator \eqn{i} is given by
#' \deqn{L(\mbox{unit }i \mbox{ wins}) = f_i(t) \times \prod_j^i [ S_j(t) ]}{L(unit i wins) = f_i(t) * prod_j<>i [ S_j(t) ]}
#' where \eqn{f(t)} is the PDF (\code{dlba_...}) and \eqn{S_j(t) = 1 - F_j(t)} is the survivor function, that is the complement of the CDF \eqn{F(t)} (\code{plba_...}) at time \eqn{t}.
#' 
#' In other words, this is just the PDF/CDF for the winning accumulator at time \eqn{t} times the probability that no other accumulators have finished at time \eqn{t}.
#' 
#' @name LBA-race
#' 
#' @example examples/examples.lba-race.R
#' 
NULL

## note, this functions does not check parameters, it is only< called internally (i.e., passed correctly).
n1PDFfixedt0 <- function(t,A,b, t0, ..., distribution, args.dist = list()) {
  # Generates defective PDF for responses on node #1.
  dots <- list(...)
  #if (is.null(names(dots))) stop("... arguments need to be named.")
  #check_single_arg(t0 = t0)
  #distribution <- match.arg(distribution)
  #browser()
  switch(distribution, 
         norm = {
           pdf <- dlba_norm
           cdf <- plba_norm
           if (any(!(c("mean_v","sd_v") %in% names(dots)))) stop("mean_v and sd_v need to be passed for distribution = \"norm\"")
         },
         gamma = {
           pdf <- dlba_gamma
           cdf <- plba_gamma
           if (!("shape_v" %in% names(dots))) stop("shape_v needs to be passed for distribution = \"gamma\"")
           if ((!("rate_v" %in% names(dots))) & (!("scale_v" %in% names(dots)))) stop("rate_v or scale_v needs to be passed for distribution = \"gamma\"")
         },
         frechet = {
           pdf <- dlba_frechet
           cdf <- plba_frechet
           if (any(!(c("shape_v","scale_v") %in% names(dots)))) stop("shape_v and scale_v need to be passed for distribution = \"frechet\"")
         },
         lnorm = {
           pdf <- dlba_lnorm
           cdf <- plba_lnorm
           if (any(!(c("meanlog_v","sdlog_v") %in% names(dots)))) stop("meanlog_v and sdlog_v need to be passed for distribution = \"lnorm\"")
         }
  )
  n_v <- max(vapply(dots, length, 0))  # Number of responses
  if (n_v>2) {
    tmp=array(dim=c(length(t),n_v-1))
    for (i in 2:n_v) tmp[,i-1] <- do.call(cdf, args = c(t=list(t), A=if(is.list(A)) A[i] else list(A), b=if(is.list(b)) b[i] else list(b), t0 = t0, sapply(dots, "[[", i = i, simplify = FALSE), args.dist = args.dist))
    G <- apply(1-tmp,1,prod)
  } else {
    G <- 1-do.call(cdf, args = c(t=list(t), A=if(is.list(A)) A[2] else list(A), b=if(is.list(b)) b[2] else list(b), t0 = t0, sapply(dots, "[[", i = 2, simplify = FALSE), args.dist = args.dist))
  }
  G*do.call(pdf, args = c(t=list(t), A=if(is.list(A)) A[1] else list(A), b=if(is.list(b)) b[1] else list(b), t0 = t0, sapply(dots, "[[", i = 1, simplify = FALSE), args.dist = args.dist))
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
  nn <- length(t)
  
  if (!is.list(A)) {
    A <- rep(A, length.out=nn)
  } else {
    if (length(A) != n_v) stop("if A is a list, its length needs to correspond to the number of accumulators.")
    for (i in seq_along(A)) {
      A[[i]] <- rep(A[[i]], length.out=nn)
    }
  }
  if (!is.list(b)) {
    b <- rep(b, length.out=nn)
  } else {
    if (length(b) != n_v) stop("if b is a list, its length needs to correspond to the number of accumulators.")
    for (i in seq_along(b)) {
      b[[i]] <- rep(b[[i]], length.out=nn)
    }
  }
  #if (is.null(dim(A)) && length(A)== 1) A <- rep(A,length.out=n_v)
  #if (length(b)<n_v) b <- rep(b,length.out=n_v)
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
    for (i in 1:length(outs)) {
      tmp <- do.call(integrate, args=c(f=tmpf, lower=t[i]-t0[1]-st0, upper=t[i]-t0[1], A=list(A), b=list(b), t0=list(0), dots, distribution = distribution, args.dist = args.dist, stop.on.error = FALSE))
      if (tmp$message != "OK") warning(paste("n1PDF:", tmp$message))
      outs[i] <- tmp$value
    }
    return(outs)
  }
}

# t = time, A=x0max, b=chi, v=drift, sv=sdI
#' @rdname LBA-race
#' @export
n1CDF <- function(t,A,b, t0, ..., st0=0, distribution = c("norm", "gamma", "frechet", "lnorm"), args.dist = list()) {  #, browser=FALSE
  # Generates defective CDF for responses on node #1. 
  dots <- list(...)
  if (is.null(names(dots))) stop("... arguments need to be named.")
  check_single_arg(t0 = t0)
  n_v <- max(vapply(dots, length, 0))  # Number of responses
  if (n_v < 2) stop("There need to be at least two accumulators/drift rates.")
  distribution <- match.arg(distribution)
  nn <- length(t)
  if (!is.list(A)) {
    A <- rep(A, length.out=nn)
  } else {
    if (length(A) != n_v) stop("if A is a list, its length needs to correspond to the number of accumulators.")
    for (i in seq_along(A)) {
      A[[i]] <- rep(A[[i]], length.out=nn)
    }
  }
  if (!is.list(b)) {
    b <- rep(b, length.out=nn)
  } else {
    if (length(b) != n_v) stop("if b is a list, its length needs to correspond to the number of accumulators.")
    for (i in seq_along(b)) {
      b[[i]] <- rep(b[[i]], length.out=nn)
    }
  }
  
  #if (length(t0)<n_v) t0 <- rep(t0,length.out=n_v)
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
      #if (i==1 && browser) browser()
      tmp_obj <- do.call(integrate, args=c(f=n1PDF,lower=bounds[i],upper=bounds[i+1],subdivisions=1000, A=list(A), b=list(b), t0 = list(t0), st0 = list(st0), dots, distribution = distribution, stop.on.error = FALSE, args.dist = args.dist))
      #browser()
      if (tmp_obj$message != "OK") {
        #browser()
        warning(tmp_obj$message)
      }
      tmp <- tmp_obj$value
      
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


