#' LBA race functions: Likelihood for first accumulator to win.
#' 
#' n1PDF and n1CDF take RTs, the distribution functions of the \link{LBA}, and corresponding parameter values and put them throughout the race equations and return the likelihood for the first accumulator winning (hence n1) in a set of accumulators.  
#'
#' @param rt a vector of RTs.
#' @param A,b,t0 LBA parameters, see \code{\link{LBA}}. Can either be a single numeric vector (which will be recycled to reach \code{length(rt)} for trialwise parameters) \emph{or} a \code{list} of such vectors in which each list element corresponds to the parameters for this accumulator (i.e., the list needs to be of the same length as there are accumulators). Each list will also be recycled to reach \code{length(rt)} for trialwise parameters per accumulator.
#' @param st0 parameter specifying the variability of \code{t0} (which varies uniformly from \code{t0} to \code{t0} + \code{st0}). Can be trialwise, and will be recycled to length of \code{rt}.
#' @param ... two \emph{named} drift rate parameters depending on \code{distribution} (e.g., \code{mean_v} and \code{sd_v} for \code{distribution=="norm"}). The parameters can either be given as a numeric vector or a list. If a numeric vector is passed each element of the vector corresponds to one accumulator. If a list is passed each list element corresponds to one accumulator allowing again trialwise driftrates. The shorter parameter will be recycled as necessary (and also the elements of the list to match the length of \code{rt}). See examples.
#' @param distribution character specifying the distribution of the drift rate. Possible values are \code{c("norm", "gamma", "frechet", "lnorm")}, default is \code{"norm"}.
#' @param args.dist list of optional further arguments to the distribution functions (i.e., \code{posdrift} or \code{robust} for \code{distribution=="norm"}).
#' @param silent logical. Should the number of accumulators used be suppressed? Default is \code{FALSE} which prints the number of accumulators.
#' 
#' 
#' @details For a set of \eqn{N} independent accumulators \eqn{i = 1...N}, the race likelihood for a given accumulator \eqn{i} is given by
#' \deqn{L(\mbox{unit }i \mbox{ wins}) = f_i(t) \times \prod_{j<>i} [ S_j(t) ]}{L(unit i wins) = f_i(t) * prod_j<>i [ S_j(t) ]}
#' where \eqn{f(t)} is the PDF (\code{dlba_...}) and \eqn{S_j(t) = 1 - F_j(t)} is the survivor function, that is the complement of the CDF \eqn{F(t)} (\code{plba_...}) at time \eqn{t}.
#' 
#' In other words, this is just the PDF/CDF for the winning accumulator at time \eqn{t} times the probability that no other accumulators have finished at time \eqn{t}.
#' 
#' @seealso For more user-friendly functions that return the PDF or CDF for the corresponding (and not first) accumulator winning see /code{/link{LBA}}.
#' 
#' @name LBA-race
#' @importFrom stats integrate
#' 
#' @example examples/examples.lba-race.R
#' 
NULL

## note, this functions does not check parameters, it is only called internally (i.e., passed correctly).
n1PDFfixedt0 <- function(rt,A,b, t0, ..., pdf, cdf, args.dist = list()) {
  # Generates defective PDF for responses on node #1.
  dots <- list(...)
  nn <- length(rt)
  #if (length(A) != nn) browser()
  #browser()
  n_v <- max(vapply(dots, length, 0))  # Number of responses
  if (n_v>2) {
    tmp=array(dim=c(length(rt),n_v-1))
    for (i in 2:n_v) tmp[,i-1] <- do.call(cdf, args = c(rt=list(rt), A=if(is.list(A)) A[i] else list(A), b=if(is.list(b)) b[i] else list(b), t0 = if(is.list(t0)) t0[i] else list(t0), sapply(dots, "[[", i = i, simplify = FALSE), args.dist = args.dist, nn=nn))
    G <- apply(1-tmp,1,prod)
  } else {
    G <- 1-do.call(cdf, args = c(rt=list(rt), A=if(is.list(A)) A[2] else list(A), b=if(is.list(b)) b[2] else list(b), t0 = if(is.list(t0)) t0[2] else list(t0), sapply(dots, "[[", i = 2, simplify = FALSE), args.dist, nn=nn))
  }
  G*do.call(pdf, args = c(rt=list(rt), A=if(is.list(A)) A[1] else list(A), b=if(is.list(b)) b[1] else list(b), t0 = if(is.list(t0)) t0[1] else list(t0), sapply(dots, "[[", i = 1, simplify = FALSE), args.dist, nn=nn))
}

#sapply(dots, rep_dots, which = 1, nn = nn, simplify = FALSE)
rep_dots <- function(arg, which, nn) {
  rep(arg[[which]], length.out=nn)
}

## functions which checks if argument is numeric and 
check_n1_arguments <- function(arg, nn, n_v, dots = FALSE) {
  mc <- match.call()
  varname <- sub("dots$", "", deparse(mc[["arg"]]), fixed = TRUE)
  if (!is.list(arg)) {
    if ((!is.vector(arg, "numeric")) || (length(arg) < 1)) stop(paste(varname, "needs to be a numeric vector of length >= 1!"))
    if (dots) {
      arg <- as.list(arg)
      arg <- lapply(arg, rep, length.out=nn)
    } else arg <- rep(arg, length.out=nn)
  } else {
    if (!dots && (length(arg) != n_v)) stop(paste("if", varname, "is a list, its length needs to correspond to the number of accumulators."))
    for (i in seq_along(arg)) {
      if ((!is.vector(arg[[i]], "numeric")) || (length(arg[[i]]) < 1)) stop(paste0(varname, "[[", i, "]] needs to be a numeric vector of length >= 1!"))
      arg[[i]] <- rep(arg[[i]], length.out=nn)
    }
  }
  return(arg)
}

#' @rdname LBA-race
#' @export
n1PDF <- function(rt, A, b, t0, ..., st0=0, distribution = c("norm", "gamma", "frechet", "lnorm"), args.dist = list(), silent = FALSE) {
  dots <- list(...)
  #browser()
  if (is.null(names(dots))) stop("... arguments need to be named.")
  if (any(names(dots) == "")) stop("all ... arguments need to be named.")
  
  n_v <- max(vapply(dots, length, 0))  # Number of responses
  if(!silent) message(paste("Results based on", n_v, "accumulators/drift rates."))
  if (n_v < 2) stop("There need to be at least two accumulators/drift rates.")
  distribution <- match.arg(distribution)
  #check_single_arg(t0 = t0)
  nn <- length(rt)
  #browser()
  A <- check_n1_arguments(A, nn=nn, n_v=n_v)
  b <- check_n1_arguments(b, nn=nn, n_v=n_v)
  t0 <- check_n1_arguments(t0, nn=nn, n_v=n_v)
  st0 <- rep(unname(st0), length.out = nn)
  switch(distribution, 
         norm = {
           pdf <- dlba_norm_core
           cdf <- plba_norm_core
           if (any(!(c("mean_v","sd_v") %in% names(dots)))) stop("mean_v and sd_v need to be passed for distribution = \"norm\"")
           dots$mean_v <- check_n1_arguments(dots$mean_v, nn=nn, n_v=n_v, dots = TRUE)
           dots$sd_v <- check_n1_arguments(dots$sd_v, nn=nn, n_v=n_v, dots = TRUE)
           dots <- dots[c("mean_v","sd_v")]
         },
         gamma = {
           pdf <- dlba_gamma_core
           cdf <- plba_gamma_core
           if (!("shape_v" %in% names(dots))) stop("shape_v needs to be passed for distribution = \"gamma\"")
           if ((!("rate_v" %in% names(dots))) & (!("scale_v" %in% names(dots)))) stop("rate_v or scale_v need to be passed for distribution = \"gamma\"")
           dots$shape_v <- check_n1_arguments(dots$shape_v, nn=nn, n_v=n_v, dots = TRUE)
           if ("scale_v" %in% names(dots)) {
             dots$scale_v <- check_n1_arguments(dots$scale_v, nn=nn, n_v=n_v, dots = TRUE)
             if (is.list(dots$scale_v)) {
               dots$rate_v <- lapply(dots$scale_v, function(x) 1/x)
             } else dots$rate_v <- 1/dots$scale_v
           } else dots$rate_v <- check_n1_arguments(dots$rate_v, nn=nn, n_v=n_v, dots = TRUE)
           dots <- dots[c("shape_v","rate_v")]
         },
         frechet = {
           pdf <- dlba_frechet_core
           cdf <- plba_frechet_core
           if (any(!(c("shape_v","scale_v") %in% names(dots)))) stop("shape_v and scale_v need to be passed for distribution = \"frechet\"")
           dots$shape_v <- check_n1_arguments(dots$shape_v, nn=nn, n_v=n_v, dots = TRUE)
           dots$scale_v <- check_n1_arguments(dots$scale_v, nn=nn, n_v=n_v, dots = TRUE)
           dots <- dots[c("shape_v","scale_v")]
         },
         lnorm = {
           pdf <- dlba_lnorm_core
           cdf <- plba_lnorm_core
           if (any(!(c("meanlog_v","sdlog_v") %in% names(dots)))) stop("meanlog_v and sdlog_v need to be passed for distribution = \"lnorm\"")
           dots$meanlog_v <- check_n1_arguments(dots$meanlog_v, nn=nn, n_v=n_v, dots = TRUE)
           dots$sdlog_v <- check_n1_arguments(dots$sdlog_v, nn=nn, n_v=n_v, dots = TRUE)
           dots <- dots[c("meanlog_v","sdlog_v")]
         }
  )
  #browser()
  for (i in seq_len(length(dots))) {
    if (length(dots[[i]]) < n_v) dots[[i]] <- rep(dots[[i]],length.out=n_v)
  }
#   if (length(st0)>1) {
#     warning("st0 set to st0[1]. Only one non-decision time variability permitted.")
#     st0 <- st0[1] # Only ONE non-decision time.
#   }
  #browser()
  do.call(n1PDF_core, args = c(rt=list(rt), A=list(A), b=list(b), t0 = list(t0), st0 = list(st0), dots, pdf=pdf, cdf=cdf, args.dist = list(args.dist)))
}


n1PDF_core <- function(rt, A, b, t0, ..., st0, pdf, cdf, args.dist = list()) {
  dots <- list(...)
  #browser()
  if (all(st0==0)) return(do.call(n1PDFfixedt0, args = c(rt=list(rt), A=list(A), b=list(b), t0 = list(t0), dots, pdf=pdf, cdf=cdf, args.dist = list(args.dist))))
  else {
    tmpf <- function(rt, A, b, t0, st0, ..., pdf, cdf, args.dist = list()) {
      #browser()
      dots2 <- list(...)
      do.call(n1PDFfixedt0, args = c(rt=list(rt), A=list(A), b=list(b), t0 = list(t0), dots2, pdf=pdf, cdf=cdf, args.dist = list(args.dist)))/st0
      #rt=list(pmax(rt-t0, 0))
    }
    outs <- vector("numeric", length = length(rt))
    if (length(st0) == 1) st0 <- rep(st0, length.out = length(rt))
    for (i in 1:length(rt)) {
      if (st0[i] != 0) {
        tmp <- do.call(integrate, args=c(f=tmpf, lower=unname(rt[i]-st0[i]), upper=unname(rt[i]), A=ret_arg(A, i), b=ret_arg(b, i), t0=ret_arg(t0, i), sapply(dots, function(z, i) sapply(z, ret_arg2, which = i, simplify=FALSE), i=i, simplify=FALSE), pdf=pdf, cdf=cdf, args.dist = list(args.dist), stop.on.error = FALSE, st0 = list(st0[i])))
        if (tmp$message != "OK") warning(paste("n1PDF:", tmp$message))
        outs[i] <- tmp$value
      } else outs[i] <- do.call(n1PDFfixedt0, args = c(rt=list(rt[i]), A=ret_arg(A, i), b=ret_arg(b, i), t0=ret_arg(t0, i), sapply(dots, function(z, i) sapply(z, ret_arg2, which = i, simplify=FALSE), i=i, simplify=FALSE), pdf=pdf, cdf=cdf, args.dist = list(args.dist)))
    }
    return(outs)
  }
}

# n1PDF_single <- function(rt, A, b, t0, ..., st0, pdf, cdf, args.dist = list()) {
#   dots <- list(...)
#   #browser()
#   if (st0==0) return(do.call(n1PDFfixedt0, args = c(rt=list(rt), A=list(A), b=list(b), t0 = list(t0), dots, pdf=pdf, cdf=cdf, args.dist = args.dist)))
#   else {
#     tmpf <- function(rt, A, b, t0, ..., pdf, cdf, args.dist = list()) {
#       #browser()
#       do.call(n1PDFfixedt0, args = c(rt=list(pmax(rt-t0, 0)), A=list(A), b=list(b), t0 = list(0), dots, pdf=pdf, cdf=cdf, args.dist = args.dist))/st0
#     }
#     outs=numeric(length(rt))
#     #browser()
#     for (i in 1:length(outs)) {
#       tmp <- do.call(integrate, args=c(f=tmpf, lower=rt[i]-t0[1]-st0, upper=rt[i]-t0[1], A=list(A), b=list(b), t0=list(0), dots, pdf=pdf, cdf=cdf, args.dist = args.dist, stop.on.error = FALSE))
#       if (tmp$message != "OK") warning(paste("n1PDF:", tmp$message))
#       outs[i] <- tmp$value
#     }
#     return(outs)
#   }
# }

ret_arg <- function(arg, which) {
  list(if(is.list(arg)) {
    if (which <= min(sapply(arg, length))) sapply(arg, "[[", i = which, simplify = FALSE) else arg 
    } else {
      if (which <= length(arg)) arg[which]
      else arg
      })
}

ret_arg2 <- function(arg, which) {
  if(is.list(arg)) {
    if (which <= min(sapply(arg, length))) sapply(arg, "[[", i = which, simplify = FALSE) else arg 
    } else {
      if (which <= length(arg)) arg[which]
      else arg
      }
}

# rt = time, A=x0max, b=chi, v=drift, sv=sdI
#' @rdname LBA-race
#' @export
n1CDF <- function(rt,A,b, t0, ..., st0=0, distribution = c("norm", "gamma", "frechet", "lnorm"), args.dist = list(), silent = FALSE) {  #, browser=FALSE
  # Generates defective CDF for responses on node #1. 
  dots <- list(...)
  if (is.null(names(dots))) stop("... arguments need to be named.")
  
  n_v <- max(vapply(dots, length, 0))  # Number of responses
  if(!silent) message(paste("Results based on", n_v, "accumulators/drift rates."))
  if (n_v < 2) stop("There need to be at least two accumulators/drift rates.")
  distribution <- match.arg(distribution)
  #check_single_arg(t0 = t0)
  nn <- length(rt)
  #browser()
  A <- check_n1_arguments(A, nn=nn, n_v=n_v)
  b <- check_n1_arguments(b, nn=nn, n_v=n_v)
  t0 <- check_n1_arguments(t0, nn=nn, n_v=n_v)
  st0 <- rep(unname(st0), length.out = nn)
  switch(distribution, 
         norm = {
           pdf <- dlba_norm_core
           cdf <- plba_norm_core
           if (any(!(c("mean_v","sd_v") %in% names(dots)))) stop("mean_v and sd_v need to be passed for distribution = \"norm\"")
           dots$mean_v <- check_n1_arguments(dots$mean_v, nn=nn, n_v=n_v, dots = TRUE)
           dots$sd_v <- check_n1_arguments(dots$sd_v, nn=nn, n_v=n_v, dots = TRUE)
           dots <- dots[c("mean_v","sd_v")]
         },
         gamma = {
           pdf <- dlba_gamma_core
           cdf <- plba_gamma_core
           if (!("shape_v" %in% names(dots))) stop("shape_v needs to be passed for distribution = \"gamma\"")
           if ((!("rate_v" %in% names(dots))) & (!("scale_v" %in% names(dots)))) stop("rate_v or scale_v need to be passed for distribution = \"gamma\"")
           dots$shape_v <- check_n1_arguments(dots$shape_v, nn=nn, n_v=n_v, dots = TRUE)
           if ("scale_v" %in% names(dots)) {
             dots$scale_v <- check_n1_arguments(dots$scale_v, nn=nn, n_v=n_v, dots = TRUE)
             if (is.list(dots$scale_v)) {
               dots$rate_v <- lapply(dots$scale_v, function(x) 1/x)
             } else dots$rate_v <- 1/dots$scale_v
           } else dots$rate_v <- check_n1_arguments(dots$rate_v, nn=nn, n_v=n_v, dots = TRUE)
           dots <- dots[c("shape_v","rate_v")]
         },
         frechet = {
           pdf <- dlba_frechet_core
           cdf <- plba_frechet_core
           if (any(!(c("shape_v","scale_v") %in% names(dots)))) stop("shape_v and scale_v need to be passed for distribution = \"frechet\"")
           dots$shape_v <- check_n1_arguments(dots$shape_v, nn=nn, n_v=n_v, dots = TRUE)
           dots$scale_v <- check_n1_arguments(dots$scale_v, nn=nn, n_v=n_v, dots = TRUE)
           dots <- dots[c("shape_v","scale_v")]
         },
         lnorm = {
           pdf <- dlba_lnorm_core
           cdf <- plba_lnorm_core
           if (any(!(c("meanlog_v","sdlog_v") %in% names(dots)))) stop("meanlog_v and sdlog_v need to be passed for distribution = \"lnorm\"")
           dots$meanlog_v <- check_n1_arguments(dots$meanlog_v, nn=nn, n_v=n_v, dots = TRUE)
           dots$sdlog_v <- check_n1_arguments(dots$sdlog_v, nn=nn, n_v=n_v, dots = TRUE)
           dots <- dots[c("meanlog_v","sdlog_v")]
         }
  )
  
  for (i in seq_len(length(dots))) {
    if (length(dots[[i]]) < n_v) dots[[i]] <- rep(dots[[i]],length.out=n_v)
  }
#   if (length(st0)>1) {
#     warning("st0 set to st0[1]. Only one non-decision time variability permitted.")
#     st0 <- st0[1] # Only ONE non-decision time.
#   }
  if (any(st0<1e-6)) {
    if (any(sapply(st0[st0<1e-6], function(x) !isTRUE(all.equal(x, 0))))) warning("st0 set to 0 for values < 1e-6. Integral can fail for small st0.")
    st0[st0<1e-6] <- 0
  } # 
  outs <- numeric(length(rt))
  bounds <- c(0,rt)
  #browser()
  for (i in 1:length(rt)) {
    if (bounds[i]>=bounds[i+1]) {
      outs[i]=0
      next
    }
    tmp_obj <- do.call(integrate, args=c(f=n1PDF_core,lower=bounds[i],upper=bounds[i+1],subdivisions=1000, A=ret_arg(A, i), b=ret_arg(b, i), t0 = ret_arg(t0, i), st0 = list(st0[i]), sapply(dots, function(z, i) sapply(z, "[[", i = i, simplify=FALSE), i=i, simplify=FALSE), pdf = pdf, cdf = cdf, stop.on.error = FALSE, args.dist = list(args.dist)))
    if (tmp_obj$message != "OK") {
      warning(tmp_obj$message)
    }
    outs[i] <- tmp_obj$value
  }
  cumsum(outs)
}
