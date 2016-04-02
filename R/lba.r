#' The Linear Ballistic Accumulator (LBA)
#' 
#' Density, distribution function, quantile function, and random generation for the LBA model with the following parameters: \code{A} (upper value of starting point), \code{b} (response threshold), \code{t0} (non-decision time), and driftrate (\code{v}). All functions are available with different distributions underlying the drift rate: Normal (\code{norm}), Gamma (\code{gamma}), Frechet (\code{frechet}), and log normal (\code{lnorm}). The functions return their values conditional on the accumulator given in the response argument winning.  
#' 
#' @param rt vector of RTs.
#' @param response integer vector of winning accumulators/responses corresponding to the vector of RTs/p (i.e., used for specifying the response for a given RT/probability). Will be recycled if necessary. Cannot contain values larger than the number of accumulators. First response/accumulator must receive value 1, second 2, and so forth.
#' @param p vector of probabilities.
#' @param n desired number of observations (scalar integer).
#' @param A start point interval or evidence in accumulator before beginning of decision process. Start point varies from trial to trial in the interval [0, \code{A}] (uniform distribution). Average amount of evidence before evidence accumulation across trials is \code{A}/2.
#' @param b response threshold. (\code{b} - \code{A}/2) is a measure of "response caution". 
#' @param t0 non-decision time or response time constant (in seconds). Lower bound for the duration of all non-decisional processes (encoding and response execution).
#' @param st0 variability of non-decision time, such that \code{t0} is uniformly distributed between \code{t0} and \code{t0} + \code{st0}. Default is 0. Can be trialwise, and will be recycled to length of \code{rt}.
#' @param ... two \emph{named} drift rate parameters depending on \code{distribution} (e.g., \code{mean_v} and \code{sd_v} for \code{distribution=="norm"}). The parameters can either be given as a numeric vector or a list. If a numeric vector is passed each element of the vector corresponds to one accumulator. If a list is passed each list element corresponds to one accumulator allowing again trialwise driftrates. The shorter parameter will be recycled as necessary (and also the elements of the list to match the length of \code{rt}). See details.
#' @param distribution character specifying the distribution of the drift rate. Possible values are \code{c("norm", "gamma", "frechet", "lnorm")}, default is \code{"norm"}.
#' @param args.dist list of optional further arguments to the distribution functions (i.e., \code{posdrift} or \code{robust} for \code{distribution=="norm"}, see \code{\link{single-LBA}}).
#' @param silent logical. Should the number of accumulators used be suppressed? Default is \code{FALSE} which prints the number of accumulators.
#' @param interval a vector containing the end-points of the interval to be searched for the desired quantiles (i.e., RTs) in \code{qLBA}. Default is \code{c(0, 10)}.
#' 
#' @details 
#' \subsection{Parameters}{
#' The following arguments are allowed as \code{...} drift rate parameters:
#' \itemize{
#' \item \code{mean_v,sd_v} mean and standard deviation of normal distribution for drift rate (\code{norm}). See \code{\link{Normal}}
#' \item \code{shape_v,rate_v,scale_v} shape, rate, and scale of gamma (\code{gamma}) and scale and shape of Frechet (\code{frechet}) distributions for drift rate. See \code{\link{GammaDist}} or \code{\link[evd]{frechet}}. For Gamma, scale = 1/shape and shape = 1/scale.
#' \item \code{meanlog_v,sdlog_v} mean and standard deviation of lognormal distribution on the log scale for drift rate (\code{lnorm}). See \code{\link{Lognormal}}.
#' }
#' 
#' As described above, the accumulator parameters can either be given as a numeric vector or a list. If a numeric vector is passed each element of the vector corresponds to one accumulator. If a list is passed each list element corresponds to one accumulator allowing trialwise driftrates. The shorter parameter will be recycled as necessary (and also the elements of the list to match the length of \code{rt}).
#' 
#' The other LBA parameters (i.e., \code{A}, \code{b}, and \code{t0}, with the exception of \code{st0}) can either be a single numeric vector (which will be recycled to reach \code{length(rt)} or \code{length(n)} for trialwise parameters) \emph{or} a \code{list} of such vectors in which each list element corresponds to the parameters for this accumulator (i.e., the list needs to be of the same length as there are accumulators). Each list will also be recycled to reach \code{length(rt)} for trialwise parameters per accumulator.
#' 
#' To make the difference between both paragraphs clear: Whereas for the accumulators both a single vector or a list corresponds to different accumulators, only the latter is true for the other parameters. For those (i.e., \code{A}, \code{b}, and \code{t0}) a single vector always corresponds to trialwise values and a list must be used for accumulator wise values.
#' 
#' \code{st0} can only vary trialwise (via a vector). And it should be noted that \code{st0} not equal to zero will considerably slow done everything.
#' }
#' 
#' \subsection{Quantile Function}{
#' Due to the bivariate nature of the LBA, single accumulators only return defective CDFs that do not reach 1. Only the sum of all accumulators reaches 1. Therefore, \code{qLBA} can only return quantiles/RTs for any accumulator up to the maximal probability of that accumulator's CDF. This can be obtained by evaluating the CDF at \code{Inf} (see examples). 
#' 
#' Also note that quantiles (i.e., predicted RTs) are obtained by numerically minimizing the absolute difference between desired probabiliy and the value returned from \code{pLBA} using \code{\link{optimize}}. If the difference between the desired probability and probability corresponding to the returned quantile is above a certain threshold (currently 0.0001) no quantile is returned but \code{NA}. This can be either because the desired quantile is above the maximal probability for this accumulator or because the limits for the numerical integration are too small (default is \code{c(0, 10)}).
#' }
#' 
#' \subsection{RNG}{
#' For random number generation at least one of the distribution parameters (i.e., \code{mean_v}, \code{sd_v}, \code{shape_v}, \code{scale_v}, \code{rate_v}, \code{meanlog_v}, and \code{sdlog_v}) should be of length > 1 to receive RTs from multiple responses. Shorter vectors are recycled as necessary.\cr
#' Note that for random number generation from a normal distribution for the driftrate the number of returned samples may be less than the number of requested samples if \code{posdrifts==FALSE}.
#' }
#' 
#' 
#' @return \code{dLBA} returns the density (PDF), \code{pLBA} returns the distribution function (CDF), \code{qLBA} returns the quantile/RT, \code{rLBA} return random response times and responses (in a \code{data.frame}).
#' 
#' @note These are the top-level functions intended for end-users. To obtain the density and cumulative density the race functions are called for each response time with the corresponding winning accumulator as first accumulator (see \code{\link{LBA-race}}). 
#' 
#' @references 
#' 
#' Brown, S. D., & Heathcote, A. (2008). The simplest complete model of choice response time: Linear ballistic accumulation. \emph{Cognitive Psychology}, 57(3), 153-178. doi:10.1016/j.cogpsych.2007.12.002
#' 
#' Donkin, C., Averell, L., Brown, S., & Heathcote, A. (2009). Getting more from accuracy and response time data: Methods for fitting the linear ballistic accumulator. \emph{Behavior Research Methods}, 41(4), 1095-1110. doi:10.3758/BRM.41.4.1095
#' 
#' Heathcote, A., & Love, J. (2012). Linear deterministic accumulator models of simple choice. \emph{Frontiers in Psychology}, 3, 292. doi:10.3389/fpsyg.2012.00292
#' 
#' @name LBA
#' @importFrom stats optimize
#' 
#' @example examples/examples.lba.R
#' 
NULL

# this functions takes the argument as entered by the user and always returns a list of length n_v in which each element is of length nn. 
# This ensures that for each parameter we have a list of length equal to the number of accumulators with each element of length equal to the number of trials. This consistency will be exploited when passing to the n1functions. 
check_i_arguments <- function(arg, nn, n_v, dots = FALSE) {
  mc <- match.call()
  varname <- sub("dots$", "", deparse(mc[["arg"]]), fixed = TRUE)
  if (!is.list(arg)) {
    if ((!is.vector(arg, "numeric")) || (length(arg) < 1)) stop(paste(varname, "needs to be a numeric vector of length >= 1!"))
    if (dots) {
      arg <- as.list(arg)
      arg <- lapply(arg, rep, length.out=nn)
    } else arg <- lapply(seq_len(n_v), function(x) rep(arg, length.out=nn))
  } else {
    if (!dots && (length(arg) != n_v)) stop(paste("if", varname, "is a list, its length needs to correspond to the number of accumulators."))
    for (i in seq_along(arg)) {
      if ((!is.vector(arg[[i]], "numeric")) || (length(arg[[i]]) < 1)) stop(paste0(varname, "[[", i, "]] needs to be a numeric vector of length >= 1!"))
      arg[[i]] <- rep(arg[[i]], length.out=nn)
    }
  }
  #if (length(arg) != n_v) stop(paste("size of", varname, "does not correspond to number of accumulators."))
  return(arg)
}

#' @rdname LBA
#' @export 
dLBA <-  function(rt, response, A, b, t0, ..., st0=0, distribution = c("norm", "gamma", "frechet", "lnorm"), args.dist = list(), silent = FALSE) {
  dots <- list(...)
  if (is.null(names(dots))) stop("... arguments need to be named.")
  
  nn <- length(rt)
  n_v <- max(vapply(dots, length, 0))  # Number of responses
  if(!silent) message(paste("Results based on", n_v, "accumulators/drift rates."))
  if (!is.numeric(response) || max(response) > n_v) stop("response needs to be a numeric vector of integers up to number of accumulators.")
  if (any(response < 1)) stop("the first response/accumulator must have value 1.")
  if (n_v < 2) stop("There need to be at least two accumulators/drift rates.")
  distribution <- match.arg(distribution)
  response <- rep(response, length.out = nn)
  A <- check_i_arguments(A, nn=nn, n_v=n_v)
  b <- check_i_arguments(b, nn=nn, n_v=n_v)
  t0 <- check_i_arguments(t0, nn=nn, n_v=n_v)
  switch(distribution, 
         norm = {
           if (any(!(c("mean_v","sd_v") %in% names(dots)))) stop("mean_v and sd_v need to be passed for distribution = \"norm\"")
           dots$mean_v <- check_i_arguments(dots$mean_v, nn=nn, n_v=n_v, dots = TRUE)
           dots$sd_v <- check_i_arguments(dots$sd_v, nn=nn, n_v=n_v, dots = TRUE)
           dots <- dots[c("mean_v","sd_v")]
         },
         gamma = {
           if (!("shape_v" %in% names(dots))) stop("shape_v needs to be passed for distribution = \"gamma\"")
           if ((!("rate_v" %in% names(dots))) & (!("scale_v" %in% names(dots)))) stop("rate_v or scale_v need to be passed for distribution = \"gamma\"")
           dots$shape_v <- check_i_arguments(dots$shape_v, nn=nn, n_v=n_v, dots = TRUE)
           if ("scale_v" %in% names(dots)) {
             dots$scale_v <- check_i_arguments(dots$scale_v, nn=nn, n_v=n_v, dots = TRUE)
             if (is.list(dots$scale_v)) {
               dots$rate_v <- lapply(dots$scale_v, function(x) 1/x)
             } else dots$rate_v <- 1/dots$scale_v
           } else dots$rate_v <- check_i_arguments(dots$rate_v, nn=nn, n_v=n_v, dots = TRUE)
           dots <- dots[c("shape_v","rate_v")]
         },
         frechet = {
           if (any(!(c("shape_v","scale_v") %in% names(dots)))) stop("shape_v and scale_v need to be passed for distribution = \"frechet\"")
           dots$shape_v <- check_i_arguments(dots$shape_v, nn=nn, n_v=n_v, dots = TRUE)
           dots$scale_v <- check_i_arguments(dots$scale_v, nn=nn, n_v=n_v, dots = TRUE)
           dots <- dots[c("shape_v","scale_v")]
         },
         lnorm = {
           if (any(!(c("meanlog_v","sdlog_v") %in% names(dots)))) stop("meanlog_v and sdlog_v need to be passed for distribution = \"lnorm\"")
           dots$meanlog_v <- check_i_arguments(dots$meanlog_v, nn=nn, n_v=n_v, dots = TRUE)
           dots$sdlog_v <- check_i_arguments(dots$sdlog_v, nn=nn, n_v=n_v, dots = TRUE)
           dots <- dots[c("meanlog_v","sdlog_v")]
         }
  )
  #browser()
  for (i in seq_len(length(dots))) {
    if (length(dots[[i]]) < n_v) dots[[i]] <- rep(dots[[i]],length.out=n_v)
  }
  out <- vector("numeric", nn)
  for (i in unique(response)) {
    sel <- response == i
    out[sel] <- do.call(n1PDF, args = c(rt=list(rt[sel]), A = list(A[c(i, seq_len(n_v)[-i])]), b = list(b[c(i, seq_len(n_v)[-i])]), t0 = list(t0[c(i, seq_len(n_v)[-i])]), lapply(dots, function(x) x[c(i, seq_len(n_v)[-i])]), distribution=distribution, args.dist=args.dist, silent=TRUE, st0 = list(st0)))
  }
  return(out)
}


#' @rdname LBA
#' @export 
pLBA <-  function(rt, response, A, b, t0, ..., st0=0, distribution = c("norm", "gamma", "frechet", "lnorm"), args.dist = list(), silent = FALSE) {
  dots <- list(...)
  if (is.null(names(dots))) stop("... arguments need to be named.")
  
  nn <- length(rt)
  n_v <- max(vapply(dots, length, 0))  # Number of responses
  if(!silent) message(paste("Results based on", n_v, "accumulators/drift rates."))
  if (!is.numeric(response) || max(response) > n_v) stop("response needs to be a numeric vector of integers up to number of accumulators.")
  if (n_v < 2) stop("There need to be at least two accumulators/drift rates.")
  distribution <- match.arg(distribution)
  response <- rep(response, length.out = nn)
  A <- check_i_arguments(A, nn=nn, n_v=n_v)
  b <- check_i_arguments(b, nn=nn, n_v=n_v)
  t0 <- check_i_arguments(t0, nn=nn, n_v=n_v)
  switch(distribution, 
         norm = {
           if (any(!(c("mean_v","sd_v") %in% names(dots)))) stop("mean_v and sd_v need to be passed for distribution = \"norm\"")
           dots$mean_v <- check_i_arguments(dots$mean_v, nn=nn, n_v=n_v, dots = TRUE)
           dots$sd_v <- check_i_arguments(dots$sd_v, nn=nn, n_v=n_v, dots = TRUE)
           dots <- dots[c("mean_v","sd_v")]
         },
         gamma = {
           if (!("shape_v" %in% names(dots))) stop("shape_v needs to be passed for distribution = \"gamma\"")
           if ((!("rate_v" %in% names(dots))) & (!("scale_v" %in% names(dots)))) stop("rate_v or scale_v need to be passed for distribution = \"gamma\"")
           dots$shape_v <- check_i_arguments(dots$shape_v, nn=nn, n_v=n_v, dots = TRUE)
           if ("scale_v" %in% names(dots)) {
             dots$scale_v <- check_i_arguments(dots$scale_v, nn=nn, n_v=n_v, dots = TRUE)
             if (is.list(dots$scale_v)) {
               dots$rate_v <- lapply(dots$scale_v, function(x) 1/x)
             } else dots$rate_v <- 1/dots$scale_v
           } else dots$rate_v <- check_i_arguments(dots$rate_v, nn=nn, n_v=n_v, dots = TRUE)
           dots <- dots[c("shape_v","rate_v")]
         },
         frechet = {
           if (any(!(c("shape_v","scale_v") %in% names(dots)))) stop("shape_v and scale_v need to be passed for distribution = \"frechet\"")
           dots$shape_v <- check_i_arguments(dots$shape_v, nn=nn, n_v=n_v, dots = TRUE)
           dots$scale_v <- check_i_arguments(dots$scale_v, nn=nn, n_v=n_v, dots = TRUE)
           dots <- dots[c("shape_v","scale_v")]
         },
         lnorm = {
           if (any(!(c("meanlog_v","sdlog_v") %in% names(dots)))) stop("meanlog_v and sdlog_v need to be passed for distribution = \"lnorm\"")
           dots$meanlog_v <- check_i_arguments(dots$meanlog_v, nn=nn, n_v=n_v, dots = TRUE)
           dots$sdlog_v <- check_i_arguments(dots$sdlog_v, nn=nn, n_v=n_v, dots = TRUE)
           dots <- dots[c("meanlog_v","sdlog_v")]
         }
  )
  #browser()
  for (i in seq_len(length(dots))) {
    if (length(dots[[i]]) < n_v) dots[[i]] <- rep(dots[[i]],length.out=n_v)
  }
  out <- vector("numeric", nn)
  for (i in unique(response)) {
    sel <- response == i
    if(!all(rt[sel] == sort(rt[sel])))  stop("rt needs to be sorted (per response)")
    out[sel] <- do.call(n1CDF, args = c(rt=list(rt[sel]), A = list(A[c(i, seq_len(n_v)[-i])]), b = list(b[c(i, seq_len(n_v)[-i])]), t0 = list(t0[c(i, seq_len(n_v)[-i])]), lapply(dots, function(x) x[c(i, seq_len(n_v)[-i])]), distribution=distribution, args.dist=args.dist, silent=TRUE, st0 = list(st0)))
  }
  return(out)
}

# rt, response, A, b, t0, ..., st0=0, distribution = c("norm", "gamma", "frechet", "lnorm"), args.dist = list(), silent = FALSE
inv_cdf_lba <- function(x, response, A, b, t0, ..., st0, distribution, args.dist, value) {
  abs(value - pLBA(rt=x, response=response, A=A, b = b, t0 = t0, ..., st0=st0, distribution=distribution, args.dist=args.dist, silent=TRUE))
}


#' @rdname LBA
#' @export 
qLBA <-  function(p, response, A, b, t0, ..., st0=0, distribution = c("norm", "gamma", "frechet", "lnorm"), args.dist = list(), silent = FALSE, interval = c(0, 10)) {
  dots <- list(...)
  if (is.null(names(dots))) stop("... arguments need to be named.")
  
  nn <- length(p)
  n_v <- max(vapply(dots, length, 0))  # Number of responses
  if(!silent) message(paste("Results based on", n_v, "accumulators/drift rates."))
  if (!is.numeric(response) || max(response) > n_v) stop("response needs to be a numeric vector of integers up to number of accumulators.")
  if (any(response < 1)) stop("the first response/accumulator must have value 1.")
  if (n_v < 2) stop("There need to be at least two accumulators/drift rates.")
  distribution <- match.arg(distribution)
  response <- rep(response, length.out = nn)
  A <- check_i_arguments(A, nn=nn, n_v=n_v)
  b <- check_i_arguments(b, nn=nn, n_v=n_v)
  t0 <- check_i_arguments(t0, nn=nn, n_v=n_v)
  st0 <- rep(st0, length.out = nn)
  out <- vector("numeric", nn)
  for (i in seq_len(nn)) {
    tmp <- do.call(optimize, args = c(f=inv_cdf_lba, interval = list(interval), response = list(response[i]), A=ret_arg(A, i), b=ret_arg(b, i), t0=ret_arg(t0, i), sapply(dots, function(z, i) sapply(z, ret_arg2, which = i, simplify=FALSE), i=i, simplify=FALSE), args.dist = list(args.dist), distribution=distribution, st0 = list(st0[i]), value =p[i], tol = .Machine$double.eps^0.5))
    #browser()
    if (tmp$objective > 0.0001) {
      warning("Cannot obtain RT that is less than 0.0001 away from desired p = ", p[i], ".\nIncrease interval or obtain for different response.", call. = FALSE)
      out[i] <- NA
    } else out[i] <- tmp[[1]]
  }
  return(out)
}



#' @rdname LBA
#' @export
rLBA <- function(n,A,b,t0, ..., st0=0, distribution = c("norm", "gamma", "frechet", "lnorm"), args.dist = list(), silent = FALSE) {
  dots <- list(...)
  if (is.null(names(dots))) stop("... arguments need to be named.")
  if (any(names(dots) == "")) stop("all ... arguments need to be named.")
  n_v <- max(vapply(dots, length, 0))  # Number of responses
  if(!silent) message(paste("Results based on", n_v, "accumulators/drift rates."))
  #if (n_v < 2) stop("There need to be at least two accumulators/drift rates.")
  nn <- n
  distribution <- match.arg(distribution)
  A <- check_i_arguments(A, nn=nn, n_v=n_v)
  b <- check_i_arguments(b, nn=nn, n_v=n_v)
  t0 <- check_i_arguments(t0, nn=nn, n_v=n_v)
  st0 <- rep(unname(st0), length.out = nn)
  switch(distribution, 
         norm = {
           rng <- rlba_norm
           if (any(!(c("mean_v","sd_v") %in% names(dots)))) stop("mean_v and sd_v need to be passed for distribution = \"norm\"")
           dots$mean_v <- check_i_arguments(dots$mean_v, nn=nn, n_v=n_v, dots = TRUE)
           dots$sd_v <- check_i_arguments(dots$sd_v, nn=nn, n_v=n_v, dots = TRUE)
           dots <- dots[c("mean_v","sd_v")]
         },
         gamma = {
           rng <- rlba_gamma
           if (!("shape_v" %in% names(dots))) stop("shape_v needs to be passed for distribution = \"gamma\"")
           if ((!("rate_v" %in% names(dots))) & (!("scale_v" %in% names(dots)))) stop("rate_v or scale_v need to be passed for distribution = \"gamma\"")
           dots$shape_v <- check_i_arguments(dots$shape_v, nn=nn, n_v=n_v, dots = TRUE)
           if ("scale_v" %in% names(dots)) {
             dots$scale_v <- check_i_arguments(dots$scale_v, nn=nn, n_v=n_v, dots = TRUE)
             if (is.list(dots$scale_v)) {
               dots$rate_v <- lapply(dots$scale_v, function(x) 1/x)
             } else dots$rate_v <- 1/dots$scale_v
           } else dots$rate_v <- check_i_arguments(dots$rate_v, nn=nn, n_v=n_v, dots = TRUE)
           dots <- dots[c("shape_v","rate_v")]
         },
         frechet = {
           rng <- rlba_frechet
           if (any(!(c("shape_v","scale_v") %in% names(dots)))) stop("shape_v and scale_v need to be passed for distribution = \"frechet\"")
           dots$shape_v <- check_i_arguments(dots$shape_v, nn=nn, n_v=n_v, dots = TRUE)
           dots$scale_v <- check_i_arguments(dots$scale_v, nn=nn, n_v=n_v, dots = TRUE)
           dots <- dots[c("shape_v","scale_v")]
         },
         lnorm = {
           pdf <- rlba_lnorm
           if (any(!(c("meanlog_v","sdlog_v") %in% names(dots)))) stop("meanlog_v and sdlog_v need to be passed for distribution = \"lnorm\"")
           dots$meanlog_v <- check_i_arguments(dots$meanlog_v, nn=nn, n_v=n_v, dots = TRUE)
           dots$sdlog_v <- check_i_arguments(dots$sdlog_v, nn=nn, n_v=n_v, dots = TRUE)
           dots <- dots[c("meanlog_v","sdlog_v")]
         }
  )
  for (i in seq_len(length(dots))) {
    if (length(dots[[i]]) < n_v) dots[[i]] <- rep(dots[[i]],length.out=n_v)
  }
  
  tmp_acc <- as.data.frame(dots, optional = TRUE)
  colnames(tmp_acc) <- sub("\\.c\\(.+", "", colnames(tmp_acc))
  unique_acc <- unique(tmp_acc)
  
  out <- data.frame(rt = rep(0, n), response = 0)  
  for (i in seq_len(nrow(unique_acc))) {
    ok_rows <- apply(tmp_acc, 1, identical, y = as.matrix(unique_acc)[i,])
    tmp_dots <- lapply(dots, function(x) sapply(x, "[[", i = which(ok_rows)[1]))
    out[ok_rows,] <- do.call(rng, args = c(n=list(sum(ok_rows)), A = list(sapply(A, "[", i = ok_rows)), b = list(sapply(b, "[", i = ok_rows)), t0 = list(sapply(t0, "[", i = ok_rows)), st0 = list(st0[ok_rows]), tmp_dots, args.dist=args.dist))
    #print(str(c(n=list(sum(ok_rows)), A = list(sapply(A, "[", i = ok_rows)), b = list(sapply(b, "[", i = ok_rows)), t0 = list(sapply(t0, "[", i = ok_rows)), st0 = list(st0[ok_rows]), tmp_dots, args.dist=args.dist)))
  }
  
  out
}