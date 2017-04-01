#' The Ratcliff Diffusion Model
#' 
#' Density, distribution function, quantile function, and random generation for the Ratcliff diffusion model with following parameters: \code{a} (threshold separation), \code{z} (starting point), \code{v} (drift rate), \code{t0} (non-decision time/response time constant), \code{d} (differences in speed of response execution), \code{sv} (inter-trial-variability of drift), \code{st0} (inter-trial-variability of non-decisional components), \code{sz} (inter-trial-variability of relative starting point), and \code{s} (diffusion constant). \strong{Note that the parameterization or defaults of non-decision time variability \code{st0} and diffusion constant \code{s} differ from what is often found in the literature and that the parameterization of \code{z} and \code{sz} have changed compared to previous versions (now absolute and not relative).}
#'
#' @param rt a vector of RTs. Or for convenience also a \code{data.frame} with columns \code{rt} and \code{response} (such as returned from \code{rdiffusion} or \code{\link{rLBA}}). See examples.
#' @param n is a desired number of observations.
#' @param response character vector. Which response boundary should be tested? Possible values are \code{c("upper", "lower")}, possibly abbreviated and \code{"upper"} being the default. Alternatively, a numeric vector with values 1=lower and 2=upper. For convenience, \code{response} is converted via \code{as.numeric} also allowing factors (see examples). Ignored if the first argument is a \code{data.frame}.
#' @param p vector of probabilities. Or for convenience also a \code{data.frame} with columns \code{p} and \code{response}. See examples.
#' 
#' @param a threshold separation. Amount of information that is considered for a decision. Large values indicate a conservative decisional style. Typical range: 0.5 < \code{a} < 2
#' @param v drift rate. Average slope of the information accumulation process. The drift gives information about the speed and direction of the accumulation of information. Large (absolute) values of drift indicate a good performance. If received information supports the response linked to the upper threshold the sign will be positive and vice versa. Typical range: -5 < \code{v} < 5
#' @param t0 non-decision time or response time constant (in seconds). Lower bound for the duration of all non-decisional processes (encoding and response execution). Typical range: 0.1 < \code{t0} < 0.5
#' @param z starting point. Indicator of an a priori bias in decision making. When the relative starting point \code{z} deviates from \code{0.5*a}, the amount of information necessary for a decision differs between response alternatives. Default is \code{0.5*a} (i.e., no bias).
#' @param d differences in speed of response execution (in seconds). Positive values indicate that response execution is faster for responses linked to the upper threshold than for responses linked to the lower threshold. Typical range: -0.1 < \code{d} < 0.1. Default is 0.
#' @param sz inter-trial-variability of starting point. Range of a uniform distribution with mean \code{z} describing the distribution of actual starting points from specific trials. Values different from 0 can predict fast errors (but can slow computation considerably). Typical range: 0 < \code{sz} < 0.5. Default is 0.
#' @param sv inter-trial-variability of drift rate. Standard deviation of a normal distribution with mean \code{v} describing the distribution of actual drift rates from specific trials. Values different from 0 can predict slow errors. Typical range: 0 < \code{sv} < 2. Default is 0.
#' @param st0 inter-trial-variability of non-decisional components. Range of a uniform distribution with mean \code{t0 + st0/2} describing the distribution of actual \code{t0} values across trials. Accounts for response times below \code{t0}. Reduces skew of predicted RT distributions. Values different from 0 can slow computation considerably. Typical range: 0 < \code{st0} < 0.2. Default is 0.
#' @param s diffusion constant; standard deviation of the random noise of the diffusion process (i.e., within-trial variability), scales \code{a}, \code{v}, and \code{sv}. Needs to be fixed to a constant in most applications. Default is 1. Note that the default used by Ratcliff and in other applications is often 0.1. 
#' 
#' @param precision \code{numerical} scalar value. Precision of calculation. Corresponds roughly to the number of decimals of the predicted CDFs that are calculated accurately. Default is 3.
#' @param maxt maximum \code{rt} allowed, used to stop integration problems. \emph{Currently ignored as numerical integration works fine with Inf.}
#' @param interval a vector containing the end-points of the interval to be searched for the desired quantiles (i.e., RTs) in \code{qdiffusion}. Default is \code{c(0, 10)}.
#' @param scale_p logical. Should entered probabilities automatically be scaled by maximally predicted probability? Default is \code{FALSE}. Convenience argument for obtaining predicted quantiles. Can be slow as the maximally predicted probability is calculated individually for each \code{p}.
#' @param scale_max numerical scalar. Value at which maximally predicted RT should be calculated if \code{scale_p} is \code{TRUE}. 
#'
#' @return \code{ddiffusion} gives the density, \code{pdiffusion} gives the distribution function, \code{qdiffusion} gives the quantile function (i.e., predicted RTs), and \code{rdiffusion} generates random response times and decisions (returning a \code{data.frame} with columns \code{rt} (numeric) and \code{response} (factor)).
#' 
#' The length of the result is determined by \code{n} for \code{rdiffusion}, equal to the length of \code{rt} for \code{ddiffusion} and \code{pdiffusion}, and equal to the length of \code{p} for \code{qdiffusion}.
#' 
#' The distribution parameters (as well as \code{response}) are recycled to the length of the result. In other words, the functions are completely vectorized for all parameters and even the response boundary.
#'
#' @details The Ratcliff diffusion model (Ratcliff, 1978) is a mathematical model for two-choice discrimination tasks. It is based on the assumption that information is accumulated continuously until one of two decision thresholds is hit. For introductions see Ratcliff and McKoon (2008), Voss, Rothermund, and Voss (2004), Voss, Nagler, and Lerche (2013), or Wagenmakers (2009).
#' 
#' All functions are fully vectorized across all parameters as well as the response. This allows for trialwise parameters for each model parameter. 
#' 
#' For convenience, all functions (with the exception of \code{rdiffusion}) allow that the first argument is a \code{data.frame} containing the information of the first and second argument in two columns (i.e., \code{rt}/\code{p} and \code{response}). Other columns (as well as passing \code{response} separately argument) will be ignored. This allows, for example, to pass the \code{data.frame} generated by \code{rdiffusion} directly to \code{pdiffusion}. See examples.
#' 
#' \subsection{Distribution Function (CDF)}{
#' For the time being, the cumulative probability function (\code{pdiffusion}) uses numerical integration of \code{ddiffusion} via \code{\link{integral}} from package \pkg{pracma}. While this is somewhat slow, it reliably produces correct results (in most cases). Please report all issues with the CDF to the package maintainer.
#' }
#' 
#' \subsection{Quantile Function}{
#' Due to the bivariate nature of the diffusion model, the diffusion processes reaching each response boundary only return the defective CDF that does not reach 1. Only the sum of the CDF for both boundaries reaches 1. Therefore, \code{qdiffusion} can only return quantiles/RTs for any accumulator up to the maximal probability of that accumulator's CDF. This can be obtained by evaluating the CDF at \code{Inf}. 
#' 
#' As a convenience for the user, if \code{scale_p = TRUE} in the call to \code{qdiffusion} the desired probabilities are automatically scaled by the maximal probability for the corresponding response. Note that this can be slow as the maximal probability is calculated separately for each desired probability. See examples.
#' 
#' Also note that quantiles (i.e., predicted RTs) are obtained by numerically minimizing the absolute difference between desired probability and the value returned from \code{pdiffusion} using \code{\link{optimize}}. If the difference between the desired probability and probability corresponding to the returned quantile is above a certain threshold (currently 0.0001) no quantile is returned but \code{NA}. This can be either because the desired quantile is above the maximal probability for this accumulator or because the limits for the numerical integration are too small (default is \code{c(0, 10)}).
#' }
#' 
#' @note The parameterization of the non-decisional components, \code{t0} and \code{st0}, differs from the parameterization used by, for example, Andreas Voss or Roger Ratcliff. In the present case \code{t0} is the lower bound of the uniform distribution of length \code{st0}, but \emph{not} its midpoint. The parameterization employed here is in line with the parametrization for the \link{LBA} code (where \code{t0} is also the lower bound).
#' 
#' The default diffusion constant \code{s} is 1 and not 0.1 as in most applications of Roger Ratcliff and others.
#' 
#' We have changed the parameterization of the start point \code{z} which is now the absolute start point in line with most published literature (it was the relative start point in previous versions of \pkg{rtdists}). 
#' 
#' RTs need to be sorted (in increasing order) for \code{pdiffusion}.
#' 
#' @references Ratcliff, R. (1978). A theory of memory retrieval. \emph{Psychological Review}, 85(2), 59-108.
#' 
#' Ratcliff, R., & McKoon, G. (2008). The diffusion decision model: Theory and data for two-choice decision tasks. \emph{Neural Computation}, 20(4), 873-922.
#' 
#' Voss, A., Rothermund, K., & Voss, J. (2004). Interpreting the parameters of the diffusion model: An empirical validation. \emph{Memory & Cognition}. Vol 32(7), 32, 1206-1220.
#' 
#' Voss, A., Nagler, M., & Lerche, V. (2013). Diffusion Models in Experimental Psychology: A Practical Introduction. \emph{Experimental Psychology}, 60(6), 385-402. doi:10.1027/1618-3169/a000218
#' 
#' Wagenmakers, E.-J., van der Maas, H. L. J., & Grasman, R. P. P. P. (2007). An EZ-diffusion model for response time and accuracy. \emph{Psychonomic Bulletin & Review}, 14(1), 3-22.
#' 
#' Wagenmakers, E.-J. (2009). Methodological and empirical developments for the Ratcliff diffusion model of response times and accuracy. \emph{European Journal of Cognitive Psychology}, 21(5), 641-671.
#' 
#' 
#' @author Underlying C code by Jochen Voss and Andreas Voss. Porting and R wrapping by Matthew Gretton, Andrew Heathcote, Scott Brown, and Henrik Singmann. \code{qdiffusion} by Henrik Singmann.
#'
#' @useDynLib rtdists
#'
#' @name Diffusion
#' @importFrom utils head
#' @importFrom stats optimize uniroot
#' @importFrom pracma integral
#' @aliases diffusion
#' @importFrom Rcpp evalCpp
#' 
#' @example examples/examples.diffusion.R
#' 


# [MG 20150616]
# In line with LBA, adjust t0 to be the lower bound of the non-decision time distribution rather than the average 
# Called from prd, drd, rrd 
recalc_t0 <- function (t0, st0) { t0 <- t0 + st0/2 }


#' @rdname Diffusion
#' @export
ddiffusion <- function (rt, response = "upper", 
                 a, v, t0, z = 0.5*a, d = 0, sz = 0, sv = 0, st0 = 0, s = 1,
                 precision = 3, stop_on_error = TRUE)
{
  if(any(missing(a), missing(v), missing(t0))) stop("a, v, and/or t0 must be supplied")
  
  # for convenience accept data.frame as first argument.
  if (is.data.frame(rt)) {
    response <- rt$response
    rt <- rt$rt
  }
  
  nn <- length(rt)
  # Build parameter matrix  
  # Convert boundaries to numeric if necessary
  if (is.character(response)) {
    response <- match.arg(response, choices=c("upper", "lower"),several.ok = TRUE)
    numeric_bounds <- ifelse(response == "upper", 2L, 1L)
    }
  else {
    response <- as.numeric(response)
    if(any(!(response %in% 1:2))) stop("response needs to be either 'upper', 'lower', or as.numeric(response) %in% 1:2!")
    numeric_bounds <- as.integer(response)
  }
  numeric_bounds <- rep(numeric_bounds, length.out = nn)
  # all parameters brought to length of rt
  s <- rep(s, length.out = nn)
  a <- rep(a, length.out = nn)
  v <- rep(v, length.out = nn)
  t0 <- rep(t0, length.out = nn)
  z <- rep(z, length.out = nn)
  z <- z/a  # transform z from absolute to relative scale (which is currently required by the C code)
  d <- rep(d, length.out = nn)
  sz <- rep(sz, length.out = nn)
  sz <- sz/a # transform sz from absolute to relative scale (which is currently required by the C code)
  sv <- rep(sv, length.out = nn)
  st0 <- rep(st0, length.out = nn)
  t0 <- recalc_t0 (t0, st0) 
  
  # Build parameter matrix (and divide a, v, and sv, by s)
  params <- cbind (a/s, v/s, t0, d, sz, sv/s, st0, z, numeric_bounds)
  
  # Check for illegal parameter values
  if(ncol(params)<9) stop("Not enough parameters supplied: probable attempt to pass NULL values?")
  if(!is.numeric(params)) stop("Parameters need to be numeric.")
  if (any(is.na(params)) || !all(is.finite(params))) stop("Parameters need to be numeric and finite.")
  
  uniques <- unique(params)
  densities <- vector("numeric",length=length(rt))  
  for (i in seq_len(nrow(uniques))) {
    ok_rows <- apply(params, 1, identical, y = uniques[i,])

    densities[ok_rows] <- d_fastdm (rt[ok_rows], uniques[i,1:8], precision, uniques[i,9], stop_on_error)
  }
  abs(densities)
}


## @param stop.on.error logical. If true (the default) an error stops the \code{integration} of \code{pdiffusion}. If false some errors will give a result with a warning in the message component.
#' @rdname Diffusion
#' @export
pdiffusion <- function (rt, response = "upper",
                 a, v, t0, z = 0.5*a, d = 0, sz = 0, sv = 0, st0 = 0, s = 1,
                 precision = 3, maxt = 20, stop_on_error = TRUE, use_precise = TRUE)
{
  if(any(missing(a), missing(v), missing(t0))) stop("a, v, and/or t0 must be supplied")
  # for convenience accept data.frame as first argument.
  if (is.data.frame(rt)) {
    response <- rt$response
    rt <- rt$rt
  }
  
  rt[rt>maxt] <- maxt
  if(!all(rt == sort(rt)))  stop("rt needs to be sorted")
  
  # Convert boundaries to numeric
  nn <- length(rt)
  # Build parameter matrix  
  # Convert boundaries to numeric
  if (is.character(response)) {
    response <- match.arg(response, choices=c("upper", "lower"),several.ok = TRUE)
    numeric_bounds <- ifelse(response == "upper", 2L, 1L)
  }
  else {
    response <- as.numeric(response)
    if(any(!(response %in% 1:2))) stop("response needs to be either 'upper', 'lower', or as.numeric(response) %in% 1:2!")
    numeric_bounds <- as.integer(response)
  }
  numeric_bounds <- rep(numeric_bounds, length.out = nn)
  # all parameters brought to length of rt
  s <- rep(s, length.out = nn)
  a <- rep(a, length.out = nn)
  v <- rep(v, length.out = nn)
  t0 <- rep(t0, length.out = nn)
  z <- rep(z, length.out = nn)
  z <- z/a  # transform z from absolute to relative scale (which is currently required by the C code)
  d <- rep(d, length.out = nn)
  sz <- rep(sz, length.out = nn)
  sz <- sz/a # transform sz from absolute to relative scale (which is currently required by the C code)
  sv <- rep(sv, length.out = nn)
  st0 <- rep(st0, length.out = nn)
  t0 <- recalc_t0 (t0, st0) 
  
  # Build parameter matrix (and divide a, v, and sv, by s)
  params <- cbind (a/s, v/s, t0, d, sz, sv/s, st0, z, numeric_bounds)
  
  # Check for illegal parameter values
  if(ncol(params)<9) stop("Not enough parameters supplied: probable attempt to pass NULL values?")
  if(!is.numeric(params)) stop("Parameters need to be numeric.")
  if (any(is.na(params)) || !all(is.finite(params))) stop("Parameters need to be numeric and finite.")
  
  pvalues <- vector("numeric", length=length(rt))
  uniques <- unique(params)
  
  if (use_precise) {
    for (i in seq_len(nrow(uniques))) {
      ok_rows <- apply(params, 1, identical, y = uniques[i,])
      pvalues[ok_rows] <- p_precise_fastdm (rt[ok_rows], uniques[i,1:8], precision, uniques[i,9], stop_on_error)
    }
  } else {
    for (i in seq_len(nrow(uniques))) {
      ok_rows <- apply(params, 1, identical, y = uniques[i,])
      pvalues[ok_rows] <- p_fastdm (rt[ok_rows], uniques[i,1:8], precision, uniques[i,9], stop_on_error)
    }
  }
  pvalues
}


inv_cdf_diffusion <- function(x, response, a, v, t0, z, d, sz, sv, st0, s, precision, maxt, value, abs = TRUE, stop_on_error = TRUE) {
  if (abs) abs(value - pdiffusion(rt=x, response=response, a=a, v=v, t0=t0, z=z, d=d, sz=sz, sv=sv, s=s, st0=st0, precision=precision, maxt=maxt, stop_on_error))
  else (value - pdiffusion(rt=x, response=response, a=a, v=v, t0=t0, z=z, d=d, sz=sz, sv=sv, st0=st0, s=s, precision=precision, maxt=maxt, stop_on_error))
}

#' @rdname Diffusion
#' @export
qdiffusion <- function (p, response = "upper", 
                 a, v, t0, z = 0.5*a, d = 0, sz = 0, sv = 0, st0 = 0, s = 1,
                 precision = 3, maxt = 20, interval = c(0, 10),
                 scale_p = FALSE, scale_max = Inf, stop_on_error = TRUE)
{
  if(any(missing(a), missing(v), missing(t0))) stop("a, v, and t0 must be supplied")

  # for convenience accept data.frame as first argument.
  if (is.data.frame(p)) {
    response <- p$response
    p <- p$p
  }
  
  nn <- length(p)
  response <- rep(unname(response), length.out = nn)
  s <- rep(s, length.out = nn) # pass s to other functions for correct handling! No division here.
  a <- rep(unname(a), length.out = nn)
  v <- rep(unname(v), length.out = nn)
  t0 <- rep(unname(t0), length.out = nn)
  z <- rep(unname(z), length.out = nn)
  d <- rep(unname(d), length.out = nn)
  sz <- rep(unname(sz), length.out = nn)
  sv <- rep(unname(sv), length.out = nn)
  st0 <- rep(unname(st0), length.out = nn)
  p <- unname(p)
  
  out <- vector("numeric", nn)
  for (i in seq_len(nn)) {
    if (scale_p) max_p <- pdiffusion(scale_max, response=response[i], a=a[i], v=v[i], t0=t0[i], z=z[i], d=d[i], sz=sz[i], sv=sv[i], st0=st0[i], s=s[i], precision=precision, maxt=maxt, stop_on_error=stop_on_error)
    else max_p <- 1
    tmp <- do.call(optimize, args = c(f=inv_cdf_diffusion, interval = list(interval), response=response[i], a=a[i], v=v[i], t0=t0[i], z=z[i], d=d[i], sz=sz[i], sv=sv[i], st0=st0[i], s=s[i], precision=precision, maxt=maxt, stop_on_error=stop_on_error, value =p[i]*max_p, tol = .Machine$double.eps^0.5))
    if (tmp$objective > 0.0001) {
      tmp <- do.call(optimize, args = c(f=inv_cdf_diffusion, interval = list(c(min(interval),max(interval)/2)), response=response[i], a=a[i], v=v[i], t0=t0[i], z=z[i], d=d[i], sz=sz[i], sv=sv[i], st0=st0[i], s=s[i], precision=precision, maxt=maxt, stop_on_error=stop_on_error, value =p[i]*max_p, tol = .Machine$double.eps^0.5))
    }
    if (tmp$objective > 0.0001) {
      try({
        uni_tmp <- do.call(uniroot, args = c(f=inv_cdf_diffusion, interval = list(interval), response=response[i], a=a[i], v=v[i], t0=t0[i], z=z[i], d=d[i], sz=sz[i], sv=sv[i], st0=st0[i], precision=precision, s=s[i], maxt=maxt, stop_on_error=stop_on_error, value =p[i]*max_p, tol = .Machine$double.eps^0.5, abs = FALSE))
      tmp$objective <- uni_tmp$f.root
      tmp$minimum <- uni_tmp$root
      }, silent = TRUE)
    }
    if (tmp$objective > 0.0001) {
      tmp[["minimum"]] <- NA
      warning("Cannot obtain RT that is less than 0.0001 away from desired p = ", p[i], ".\nIncrease/decrease interval or obtain for different response.", call. = FALSE)
    }
    out[i] <- tmp[["minimum"]]
  }
  return(out)
}

## When given vectorised parameters, n is the number of replicates for each parameter set
#' @rdname Diffusion
#' @export
rdiffusion <- function (n, 
                 a, v, t0, z = 0.5*a, d = 0, sz = 0, sv = 0, st0 = 0, s = 1,
                 precision = 3, stop_on_error = TRUE)
{
  if(any(missing(a), missing(v), missing(t0))) stop("a, v, and/or t0 must be supplied")
  
  s <- rep(s, length.out = n)
  a <- rep(a, length.out = n)
  v <- rep(v, length.out = n)
  t0 <- rep(t0, length.out = n)
  z <- rep(z, length.out = n)
  z <- z/a  # transform z from absolute to relative scale (which is currently required by the C code)
  d <- rep(d, length.out = n)
  sz <- rep(sz, length.out = n)
  sz <- sz/a # transform sz from absolute to relative scale (which is currently required by the C code)
  sv <- rep(sv, length.out = n)
  st0 <- rep(st0, length.out = n)
  t0 <- recalc_t0 (t0, st0) 
  
  # Build parameter matrix (and divide a, v, and sv, by s)
  params <- cbind (a/s, v/s, t0, d, sz, sv/s, st0, z)
  
  # Check for illegal parameter values
  if(ncol(params)<8) stop("Not enough parameters supplied: probable attempt to pass NULL values?")
  if(!is.numeric(params)) stop("Parameters need to be numeric.")
  if (any(is.na(params)) || !all(is.finite(params))) stop("Parameters need to be numeric and finite.")
  
  randRTs    <- vector("numeric",length=n)
  randBounds <- vector("numeric",length=n)

  uniques <- unique(params)
  for (i in seq_len(nrow(uniques))) {
      ok_rows <- apply(params, 1, identical, y = uniques[i,])
    
      # Calculate n for this row
      current_n <- sum(ok_rows)
      
      out <- r_fastdm (current_n, uniques[i,1:8], precision, stop_on_error=stop_on_error)
      
      randRTs[ok_rows]    <- out$rt       
      randBounds[ok_rows] <- out$boundary 
  }
  response <- factor(randBounds, levels = 0:1, labels = c("lower", "upper"))
  data.frame(rt = randRTs, response)
}

