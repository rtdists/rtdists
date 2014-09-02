#' The Ratcliff Diffusion Model
#' 
#' Density, distribution function, and random generation for the Ratcliff diffusion model with eight parameters: \code{a} (threshold separation), \code{z} (relative starting point), \code{v} (drift rate), \code{t0} (non-decision time/response time constant), \code{d} (differences in speed of response execution), \code{sv} (inter-trial-variability of drift), \code{st0} (inter-trial-variability of non-decisional components), and \code{sz} (inter-trial-variability of relative starting point). Parameters can either be specified individually or as a vector via \code{parameters}.
#'
#' @param t a vector of RTs.
#' @param n desired number of observations.
#' @param boundary character vector. Which boundary should be tested. Possible values are \code{c("upper", "lower")}, possibly abbreviated and \code{"upper"} being the default.
#' 
#' @param a threshold separation. Amount of information that is considered for a decision. Large values indicate a conservative decisional style. Typical range: 0.5 < \code{a} < 2
#' @param z relative starting point. Indicator of an a priori bias in decision making. When the relative starting point \code{z} deviates from 0.5, the amount of information necessary for a decision differs between response alternatives. Typical range: 0.3 < \code{z} < 0.7
#' @param v drift rate. Average slope of the information accumulation process. The drift gives information about the speed and direction of the accumulation of information. Large (absolute) values of drift indicate a good performance. If received information supports the response linked to the upper threshold the sign will be positive and vice versa. Typical range: -5 < \code{v} < 5
#' @param t0 non-decision time or response time constant (in seconds). Average duration of all non-decisional processes (encoding and response execution). Typical range: 0.1 < \code{t0} < 0.5
#' @param d differences in speed of response execution (in seconds). Positive values indicate that response execution is faster for responses linked to the upper threshold than for responses linked to the lower threshold. Typical range: -0.1 < \code{d} < 0.1
#' @param sz inter-trial-variability of (relative) starting point. Range of a uniform distribution with mean \code{z} describing the distribution of actual starting points from specific trials. Minimal impact on the RT distributions. Can be fixed to 0 in most applications. Typical range: 0 < \code{sz} < 0.5
#' @param sv inter-trial-variability of drift rate. Standard deviation of a normal distribution with mean \code{v} describing the distribution of actual drift rates from specific trials. 	Minimal impact on the RT distributions. Can be fixed to 0 in most applications. Typical range: 0 < \code{sv} < 2
#' @param st0 inter-trial-variability of non-decisional components. Range of a uniform distribution with mean \code{t0} describing the distribution of actual \code{t0} values across trials. Accounts for response times below \code{t0}. Reduces skew of predicted RT distributions. Typical range: 0 < \code{st0} < 0.2
#' @param parameters a numeric vector of paramaters, can be used as an alternative to specifying all parameters via individual arguments. If named, order is irrelevant. If unnamed, following order is necessary: \code{c("a","z","v","t0","d","sz","sv","st0")}. 
#' 
#' @param precision \code{numerical} scalar value. Precision of calculation. Corresponds roughly to the number of decimals of the predicted CDFs that are calculated accuratly. Default is 3.
#' @param maxt maximum \code{rt} allowed, used to stop integration problems (\code{prd} only).
#'
#' @return \code{drd} gives the density, \code{prd} gives the distribution function, and \code{rrd} generates random response times and decisions (returning a \code{data.frame} with columns \code{rts} (numeric) and \code{response} (factor)).
#'
#' @details The Ratcliff diffusion model (Ratcliff, 1978) is a mathematical model for two-choice discrimination tasks. It is based on the assumption that information is accumulated continuously until one of two decision thresholds is hit. For more information, see Voss, Rothermund, and Voss (2004), Voss, Nagler, and Lerche (2013), or Wagenmakers (2009).
#' 
#' @references Ratcliff, R. (1978). A theory of memory retrieval. \emph{Psychological Review}, 85(2), 59-108.
#' 
#' Voss, A., Rothermund, K., & Voss, J. (2004). Interpreting the parameters of the diffusion model: An empirical validation. \emph{Memory & Cognition}. Vol 32(7), 32, 1206-1220.
#' 
#' Voss, A., Nagler, M., & Lerche, V. (2013). Diffusion Models in Experimental Psychology: A Practical Introduction. \emph{Experimental Psychology}, 60(6), 385-402. doi:10.1027/1618-3169/a000218
#' 
#' Wagenmakers, E.-J. (2009). Methodological and empirical developments for the Ratcliff diffusion model of response times and accuracy. \emph{European Journal of Cognitive Psychology}, 21(5), 641-671.
#' 
#' 
#' @author Underlying C code by Jochen Voss and Andreas Voss. Porting and R wrapping by Matthew Gretton, Andrew Heathcote, and Henrik Singmann.
#'
#' @useDynLib rtdists, dfastdm_b, pfastdm_b, rfastdm
#'
#' @name diffusion
#' 
#' @example examples/examples.diffusion.R
#' 
NULL

get_rd_parameters <- function(a, z, v, t0, d, sz, sv, st0, pl) {
  miss_p_indiv <- c(missing(a), missing(v), missing(t0), missing(z), missing(d), missing(sz), missing(sv), missing(st0))
  miss_pl <- missing(pl)
  if(all(miss_p_indiv) & miss_pl) stop("No parameter values specified.")
  if(any(!miss_p_indiv) & !miss_pl) stop("Parameter values need to be specified either via individual arguments or via 'pl', but not both.")
  if((sum(miss_p_indiv) > 0) & miss_pl) stop("Not all parameter values specified via individual arguments.")
  if(miss_pl) {
    rd_parameters <- c(a,v,t0,d,sz,sv,st0,z)
    if(length(rd_parameters)!=8) stop("Each parameter needs to be of length 1.")
    if(!is.numeric(rd_parameters)) stop("Parameters need to be numeric.")
    if (any(is.na(rd_parameters)) || !all(is.finite(rd_parameters))) stop("Parameters need to be numeric and finite.")
    return(rd_parameters)
  } else {
    if (!is.numeric(pl) | !is.vector(pl)) stop("pl needs to be a numeric vector.")
    if(length(pl)!=8) stop("Length of pl is incorrect (must be 8).")
    if (any(is.na(pl)) || !all(is.finite(pl))) stop("Parameters need to be numeric and finite.")
    if (!is.null(names(pl))) {
      if(!all(sapply(c("a","v","t0","d","sz","sv","st0","z"), function(x) x %in% names(pl)))) stop("names(pl) not complete. Need to be: c(\"a\",\"v\",\"t0\",\"z\",\"d\",\"sz\",\"sv\",\"st0\")")
      pl <- pl[c("a","v","t0","d","sz","sv","st0","z")]
    } else {
      names(pl) <- c("a","z","v","t0","d","sz","sv","st0")
      pl <- pl[c("a","v","t0","d","sz","sv","st0","z")]
    }
    return(pl)
  }
}

#

#' @rdname diffusion
#' @export drd
drd <- function (t, boundary = c("upper", "lower"), a, z, v, t0, d, sz, sv, st0, parameters, precision = 3)
{
  # Check for illegal parameter values
  pl <- get_rd_parameters(a = a, z = z, v = v, t0 = t0, d = d, sz = sz, sv = sv, st0 = st0, pl = parameters)
  
  boundary <- match.arg(boundary)
  if (boundary == "upper") i <- 2L
  if (boundary == "lower") i <- 1L
  
  # Call the C code
  densities <- vector(length=length(t))    
  output <- .C("dfastdm_b", 
               as.integer (length(t)),                 # 1  IN:  number of densities
               as.vector  (pl),                        # 2  IN:  parameters
               as.vector  (t),                         # 3  IN:  RTs
               as.double  (precision),                 # 4  IN:  precision
               as.integer (i),                         # 5  IN:  boundart 
               as.vector  (densities, mode="numeric")  # 6 OUT:  densities
  )
  
  unlist(output[6])
  
}

#' @rdname diffusion
#' @export prd
# set maximum t value to stop integration problems
prd <- function (t, boundary = c("upper", "lower"), a, z, v, t0, d, sz, sv, st0, parameters, precision = 3, maxt = 1e4) 
{
  # Check for illegal parameter values
  pl <- get_rd_parameters(a = a, z = z, v = v, t0 = t0, d = d, sz = sz, sv = sv, st0 = st0, pl = parameters)
  t[t>maxt] <- maxt
  if(!all(t == sort(t)))  stop("t needs to be sorted")
  
  boundary <- match.arg(boundary)
  if (boundary == "upper") i <- 2L
  if (boundary == "lower") i <- 1L
  
  # Call the C code
  pvalues <- vector(length=length(t))    
  output <- .C("pfastdm_b", 
               as.integer (length(t)),               # 1  IN:  number of densities
               as.vector  (pl),                      # 2  IN:  parameters
               as.vector  (t),                       # 3  IN:  RTs
               as.double  (precision),               # 4  IN:  number of densities
               as.integer (i),                       # 5  IN:  boundary 
               as.vector  (pvalues, mode="numeric")  # 6 OUT:  pvalues
  )
  
  unlist(output[6])
  
}

#' @rdname diffusion
#' @export rrd
# Returns a matrix of 2 x n (RTs x boundaries)
rrd <- function (n, a, z, v, t0, d, sz, sv, st0, parameters, precision = 3)
{
  randRTs    <- vector(length=n)
  randBounds <- vector(length=n)
  pl <- get_rd_parameters(a = a, z = z, v = v, t0 = t0, d = d, sz = sz, sv = sv, st0 = st0, pl = parameters)

  output <- .C("rfastdm", 
               as.integer (n),                          # 1  IN:  number of densities
               as.vector  (pl),                         # 2  IN:  parameters
               as.double  (precision),                  # 3  IN:  precision
               as.vector  (randRTs, mode="numeric"),    # 4 OUT:  RTs 
               as.vector  (randBounds, mode="numeric")  # 5 OUT:  bounds 
  )
  
  randRTs <- unlist(output[4]) 
  randBounds <- unlist(output[5])
  response <- factor(randBounds, levels = 0:1, labels = c("lower", "upper"))
  data.frame(rt = randRTs, response)
}

