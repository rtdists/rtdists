#' The Ratcliff Diffusion Model
#' 
#' Density, distribution function, and random generation for the Ratcliff diffusion model with eight parameters: \code{a} (threshold separation), \code{z} (relative starting point), \code{v} (drift rate), \code{t0} (non-decision time/response time constant), \code{d} (differences in speed of response execution), \code{sv} (inter-trial-variability of drift), \code{st0} (inter-trial-variability of non-decisional components), and \code{sz} (inter-trial-variability of relative starting point). 
#'
#' @param t a vector of RTs.
#' @param n desired number of observations.
#' @param boundary character vector. Which boundary should be tested. Possible values are \code{c("upper", "lower")}, possibly abbreviated and \code{"upper"} being the default.
#' 
#' @param a threshold separation. Amount of information that is considered for a decision. Large values indicate a conservative decisional style. Typical range: 0.5 < \code{a} < 2
#' @param v drift rate. Average slope of the information accumulation process. The drift gives information about the speed and direction of the accumulation of information. Large (absolute) values of drift indicate a good performance. If received information supports the response linked to the upper threshold the sign will be positive and vice versa. Typical range: -5 < \code{v} < 5
#' @param t0 non-decision time or response time constant (in seconds). Lower bound for the duration of all non-decisional processes (encoding and response execution). Typical range: 0.1 < \code{t0} < 0.5
#' 
#' @param z relative starting point. Indicator of an a priori bias in decision making. When the relative starting point \code{z} deviates from 0.5, the amount of information necessary for a decision differs between response alternatives. Typical range: 0.3 < \code{z} < 0.7. Default is 0.5 (i.e., no bias).
#' @param d differences in speed of response execution (in seconds). Positive values indicate that response execution is faster for responses linked to the upper threshold than for responses linked to the lower threshold. Typical range: -0.1 < \code{d} < 0.1. Default is 0.
#' @param sz inter-trial-variability of (relative) starting point. Range of a uniform distribution with mean \code{z} describing the distribution of actual starting points from specific trials. Minimal impact on the RT distributions. Can be fixed to 0 in most applications. Typical range: 0 < \code{sz} < 0.5. Default is 0.
#' @param sv inter-trial-variability of drift rate. Standard deviation of a normal distribution with mean \code{v} describing the distribution of actual drift rates from specific trials. 	Minimal impact on the RT distributions. Can be fixed to 0 in most applications. Typical range: 0 < \code{sv} < 2. Default is 0.
#' @param st0 inter-trial-variability of non-decisional components. Range of a uniform distribution with mean \code{t0 + st0/2} describing the distribution of actual \code{t0} values across trials. Accounts for response times below \code{t0}. Reduces skew of predicted RT distributions. Typical range: 0 < \code{st0} < 0.2. Default is 0.
#' 
#' @param precision \code{numerical} scalar value. Precision of calculation. Corresponds roughly to the number of decimals of the predicted CDFs that are calculated accurately. Default is 3.
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
#' @name Diffusion
#' 
#' @example examples/examples.diffusion.R
#' 
NULL

#pl <- get_rd_parameters(a = a, z = z, v = v, t0 = t0, d = d, sz = sz, sv = sv, st0 = st0, pl = parameters)
# get_rd_parameters <- function(a, z, v, t0, d, sz, sv, st0, pl) {
#   miss_p_indiv <- c(missing(a), missing(v), missing(t0), missing(z), missing(d), missing(sz), missing(sv), missing(st0))
#   miss_pl <- missing(pl)
#   if(all(miss_p_indiv) & miss_pl) stop("No parameter values specified.")
#   if(any(!miss_p_indiv) & !miss_pl) stop("Parameter values need to be specified either via individual arguments or via 'pl', but not both.")
#   if((sum(miss_p_indiv) > 0) & miss_pl) stop("Not all parameter values specified via individual arguments.")
#   if(miss_pl) {
#     rd_parameters <- c(a,v,t0,d,sz,sv,st0,z)
#     if(length(rd_parameters)!=8) stop("Each parameter needs to be of length 1.")
#     if(!is.numeric(rd_parameters)) stop("Parameters need to be numeric.")
#     if (any(is.na(rd_parameters)) || !all(is.finite(rd_parameters))) stop("Parameters need to be numeric and finite.")
#     return(rd_parameters)
#   } else {
#     if (!is.numeric(pl) | !is.vector(pl)) stop("pl needs to be a numeric vector.")
#     if(length(pl)!=8) stop("Length of pl is incorrect (must be 8).")
#     if (any(is.na(pl)) || !all(is.finite(pl))) stop("Parameters need to be numeric and finite.")
#     if (!is.null(names(pl))) {
#       if(!all(sapply(c("a","v","t0","d","sz","sv","st0","z"), function(x) x %in% names(pl)))) stop("names(pl) not complete. Need to be: c(\"a\",\"v\",\"t0\",\"z\",\"d\",\"sz\",\"sv\",\"st0\")")
#       pl <- pl[c("a","v","t0","d","sz","sv","st0","z")]
#     } else {
#       names(pl) <- c("a","z","v","t0","d","sz","sv","st0")
#       pl <- pl[c("a","v","t0","d","sz","sv","st0","z")]
#     }
#     return(pl)
#   }
# }
#





# [MG 20150616]
# In line with LBA, adjust t0 to be the lower bound of the non-decision time distribution rather than the average 
# Called from prd, drd, rrd 
recalc_t0 <- function (t0, st0) { t0 <- t0 + st0/2 }


#' @rdname Diffusion
#' @export drd_batch
#' [TODO]
drd_batch <- function (t, boundary = c("upper", "lower"), 
                 a, v, t0, z = 0.5, d = 0, sz = 0, sv = 0, st0 = 0, 
                 precision = 3)
{
  t0 <- recalc_t0 (t0, st0) 
  
  # Convert boundaries to numeric (should find a less clumsy way, but 'replace' sucks)
  numeric_bounds <- vector (mode="numeric")
  for (i in 1:length(boundary)) { 
    if (boundary[i] == "upper") numeric_bounds[i] <- 2L
    if (boundary[i] == "lower") numeric_bounds[i] <- 1L
  }
  
  # Build parameter matrix
  params <- cbind (a, v, t0, z, d, sz, sv, st0, numeric_bounds)
  
  # Check for illegal parameter values
  if(any(missing(a), missing(v), missing(t0))) stop("a, v, and/ot t0 must be supplied")
  if(!is.numeric(params)) stop("Parameters need to be numeric.")
  if (any(is.na(params)) || !all(is.finite(params))) stop("Parameters need to be numeric and finite.")

  densities <- vector(length=length(t))    
  
  if (length(params[,1]) == 1)  # Treat this like a parameter vector
  {
    pl <- c(a,v,t0,d,sz,sv,st0,z)
    
    # Call the C code
    output <- .C("dfastdm_b", 
                 as.integer (length(t)),                 # 1  IN:  number of densities
                 as.vector  (pl),                        # 2  IN:  parameters
                 as.vector  (t),                         # 3  IN:  RTs
                 as.double  (precision),                 # 4  IN:  precision
                 as.integer (numeric_bounds),                  # 5  IN:  boundary
                 as.vector  (densities, mode="numeric")    # 6 OUT:  densities
    )
    
    unlist(output[6])
    
  } else
  {
    # Check matrix-specific errors
    #  Needs to be 2D, with second dimension of length 9 (incl. boundary)
    if (length(dim(params)) != 2)   { stop("Needs to be a two-dimensional parameter matrix.") }
    if (length (params[1,]) != 9)   { stop("Incorrect length of parameter rows") }
    
    uniques <- unique(params)
    unique_rows <- length(uniques[,1])
    
    # ?TODO: If we wanted to do an optimising heuristic with the percent of unique rows, we'd do it here
    #        e.g. if >95% rows are unique is probably slower to batch them than send them off individually
    
    for (i in 1:unique_rows)
    {
      # Get list of row indices equal to the i'th unique row
      # Todo: check this actually works!
      ok_rows <- which(params[,1]      == uniques[i,1])
      for (j in 2:9) { ok_rows <- ok_rows[which(params[ok_rows,j] == uniques[i,j])] }    
      
      # Select the correct RT indices and parameters for this batch of 'ok' parameter rows
      ok_params <- params[ok_rows[1],]
      
      # pl <- c(a,v,t0,d,sz,sv,st0,z) (need to put z at the end)
      pl <- c(ok_params[1], ok_params[2], ok_params[3], ok_params[5], 
              ok_params[6], ok_params[7], ok_params[8], ok_params[4])  
      
      # Call the C code
      output <- .C("dfastdm_b", 
                   as.integer (length(t[ok_rows])),        # 1  IN:  number of densities
                   as.vector  (pl),                        # 2  IN:  parameters
                   as.vector  (t[ok_rows]),                # 3  IN:  RTs
                   as.double  (precision),                 # 4  IN:  precision
                   as.integer (ok_params[9]),                # 5  IN:  boundary
                   as.vector  (densities, mode="numeric")  # 6 OUT:  densities
      )
      densities[ok_rows] <- head(unlist(output[6]), length(ok_rows))
    }
    densities
  }
  
  
}

#' @rdname Diffusion
#' @export prd_batch
# [TODO]
prd_batch <- function (t, boundary = c("upper", "lower"), 
                 a, v, t0, z = 0.5, d = 0, sz = 0, sv = 0, st0 = 0, 
                 precision = 3, maxt = 1e4) 
{
  t0 <- recalc_t0 (t0, st0) 

  t[t>maxt] <- maxt
  if(!all(t == sort(t)))  stop("t needs to be sorted")

  # Convert boundaries to numeric (should find a less clumsy way, but 'replace' sucks)
  numeric_bounds <- vector (mode="numeric")
  for (i in 1:length(boundary)) { 
    if (boundary[i] == "upper") numeric_bounds[i] <- 2L
    if (boundary[i] == "lower") numeric_bounds[i] <- 1L
  }
  
  # Build parameter matrix
  params <- cbind (a, v, t0, z, d, sz, sv, st0, numeric_bounds)
  
  # Check for illegal parameter values
  if(any(missing(a), missing(v), missing(t0))) stop("a, v, and/ot t0 must be supplied")
  if(!is.numeric(params)) stop("Parameters need to be numeric.")
  if (any(is.na(params)) || !all(is.finite(params))) stop("Parameters need to be numeric and finite.")
  
  pvalues <- vector(length=length(t))    
  
  if (length(params[,1]) == 1)  # Treat this like a parameter vector
  {
    pl <- c(a,v,t0,d,sz,sv,st0,z)
    
    # Call the C code
    output <- .C("pfastdm_b", 
                 as.integer (length(t)),                 # 1  IN:  number of densities
                 as.vector  (pl),                        # 2  IN:  parameters
                 as.vector  (t),                         # 3  IN:  RTs
                 as.double  (precision),                 # 4  IN:  precision
                 as.integer (numeric_bounds),                  # 5  IN:  boundary
                 as.vector  (pvalues, mode="numeric")    # 6 OUT:  densities
    )
    
    unlist(output[6])
    
  } else
  {
    # Check matrix-specific errors
    #  Needs to be 2D, with second dimension of length 9 (incl. boundary)
    if (length(dim(params)) != 2)   { stop("Needs to be a two-dimensional parameter matrix.") }
    if (length (params[1,]) != 9)   { stop("Incorrect length of parameter rows") }
    
    uniques <- unique(params)
    unique_rows <- length(uniques[,1])
    
    # ?TODO: If we wanted to do an optimising heuristic with the percent of unique rows, we'd do it here
    #        e.g. if >95% rows are unique is probably slower to batch them than send them off individually
    
    for (i in 1:unique_rows)
    {
      # Get list of row indices equal to the i'th unique row
      # Todo: check this actually works!
      ok_rows <- which(params[,1]      == uniques[i,1])
      for (j in 2:9) { ok_rows <- ok_rows[which(params[ok_rows,j] == uniques[i,j])] }    
      
      # Select the correct RT indices and parameters for this batch of 'ok' parameter rows
      ok_params <- params[ok_rows[1],]
      
      # pl <- c(a,v,t0,d,sz,sv,st0,z) (need to put z at the end)
      pl <- c(ok_params[1], ok_params[2], ok_params[3], ok_params[5], 
              ok_params[6], ok_params[7], ok_params[8], ok_params[4])  
      
      # Call the C code
      output <- .C("pfastdm_b", 
                   as.integer (length(t[ok_rows])),        # 1  IN:  number of densities
                   as.vector  (pl),                        # 2  IN:  parameters
                   as.vector  (t[ok_rows]),                # 3  IN:  RTs
                   as.double  (precision),                 # 4  IN:  precision
                   as.integer (ok_params[9]),                # 5  IN:  boundary
                   as.vector  (pvalues, mode="numeric")  # 6 OUT:  densities
      )
      pvalues[ok_rows] <- head(unlist(output[6]), length(ok_rows))
    }
    pvalues
  }
}


#' @rdname Diffusion
#' @export rrd_batch
#' When given vectorised parameters, n is the number of replicates for each parameter set
rrd_batch <- function (n, 
                 a, v, t0, z = 0.5, d = 0, sz = 0, sv = 0, st0 = 0, 
                 precision = 3)
{

  #t0 <- recalc_t0 (t0, st0) 
  
  # Build parameter matrix
  params <- cbind (a, v, t0, z, d, sz, sv, st0)
  
  params[3,] <- recalc_t0 (params[3,], params[8,]) 
  
  
  # Check for illegal parameter values
  if(any(missing(a), missing(v), missing(t0))) stop("a, v, and/ot t0 must be supplied")
  if(!is.numeric(params)) stop("Parameters need to be numeric.")
  if (any(is.na(params)) || !all(is.finite(params))) stop("Parameters need to be numeric and finite.")
  
  randRTs    <- vector(length=n*length(params[,1]))
  randBounds <- vector(length=n*length(params[,1]))

  if (length(params[,1]) == 1)  # Treat this like a parameter vector
  {
    pl <- c(a,v,t0,d,sz,sv,st0,z)
  
    output <- .C("rfastdm", 
                 as.integer (n),                          # 1  IN:  number of densities
                 as.vector  (pl),                         # 2  IN:  parameters
                 as.double  (precision),                  # 3  IN:  precision
                 as.vector  (randRTs, mode="numeric"),    # 4 OUT:  RTs 
                 as.vector  (randBounds, mode="numeric")  # 5 OUT:  bounds 
    )
    randRTs <- unlist(output[4]) 
    randBounds <- unlist(output[5])
  } else
  {
    # Check matrix-specific errors
    #  Needs to be 2D, with second dimension of length 9 (incl. boundary)
    if (length(dim(params)) != 2)   { stop("Needs to be a two-dimensional parameter matrix.") }
    if (length (params[1,]) != 8)   { stop("Incorrect length of parameter rows") }
    uniques <- unique(params)
    unique_rows <- length(uniques[,1])
    
    # ?TODO: If we wanted to do an optimising heuristic with the percent of unique rows, we'd do it here
    #        e.g. if >95% rows are unique is probably slower to batch them than send them off individually
  
    start_idx <-1
    for (i in 1:unique_rows)
    {
      
      # Get list of row indices equal to the i'th unique row
      # Todo: check this actually works!
      ok_rows <- which(params[,1]      == uniques[i,1])
      for (j in 2:8) { ok_rows <- ok_rows[which(params[ok_rows,j] == uniques[i,j])] }    
    
      # Calculate n for this row
      current_n <- n * length(ok_rows)
      
      # Select the correct RT indices and parameters for this batch of 'ok' parameter rows
      ok_params <- params[ok_rows[1],]
      
      # pl <- c(a,v,t0,d,sz,sv,st0,z) (need to put z at the end)
      pl <- c(ok_params[1], ok_params[2], ok_params[3], ok_params[5], 
              ok_params[6], ok_params[7], ok_params[8], ok_params[4])  
      
      # Call the C code
      output <- .C("rfastdm", 
                   as.integer (current_n),                 # 1  IN:  number of densities
                   as.vector  (pl),                        # 2  IN:  parameters
                   as.double  (precision),                 # 3  IN:  precision
                   as.vector  (randRTs, mode="numeric"),   # 4 OUT:  RTs 
                   as.vector  (randBounds, mode="numeric") # 5 OUT:  bounds 
      )
      
      end_idx <- start_idx+current_n-1              

      randRTs[start_idx:end_idx]    <- head(unlist(output[4]), current_n) 
      randBounds[start_idx:end_idx] <- head(unlist(output[5]), current_n)
      start_idx <- start_idx + current_n
    }
  }
  response <- factor(randBounds, levels = 0:1, labels = c("lower", "upper"))
  data.frame(rt = randRTs, response)
}

