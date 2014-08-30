#' The Ratcliff Diffusion Model
#' 
#' Density, distribution function, and random generation for the Ratcliff diffusion model with eight parameters: \code{a} (threshold separation), \code{v} (drift rate), \code{t0} (non-decision time/response time constant), \code{z} (relative starting point), \code{d} (differences in speed of response execution), \code{sv} (inter-trial-variability of drift), \code{st0} (inter-trial-variability of non-decisional components), and \code{sz} (inter-trial-variability of relative starting point).
#'
#'@param t a vector of RTs.
#'@param n desired number of observations.
#'@param pl a list of paramaters in the following order: \code{c("a","v","t0","d","sz","sv","st0","z")}
#'@param i the boundary to test, upper = 2, lower = 1
#'@param precision \code{numerical} scalar value. Precision of calculation. Corresponds roughly to the number of decimals of the predicted CDFs that are calculated accuratly. Default is 3.
#'@param maxt maximum \code{t0} value, used to stop integration problems (\code{prd} only).
#'
#'@return \code{drd} gives the density, \code{prd} gives the dsitribution function, and \code{rrd} generates random response times and decisions (in a \code{matrix}).
#'
#'@author Underlying C code by Jochen Voss and Andreas Voss. Porting and R wrapping by Matthew Gretton and Andrew Heathcote.
#'
#' @useDynLib rtdists, dfastdm_b, pfastdm_b, rfastdm
#'
#' @name diffusion
#' 
NULL


#' @rdname diffusion
drd <- function (t, pl, i, precision = 3)
{
  # Check for illegal parameter values
  if (any(is.na(pl)) || !all(is.finite(pl))) return(rep(NA,length=length(t)))
  
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
# set maximum t value to stop integration problems
prd <- function (t, pl, i, precision = 3, maxt = 1e4) 
{
  # Check for illegal parameter values
  pvalues <- rep(NA,length=length(t))
  if (any(is.na(pl)) || !all(is.finite(pl))) return(pvalues)
  t[t>maxt] <- maxt
  
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
# Returns a matrix of 2 x n (RTs x boundaries)
rrd <- function (n, pl, precision = 3)
{
  randRTs    <- vector(length=n)
  randBounds <- vector(length=n)
  
  output <- .C("rfastdm", 
               as.integer (n),                          # 1  IN:  number of densities
               as.vector  (pl),                         # 2  IN:  parameters
               as.double  (precision),                  # 3  IN:  precision
               as.vector  (randRTs, mode="numeric"),    # 4 OUT:  RTs 
               as.vector  (randBounds, mode="numeric")  # 5 OUT:  bounds 
  )
  
  randRTs <- unlist(output[4]) 
  randBounds <- unlist(output[5])
  
  # Note: incrementing the bounds for 0,1 to be 1,2
  matrix (c(randRTs, randBounds+1), ncol=2)
}

