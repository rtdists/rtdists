# Rfastdm Wrapper

## REMOVE THESE 3 LINES FOR PRODUCTION
#system ("R CMD SHLIB Rfastdm2.c density.c cdf.c pde.c phi.c precision.c") 


# For these functions:
#   t  - a vector of RTs
#   pl - a list of parameters
#        order: c("a","v","ter","d","sz","sv","st","z")
#   i  - the boundary to test, upper = 2, lower = 1
#   n  - number of samples to produce
drd <- function (t, pl, i, precision = 3)
{
    # Check for illegal parameter values
     if (any(is.na(pl)) || !all(is.finite(pl))) return(rep(NA,length=length(t)))
       
    # Call the C code
    densities <- vector(length=length(t))    
#    dyn.load ("Rfastdm2.so") # Load library
    output <- .C("dfastdm_b", 
                 as.integer (length(t)),                 # 1  IN:  number of densities
                 as.vector  (pl),                        # 2  IN:  parameters
                 as.vector  (t),                         # 3  IN:  RTs
                 as.double  (precision),                 # 4  IN:  precision
                 as.integer (i),                         # 5  IN:  boundart 
                 as.vector  (densities, mode="numeric")  # 6 OUT:  densities
    )
#    dyn.unload ("Rfastdm2.so") # Unload library

   unlist(output[6])
  
}

# set maximum t value to stop integration problems
prd <- function (t, pl, i, precision = 3, maxt = 1e4) 
{
  # Check for illegal parameter values
  pvalues <- rep(NA,length=length(t))
  if (any(is.na(pl)) || !all(is.finite(pl))) return(pvalues)
  t[t>maxt] <- maxt
  
  # Call the C code
#  dyn.load ("Rfastdm2.so") # Load library
  pvalues <- vector(length=length(t))    
  output <- .C("pfastdm_b", 
               as.integer (length(t)),               # 1  IN:  number of densities
               as.vector  (pl),                      # 2  IN:  parameters
               as.vector  (t),                       # 3  IN:  RTs
               as.double  (precision),               # 4  IN:  number of densities
               as.integer (i),                       # 5  IN:  boundary 
               as.vector  (pvalues, mode="numeric")  # 6 OUT:  pvalues
  )
#  dyn.unload ("Rfastdm2.so") # Unload library
  
  unlist(output[6])
  
}

# Returns a matrix of 2 x n (RTs x boundaries)
rrd <- function (n, pl, precision = 3)
{
  randRTs    <- vector(length=n)
  randBounds <- vector(length=n)
  
#  dyn.load ("Rfastdm2.so") # Load library
  output <- .C("rfastdm", 
               as.integer (n),                          # 1  IN:  number of densities
               as.vector  (pl),                         # 2  IN:  parameters
               as.double  (precision),                  # 3  IN:  precision
               as.vector  (randRTs, mode="numeric"),    # 4 OUT:  RTs 
               as.vector  (randBounds, mode="numeric")  # 5 OUT:  bounds 
  )
#  dyn.unload ("Rfastdm2.so") # Unload library
  
  randRTs <- unlist(output[4]) 
  randBounds <- unlist(output[5])
  
  # Note: incrementing the bounds for 0,1 to be 1,2
  matrix (c(randRTs, randBounds+1), ncol=2)
}

