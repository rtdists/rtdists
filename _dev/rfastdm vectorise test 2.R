# rfastdm vectorise test

library (devtools)

devtools::install_git ("https://github.com/rtdists/rtdists.git")
library (rtdists)


#' @rdname Diffusion
#' @export prd
# set maximum t value to stop integration problems
prd_test <- function (t, boundary = c("upper", "lower"), 
                 a, v, t0, z = 0.5, d = 0, sz = 0, sv = 0, st0 = 0, 
                 precision = 3, maxt = 1e4) 
{
  t0 <- recalc_t0 (t0, st0) # [MG 20150616] Call recalc_t0
  
  # Check for illegal parameter values
  if(any(missing(a), missing(v), missing(t0))) stop("a, v, and/ot t0 must be supplied")
  pl <- c(a,v,t0,d,sz,sv,st0,z)
  
  length(pl)
  
  if(length(pl)!=8) stop("Each parameter needs to be of length 1.")
  if(!is.numeric(pl)) stop("Parameters need to be numeric.")
  if (any(is.na(pl)) || !all(is.finite(pl))) stop("Parameters need to be numeric and finite.")
  
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


params <- matrix (nrow=10, ncol=8)
#                   a    z      v    t0     d     sz    sv   st0
params[ 1,] <- c( 1,   0.2,   0.3,  0.4,  0.5,  0.06,  0.7,  0.8)
params[ 2,] <- c( 1,   0.2,   0.3,  0.4,  0.5,  0.06,  0.7,  0.8)
params[ 3,] <- c( 1,   0.2,   0.3,  0.4,  0.5,  0.06,  0.7,  0.8)
params[ 4,] <- c( 2,   0.2,   0.3,  0.4,  0.5,  0.06,  0.7,  0.8)
params[ 5,] <- c( 3,   0.2,   0.3,  0.4,  0.5,  0.06,  0.7,  0.8)
params[ 6,] <- c( 3,   0.2,   0.3,  0.4,  0.5,  0.06,  0.7,  0.8)
params[ 7,] <- c( 3,   0.2,   0.3,  0.4,  0.5,  0.06,  0.7,  0.8)
params[ 8,] <- c( 4,   0.2,   0.3,  0.4,  0.5,  0.06,  0.7,  0.8)
params[ 9,] <- c( 4,   0.2,   0.3,  0.4,  0.5,  0.06,  0.7,  0.8)
params[10,] <- c( 4,   0.2,   0.3,  0.4,  0.5,  0.06,  0.7,  0.8)

rts <- seq (1, 2.9, by=0.2)

# Get known right answers
prd_vals <- vector(length=10)
drd_vals <- vector(length=10)
for (i in 1:10)
{
  prd_vals[i] <- prd(rts[i], a =params[i,1], z=params[i,2], v =params[i,3], 
                             t0=params[i,4], d=params[i,5], sz =params[i,6], 
                             sv=params[i,7], st0 = params[i,8])
  drd_vals[i] <- drd(rts[i], a =params[i,1], z=params[i,2], v =params[i,3], 
                             t0=params[i,4], d=params[i,5], sz =params[i,6], 
                             sv=params[i,7], st0 = params[i,8])
}  
  


### START OF FUNCTION

if (dim (params) == NULL)
{
   # This is just a single vector of parameters 
  output <- prd_test (rts, a =params[1], z=params[2], v =params[3], 
                           t0=params[4], d=params[5], sz =params[6], 
                           sv=params[7], st0 = params[8])
}
else
{
  # Check matrix-specific errors
  #  Needs to be 2D, with second dimension of length 8
  if (length(dim(params)) != 2)   { stop("Bugger 1.") }
  if (length (params[1,]) != 8)   { stop("Bugger 2.") }

  uniques <- unique(params)
  unique_rows <- length(uniques[,1])

  # ?TODO: If we wanted to do an optimising heuristic with the percent of unique rows, we'd do it here
 
  # Set up a vector of output values
  output <- vector (length=length(params[,1]))
#   output2 <- vector (length=length(params[,1]))
  
  for (i in 1:unique_rows)
  {
    cat ("Processing unique row", i, "\n")

    # Get list of row indices equal to the i'th unique row
    # Todo: check this actually works!
    ok_rows <- which(params[,1]      == uniques[i,1])
    for (j in 2:8) { ok_rows <- ok_rows[which(params[ok_rows,j] == uniques[i,j])] }    
    
    # Select the correct RT indices and parameters for this batch of 'ok' parameter rows
    current_params <- params[ok_rows[1],]

    cat ("Rows to process:", ok_rows, "\n")
    cat ("RTs to process:", rts[ok_rows], "\n")
    cat ("Params to use:", current_params, "\n")
    out <- prd_test (rts[ok_rows], a =current_params[1], z  =current_params[2], v =current_params[3], 
                     t0=current_params[4], d  =current_params[5], sz=current_params[6], 
                     sv=current_params[7], st0=current_params[8])
      
    output[ok_rows] <- out

  
#     out2 <- drd (rts[ok_rows], a =current_params[1], z  =current_params[2], v =current_params[3], 
#                      t0=current_params[4], d  =current_params[5], sz=current_params[6], 
#                      sv=current_params[7], st0=current_params[8])
#     
#     output2[ok_rows] <- out2
  }
  output
  
  # Test
  all.equal (prd_vals, output)
  
}

