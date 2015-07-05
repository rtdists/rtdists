# rfastdm vectorise test

library (devtools)

devtools::install_git ("https://github.com/rtdists/rtdists.git")
library (rtdists)

source ("_dev/diffusion_batch_test.R")

# VECTOR VERSION: Set up parameters
test_vec_len <- 1
vec_rts <- seq (1, 2.9, by=0.2)
vec_bounds <- "upper"

vec_a   <- 1
vec_z   <- .20 
vec_v   <- .30
vec_t0  <- .40
vec_d   <- .50
vec_sz  <- .06
vec_sv  <- .70
vec_st0 <- .80

# MATRIX VERSION: Set up parameter vectors
test_vec_len <- 10
vec_rts <- seq (1, 2.9, by=0.2)
vec_bounds <- sample (c("upper","lower"), test_vec_len, replace=TRUE)

vec_a   <- c(1,    1,    1,    2,    3,    3,    3,    4,    4,    4    )
vec_z   <- rep (.20, test_vec_len)
vec_v   <- rep (.30, test_vec_len)
vec_t0  <- rep (.40, test_vec_len)
vec_d   <- rep (.50, test_vec_len)
vec_sz  <- rep (.06, test_vec_len)
vec_sv  <- rep (.70, test_vec_len)
vec_st0 <- rep (.80, test_vec_len)


# !! DEBUG ONLY
# Get known right answers
#
correct_prd_vals <- vector(length=test_vec_len)
correct_drd_vals <- vector(length=test_vec_len)
for (i in 1:test_vec_len)
{
  correct_prd_vals[i] <- prd(vec_rts[i], boundary=vec_bounds[i], a =vec_a[i], z=vec_z[i], v =vec_v[i], 
                             t0=vec_t0[i], d=vec_d[i], sz =vec_sz[i], sv=vec_sv[i], st0 = vec_st0[i])
  correct_drd_vals[i] <- drd(vec_rts[i], boundary=vec_bounds[i], a =vec_a[i], z=vec_z[i], v =vec_v[i], 
                     t0=vec_t0[i], d=vec_d[i], sz =vec_sz[i], sv=vec_sv[i], st0 = vec_st0[i])
}  


# drd (vec_rts, boundary=vec_bounds, vec_a, vec_z, vec_v, vec_t0, vec_d, vec_sz, vec_sv, vec_st0)
t        <- vec_rts
boundary <- vec_bounds

a        <- vec_a
v        <- vec_v
t0       <- vec_t0
z        <- vec_z
d        <- vec_d
sz       <- vec_sz
sv       <- vec_sv
st0      <- vec_st0


### START OF FUNCTION

# Checks to add:
#   1. Are all vectors the same length?
#   2. If not, what do we do? 
#      ? Cycle parameter vectors over rts vector?
#      ? Allow some params to be scalar? (we'll need to use something as a canonical parameter length?)


# Build the parameter matrix

# Convert boundaries to numeric
numeric_bounds <- vector (mode="numeric")
for (i in 1:length(boundary)) { 
  if (boundary[i] == "upper") numeric_bounds[i] <- 2L
  if (boundary[i] == "lower") numeric_bounds[i] <- 1L
}

params <- cbind (a, v, t0, z, d, sz, sv, st0, numeric_bounds)

# Check for illegal parameter values
if(any(missing(a), missing(v), missing(t0))) stop("a, v, and/ot t0 must be supplied")
if(!is.numeric(params)) stop("Parameters need to be numeric.")
if (any(is.na(params)) || !all(is.finite(params))) stop("Parameters need to be numeric and finite.")

densities <- vector(length=length(t))    

if (length(params[,1]) == 1)
{
  # Call the C code
  output <- .C("dfastdm_b", 
               as.integer (length(t)),                 # 1  IN:  number of densities
               as.vector  (pl),                        # 2  IN:  parameters
               as.vector  (t),                         # 3  IN:  RTs
               as.double  (precision),                 # 4  IN:  precision
               as.integer (i),                         # 5  IN:  boundary
               as.vector  (densities, mode="numeric")  # 6 OUT:  densities
  )
  
  unlist(output[6])
  
} else
{
  cat ("Doing vector params,\n")
  # Check matrix-specific errors
  #  Needs to be 2D, with second dimension of length 9 (incl. boundary)
  if (length(dim(params)) != 2)   { stop("Needs to be a two-dimensional parameter matrix.") }
  if (length (params[1,]) != 9)   { stop("Incorrect length of parameter rows") }

  uniques <- unique(params)
  unique_rows <- length(uniques[,1])

  # ?TODO: If we wanted to do an optimising heuristic with the percent of unique rows, we'd do it here
 
  for (i in 1:unique_rows)
  {
    cat ("Processing unique row", i, "\n")

    # Get list of row indices equal to the i'th unique row
    # Todo: check this actually works!
    ok_rows <- which(params[,1]      == uniques[i,1])
    for (j in 2:9) { ok_rows <- ok_rows[which(params[ok_rows,j] == uniques[i,j])] }    
    
    # Select the correct RT indices and parameters for this batch of 'ok' parameter rows
    ok_params <- params[ok_rows[1],]

    cat ("Rows to process:", ok_rows, "\n")
    cat ("RTs to process:", t[ok_rows], "\n")
    cat ("Params to use:", ok_params, "\n")
    
    # pl <- c(a,v,t0,d,sz,sv,st0,z) (need to put z at the end)
    pl <- c(ok_params[1], ok_params[2], ok_params[3], ok_params[5], 
            ok_params[6], ok_params[7], ok_params[8], ok_params[4])  

    # Call the C code
    output <- .C("dfastdm_b", 
                 as.integer (length(rts[ok_rows])),                 # 1  IN:  number of densities
                 as.vector  (pl),                        # 2  IN:  parameters
                 as.vector  (t[ok_rows]),                         # 3  IN:  RTs
                 as.double  (precision),                 # 4  IN:  precision
                 as.integer (ok_paramsi),                         # 5  IN:  boundary
                 as.vector  (densities, mode="numeric")  # 6 OUT:  densities
    )
    densities[ok_rows] <- out

  
#     out2 <- drd (rts[ok_rows], a =current_params[1], z  =current_params[2], v =current_params[3], 
#                      t0=current_params[4], d  =current_params[5], sz=current_params[6], 
#                      sv=current_params[7], st0=current_params[8])
#     
#     output2[ok_rows] <- out2
  }
  densities
  
  # Test
  all.equal (correct_prd_vals, densities)
  
}

