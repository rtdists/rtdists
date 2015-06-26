# rfastdm vectorise test

library (devtools)

devtools::install_local ("rtdists-master (120615)/rtdists_0.3-1.tar.gz")
library (rtdists)


#' @rdname Diffusion
#' @export prd
# set maximum t value to stop integration problems
prd_batch <- function (t, boundary = c("upper", "lower"), 
                 a, v, t0, z = 0.5, d = 0, sz = 0, sv = 0, st0 = 0, 
                 precision = 3, maxt = 1e4) 
{
  t0 <- recalc_t0 (t0, st0) # [MG 20150616] Call recalc_t0
  
  
  
  
  param_order <- order (params[,1], params[,2], 
                        params[,3], params[,4], 
                        params[,5], params[,6], 
                        params[,7], params[,8])
  
  sorted_params <- params[param_order,];
  
  sorted_params <- cbind (sorted_params, duplicated (sorted_params))  # Flag duplicates
  sorted_params <- cbind (sorted_params, param_order)                 # Add order so we can re-sort back to original
  
  
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



params <- matrix (nrow=25, ncol=8)
params[ 1,] <- c(41,  42,  43, 44, 45, 46, 47, 48)
params[ 2,] <- c(51,  52,  53, 54, 55, 56, 57, 58)
params[ 3,] <- c(31,  32,  33, 34, 35, 36, 37, 38)
params[ 4,] <- c(11,  12,  13, 14, 15, 16, 17, 18)
params[ 5,] <- c(21,  22,  23, 24, 25, 26, 27, 28)
params[ 6,] <- c(61,  62,  63, 64, 65, 66, 67, 68)
params[ 7,] <- c(11, 152,  13, 14, 15, 16, 17, 18)
params[ 8,] <- c(11, 132, 135, 14, 15, 16, 17, 18)
params[ 9,] <- c(11, 132, 135, 14, 15, 16, 17, 11)
params[10,] <- c(11,  12, 135, 14, 15, 16, 17, 18)
params[11,] <- c(11,  12,  13, 14, 15, 16, 17, 18)
params[12,] <- c(11, 132, 135, 14, 15, 16, 17, 11)
params[13,] <- c(11.0423, 132.323, 135.64563542899, 14.57490767, 15.3424565676798905654, 16.34534534, 0.000000534543, 0.0000000000000000004)
params[14,] <- c(11.0423, 132.323, 135.64563542899, 14.57490767, 15.3424565676798905654, 16.34534534, 0.000000534543, 0.0000000000000000004)
params[15,] <- c(11.0423, 132.323, 135.64563542899, 14.57490767, 15.3424565676798905654, 16.34534534, 0.000000534543, 0.0000000000000000004)
params[16,] <- c(11.0423, 132.323, 135.64563542899, 14.57490767, 15.3424565676798905654, 16.34534534, 0.000000534543, 0.0000000000000000004)
params[17,] <- c(11.0423, 132.323, 135.64563542899, 14.57490767, 15.3424565676798905654, 16.34534534, 0.000000534543, 0.0000000000000000004)
params[18,] <- c(11.0423, 132.323, 135.64563542899, 14.57490767, 15.3424565676798905654, 16.34534534, 0.000000534543, 0.0000000000000000004)
params[19,] <- c(11.0423, 132.323, 135.64563542899, 14.57490767, 15.3424565676798905654, 16.34534534, 0.000000534543, 0.0000000000000000004)
params[20,] <- c(11.0423, 132.323, 135.64563542899, 14.57490767, 15.3424565676798905654, 16.34534534, 0.000000534543, 0.0000000000000000004)
params[21,] <- c(61,  62,  63, 64, 65, 66, 67, 68)
params[22,] <- c(11, 152,  13, 14, 15, 16, 17, 18)
params[23,] <- c(11, 132, 135, 14, 15, 16, 17, 18)
params[24,] <- c(11, 132, 135, 14, 15, 16, 17, 11)
params[25,] <- c(61,  62,  63, 64, 65, 66, 67, 68)
rts <- runif (25, 0.1, 3.0)

# big_params <- rbind(params, params)
# big_params <- rbind(big_params, big_params)
# big_params <- rbind(big_params, big_params)
# big_params <- rbind(big_params, big_params)
# big_params <- rbind(big_params, big_params)
# big_params <- rbind(big_params, big_params)
# big_params <- rbind(big_params, big_params)
# big_params <- rbind(big_params, big_params)
# big_params <- rbind(big_params, big_params)
# #for (i in 1:27200) { big_params <- rbind (big_params, sample (200, 8)) }
# for (i in 1:87200) { big_params <- rbind (big_params, c(61,  62,  63, 64, 65, 66, 67, 68)) }
# params <- big_params


# For checking
# param_order <- order (params[,1], params[,2], 
#                        params[,3], params[,4], 
#                        params[,5], params[,6], 
#                        params[,7], params[,8])
# sorted_params <- params[param_order,];

uniques <- unique(params)
unique_rows <- length(uniques[,1])
#duplicated (params) 

# ok <- params == params[8,]
# 
# for (i in 1:unique_rows)
# {
#   ok <- params[,1] == params[i]
#   
#   sub_params <- params[ok,]
# }
# 
# sorted_params <- cbind (sorted_params, duplicated (sorted_params))  # Flag duplicates
# sorted_params <- cbind (sorted_params, param_order)                 # Add order so we can re-sort back to original
# 
# test <- prd_batch (1, a=params[,1], z=params[,2], v=params[,3], t0=params[,4], d=params[,5], sz=params[,6], sv=params[,7], st0=params[,8])
# 
# x<-paste (c(0.453, 24.24352, 534.345, 2.245744,0.342, 0.2432), collapse='')
# 


system.time ({
#### Stepwise parameter compare
  # 100000 rows, ~90%+ Unique: 110.25 secs
  # 100000 rows, ~25%+ Unique: 45.56 secs
  # 100000 rows, ~<1%  Unique:  1.39 secs
uniques <- unique(params)
unique_rows <- length(uniques[,1])

for (i in 1:unique_rows)
{
  ok = which(params[,1] == uniques[i,1])
  ok = ok[which(params[ok,2] == uniques[i,2])]
  ok = ok[which(params[ok,3] == uniques[i,3])]
  ok = ok[which(params[ok,4] == uniques[i,4])]
  ok = ok[which(params[ok,5] == uniques[i,5])]
  ok = ok[which(params[ok,6] == uniques[i,6])]
  ok = ok[which(params[ok,7] == uniques[i,7])]
  ok = ok[which(params[ok,8] == uniques[i,8])]
  
#  cat (params[aux, ])
  if (length(ok) > 1) 
  {
    # Process this unique[i,] length(aux) times and replace in the correct locations of the output matrix
#    cat ("uniques[ i = ", i, "]:", params[aux[1],], "; duplicates:", aux, "\n")
    params[ok,]
  }
  else
  {
    # Treat as a scalar parameter set
    out <- prd_batch (rts, a=params[ok,1], z=params[ok,2], v=params[ok,3], t0=params[ok,4], d=params[ok,5], sz=params[ok,6], sv=params[ok,7], st0=params[ok,8])
  }
}
})



