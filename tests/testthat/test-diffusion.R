context("Diffusion prd and drd functions.")

# Include the old non-vectorised prd and drd funcs to ensure they haven't been broken

# [MG 20150616]
# In line with LBA, adjust t0 to be the lower bound of the non-decision time distribution rather than the average 
# Called from prd, drd, rrd 
recalc_t0 <- function (t0, st0) { t0 <- t0 + st0/2 }

old_drd <- function (t, boundary = c("upper", "lower"), 
                 a, v, t0, z = 0.5, d = 0, sz = 0, sv = 0, st0 = 0, 
                 precision = 3)
{
  t0 <- recalc_t0 (t0, st0) 
  
  # Check for illegal parameter values
  if(any(missing(a), missing(v), missing(t0))) stop("a, v, and/ot t0 must be supplied")
  pl <- c(a,v,t0,d,sz,sv,st0,z)
  if(length(pl)!=8) stop("Each parameter needs to be of length 1.")
  if(!is.numeric(pl)) stop("Parameters need to be numeric.")
  if (any(is.na(pl)) || !all(is.finite(pl))) stop("Parameters need to be numeric and finite.")
  
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
  
  abs(unlist(output[6]))
  
}

old_prd <- function (t, boundary = c("upper", "lower"), 
                 a, v, t0, z = 0.5, d = 0, sz = 0, sv = 0, st0 = 0, 
                 precision = 3, maxt = 1e4) 
{
  t0 <- recalc_t0 (t0, st0) 
  
  # Check for illegal parameter values
  if(any(missing(a), missing(v), missing(t0))) stop("a, v, and/ot t0 must be supplied")
  pl <- c(a,v,t0,d,sz,sv,st0,z)
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

test_that("ensure vectorised functions are equal to previous non-vectorised versions:", {
  
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
  vec_st0 <- rep (.08, test_vec_len)
  
  correct_prd_vals <- vector(length=test_vec_len)
  correct_drd_vals <- vector(length=test_vec_len)
  for (i in 1:test_vec_len)
  {
    correct_prd_vals[i] <- old_prd(vec_rts[i], boundary=vec_bounds[i], a =vec_a[i], z=vec_z[i], v =vec_v[i], 
                               t0=vec_t0[i], d=vec_d[i], sz =vec_sz[i], sv=vec_sv[i], st0 = vec_st0[i])
    correct_drd_vals[i] <- old_drd(vec_rts[i], boundary=vec_bounds[i], a =vec_a[i], z=vec_z[i], v =vec_v[i], 
                               t0=vec_t0[i], d=vec_d[i], sz =vec_sz[i], sv=vec_sv[i], st0 = vec_st0[i])
  }  
  
  prds <- prd (vec_rts, boundary=vec_bounds, a=vec_a, z=vec_z, v=vec_v, t0=vec_t0, d=vec_d, sz=vec_sz, sv=vec_sv, st0=vec_st0)
  drds <- drd (vec_rts, boundary=vec_bounds, a=vec_a, z=vec_z, v=vec_v, t0=vec_t0, d=vec_d, sz=vec_sz, sv=vec_sv, st0=vec_st0)
  # Note: allow a lot of tolerance for prd difference due to sampling error 
  #       (should never be as high as 1e-3, though)
  expect_that (prds, equals (correct_prd_vals, tolerance=1e-3))
  expect_that (drds, equals (correct_drd_vals))
  
  prds2 <- prd (vec_rts, boundary=vec_bounds, a=vec_a, z=vec_z[1:sample(test_vec_len, 1)], v=vec_v[1:sample(test_vec_len, 1)], t0=vec_t0[1:sample(test_vec_len, 1)], d=vec_d[1:sample(test_vec_len, 1)], sz=vec_sz[1:sample(test_vec_len, 1)], sv=vec_sv[1:sample(test_vec_len, 1)], st0=vec_st0[1:sample(test_vec_len, 1)])
  drds2 <- drd (vec_rts, boundary=vec_bounds, a=vec_a, z=vec_z[1:sample(test_vec_len, 1)], v=vec_v[1:sample(test_vec_len, 1)], t0=vec_t0[1:sample(test_vec_len, 1)], d=vec_d[1:sample(test_vec_len, 1)], sz=vec_sz[1:sample(test_vec_len, 1)], sv=vec_sv[1:sample(test_vec_len, 1)], st0=vec_st0[1:sample(test_vec_len, 1)])
  expect_that (prds2, equals (correct_prd_vals, tolerance=1e-3))
  expect_that (drds2, equals (correct_drd_vals))
})


test_that("ensure vectorised functions are equal to previous non-vectorised versions (v2):", {
  
  n_test <- 20
  rts <- rrd(n_test, a=1, z=0.5, v=2, t0=0.5, d=0, sz = 0, sv = 0, st0 = 0)
  
  correct_prd_vals <- vector("numeric", n_test)
  correct_drd_vals <- vector("numeric", n_test)
  for (i in seq(1, n_test, by = 2))
  {
    correct_prd_vals[i] <- old_prd(sort(rts[, "rt"])[i], boundary=as.character(rts[i, "response"]), a=1, z=0.5, v=2, t0=0.5, d=0, sz = 0, sv = 0, st0 = 0)
    correct_prd_vals[i+1] <- old_prd(sort(rts[, "rt"])[i+1], boundary=as.character(rts[i+1, "response"]), a=1.5, z=0.75, v=2.25, t0=0.4, d=0.1, sz = 0.1, sv = 0.1, st0 = 0.1)
    correct_drd_vals[i] <- old_drd(rts[i, "rt"], boundary=as.character(rts[i, "response"]), a=1, z=0.5, v=2, t0=0.5, d=0, sz = 0, sv = 0, st0 = 0)
    correct_drd_vals[i+1] <- old_drd(rts[i+1, "rt"], boundary=as.character(rts[i+1, "response"]), a=1.5, z=0.75, v=2.25, t0=0.4, d=0.1, sz = 0.1, sv = 0.1, st0 = 0.1)
  }  
  
  prds <- prd(sort(rts[, "rt"]), boundary=as.character(rts[, "response"]), a=c(1, 1.5), z=c(0.5, 0.75), v=c(2, 2.25), t0=c(0.5, 0.4), d=c(0,0.1), sz = c(0,0.1), sv = c(0,0.1), st0 = c(0, 0.1))
  drds <- drd (rts[, "rt"], boundary=as.character(rts[, "response"]), a=c(1, 1.5), z=c(0.5, 0.75), v=c(2, 2.25), t0=c(0.5, 0.4), d=c(0,0.1), sz = c(0,0.1), sv = c(0,0.1), st0 = c(0, 0.1))
  # Note: allow a lot of tolerance for prd difference due to sampling error 
  #       (should never be as high as 1e-3, though)
  expect_equal(prds, correct_prd_vals, tolerance=1e-3)
  expect_equal(drds, correct_drd_vals)
  
})

