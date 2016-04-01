context("Diffusion pdiffusion and ddiffusion functions.")

# Include the old non-vectorised pdiffusion and ddiffusion funcs to ensure they haven't been broken

# [MG 20150616]
# In line with LBA, adjust t0 to be the lower bound of the non-decision time distribution rather than the average 
# Called from pdiffusion, ddiffusion, rdiffusion 
recalc_t0 <- function (t0, st0) { t0 <- t0 + st0/2 }

old_ddiffusion <- function (t, boundary = c("upper", "lower"), 
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

old_pdiffusion <- function (t, boundary = c("upper", "lower"), 
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
  
  correct_pdiffusion_vals <- vector(length=test_vec_len)
  correct_ddiffusion_vals <- vector(length=test_vec_len)
  for (i in 1:test_vec_len)
  {
    correct_pdiffusion_vals[i] <- old_pdiffusion(vec_rts[i], boundary=vec_bounds[i], a =vec_a[i], z=vec_z[i], v =vec_v[i], 
                               t0=vec_t0[i], d=vec_d[i], sz =vec_sz[i], sv=vec_sv[i], st0 = vec_st0[i])
    correct_ddiffusion_vals[i] <- old_ddiffusion(vec_rts[i], boundary=vec_bounds[i], a =vec_a[i], z=vec_z[i], v =vec_v[i], 
                               t0=vec_t0[i], d=vec_d[i], sz =vec_sz[i], sv=vec_sv[i], st0 = vec_st0[i])
  }  
  
  pdiffusions <- pdiffusion (vec_rts, boundary=vec_bounds, a=vec_a, z=vec_z, v=vec_v, t0=vec_t0, d=vec_d, sz=vec_sz, sv=vec_sv, st0=vec_st0)
  ddiffusions <- ddiffusion (vec_rts, boundary=vec_bounds, a=vec_a, z=vec_z, v=vec_v, t0=vec_t0, d=vec_d, sz=vec_sz, sv=vec_sv, st0=vec_st0)
  # Note: allow a lot of tolerance for pdiffusion difference due to sampling error 
  #       (should never be as high as 1e-3, though)
  expect_that (pdiffusions, equals (correct_pdiffusion_vals, tolerance=1e-3))
  expect_that (ddiffusions, equals (correct_ddiffusion_vals))
  
  pdiffusions2 <- pdiffusion (vec_rts, boundary=vec_bounds, a=vec_a, z=vec_z[1:sample(test_vec_len, 1)], v=vec_v[1:sample(test_vec_len, 1)], t0=vec_t0[1:sample(test_vec_len, 1)], d=vec_d[1:sample(test_vec_len, 1)], sz=vec_sz[1:sample(test_vec_len, 1)], sv=vec_sv[1:sample(test_vec_len, 1)], st0=vec_st0[1:sample(test_vec_len, 1)])
  ddiffusions2 <- ddiffusion (vec_rts, boundary=vec_bounds, a=vec_a, z=vec_z[1:sample(test_vec_len, 1)], v=vec_v[1:sample(test_vec_len, 1)], t0=vec_t0[1:sample(test_vec_len, 1)], d=vec_d[1:sample(test_vec_len, 1)], sz=vec_sz[1:sample(test_vec_len, 1)], sv=vec_sv[1:sample(test_vec_len, 1)], st0=vec_st0[1:sample(test_vec_len, 1)])
  expect_that (pdiffusions2, equals (correct_pdiffusion_vals, tolerance=1e-3))
  expect_that (ddiffusions2, equals (correct_ddiffusion_vals))
})


test_that("ensure vectorised functions are equal to previous non-vectorised versions (v2):", {
  
  n_test <- 20
  rts <- rdiffusion(n_test, a=1, z=0.5, v=2, t0=0.5, d=0, sz = 0, sv = 0, st0 = 0)
  
  correct_pdiffusion_vals <- vector("numeric", n_test)
  correct_ddiffusion_vals <- vector("numeric", n_test)
  for (i in seq(1, n_test, by = 2))
  {
    correct_pdiffusion_vals[i] <- old_pdiffusion(sort(rts[, "rt"])[i], boundary=as.character(rts[i, "response"]), a=1, z=0.5, v=2, t0=0.5, d=0, sz = 0, sv = 0, st0 = 0)
    correct_pdiffusion_vals[i+1] <- old_pdiffusion(sort(rts[, "rt"])[i+1], boundary=as.character(rts[i+1, "response"]), a=1.5, z=0.75, v=2.25, t0=0.4, d=0.1, sz = 0.1, sv = 0.1, st0 = 0.1)
    correct_ddiffusion_vals[i] <- old_ddiffusion(rts[i, "rt"], boundary=as.character(rts[i, "response"]), a=1, z=0.5, v=2, t0=0.5, d=0, sz = 0, sv = 0, st0 = 0)
    correct_ddiffusion_vals[i+1] <- old_ddiffusion(rts[i+1, "rt"], boundary=as.character(rts[i+1, "response"]), a=1.5, z=0.75, v=2.25, t0=0.4, d=0.1, sz = 0.1, sv = 0.1, st0 = 0.1)
  }  
  
  pdiffusions <- pdiffusion(sort(rts[, "rt"]), boundary=as.character(rts[, "response"]), a=c(1, 1.5), z=c(0.5, 0.75), v=c(2, 2.25), t0=c(0.5, 0.4), d=c(0,0.1), sz = c(0,0.1), sv = c(0,0.1), st0 = c(0, 0.1))
  ddiffusions <- ddiffusion (rts[, "rt"], boundary=as.character(rts[, "response"]), a=c(1, 1.5), z=c(0.5, 0.75), v=c(2, 2.25), t0=c(0.5, 0.4), d=c(0,0.1), sz = c(0,0.1), sv = c(0,0.1), st0 = c(0, 0.1))
  # Note: allow a lot of tolerance for pdiffusion difference due to sampling error 
  #       (should never be as high as 1e-3, though)
  expect_equal(pdiffusions, correct_pdiffusion_vals, tolerance=1e-3)
  expect_equal(ddiffusions, correct_ddiffusion_vals)
  
})


test_that("diffusion functions work with numeric and factor boundaries", {
  n_test <- 20
  rts <- rdiffusion(n_test, a=1, z=0.5, v=2, t0=0.5, d=0, sz = 0, sv = 0, st0 = 0)
  expect_is(ddiffusion(rts$rt, boundary = rts$response, a=1, z=0.5, v=2, t0=0.5, d=0, sz = 0, sv = 0, st0 = 0), "numeric")
  expect_is(pdiffusion(sort(rts$rt), boundary = rts$response, a=1, z=0.5, v=2, t0=0.5, d=0, sz = 0, sv = 0, st0 = 0), "numeric")
  expect_is(ddiffusion(rts$rt, boundary = sample(1:2, 20, replace = TRUE), a=1, z=0.5, v=2, t0=0.5, d=0, sz = 0, sv = 0, st0 = 0), "numeric")
  expect_is(pdiffusion(sort(rts$rt), boundary = sample(1:2, 20, replace = TRUE), a=1, z=0.5, v=2, t0=0.5, d=0, sz = 0, sv = 0, st0 = 0), "numeric")
  expect_error(ddiffusion(rts$rt, sample(1:3, 20, replace = TRUE), a=1, z=0.5, v=2, t0=0.5, d=0, sz = 0, sv = 0, st0 = 0), "boundary")
  expect_error(pdiffusion(sort(rts$rt), sample(1:3, 20, replace = TRUE), a=1, z=0.5, v=2, t0=0.5, d=0, sz = 0, sv = 0, st0 = 0), "boundary")
})


test_that("qdiffusion is equivalent to manual calculation",{
  p11_fit <- structure(list(par = structure(c(1.32060063610882, 3.27271614698074, 0.338560144920614, 0.34996447540773, 0.201794924457386, 1.05516829794661), .Names = c("a", "v", "t0", "sz", "st0", "sv"))))
  q <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  prop_correct <- do.call(pdiffusion, args = c(rt = 20, as.list(p11_fit$par), boundary = "upper"))

  expect_equal(qdiffusion(q*prop_correct, boundary = "upper", a=p11_fit$par["a"], v=p11_fit$par["v"], t0=p11_fit$par["t0"], sz=p11_fit$par["sz"], st0=p11_fit$par["st0"], sv=p11_fit$par["sv"]),c(0.496360781904238, 0.573700985861292, 0.636165042003583, 0.714822449217467, 0.881706249900822))
  
  expect_equal(suppressWarnings(qdiffusion(q*prop_correct, boundary = "lower", a=p11_fit$par["a"], v=p11_fit$par["v"], t0=p11_fit$par["t0"], sz=p11_fit$par["sz"], st0=p11_fit$par["st0"], sv=p11_fit$par["sv"])),as.numeric(rep(NA, 5)))
})
