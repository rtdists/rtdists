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
vec_st0 <- rep (.08, test_vec_len)


# !! DEBUG ONLY
# For scalar:

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

prds <- prd_batch (vec_rts, boundary=vec_bounds, a=vec_a, z=vec_z, v=vec_v, t0=vec_t0, d=vec_d, sz=vec_sz, sv=vec_sv, st0=vec_st0)
drds <- drd_batch (vec_rts, boundary=vec_bounds, a=vec_a, z=vec_z, v=vec_v, t0=vec_t0, d=vec_d, sz=vec_sz, sv=vec_sv, st0=vec_st0)

all.equal (correct_prd_vals, prds)
all.equal (correct_drd_vals, drds)

rrds <- 