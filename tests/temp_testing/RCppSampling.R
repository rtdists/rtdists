# TODO: Which bits of this are worth moving into C++? 
#      Easy but probably not much difference: 
#          recalc_t0
#          scale z and sz by a

recalc_t0 <- function (t0, st0) { t0 <- t0 + st0/2 }

## When given vectorised parameters, n is the number of replicates for each parameter set
#' @rdname Diffusion
#' @export
rcpp_rdiffusion <- function (n, 
                        a, v, t0, z = 0.5*a, d = 0, sz = 0, sv = 0, st0 = 0, s = 1,
                        precision = 3)
{
  if(any(missing(a), missing(v), missing(t0))) stop("a, v, and/or t0 must be supplied")
  
  s <- rep(s, length.out = n)
  a <- rep(a, length.out = n)
  v <- rep(v, length.out = n)
  t0 <- rep(t0, length.out = n)
  z <- rep(z, length.out = n)
  z <- z/a  # transform z from absolute to relative scale (which is currently required by the C code)
  d <- rep(d, length.out = n)
  sz <- rep(sz, length.out = n)
  sz <- sz/a # transform sz from absolute to relative scale (which is currently required by the C code)
  sv <- rep(sv, length.out = n)
  st0 <- rep(st0, length.out = n)
  t0 <- recalc_t0 (t0, st0) 
  
  # Build parameter matrix (and divide a, v, and sv, by s)
  params <- cbind (a/s, v/s, t0, d, sz, sv/s, st0, z)
  
  # Check for illegal parameter values
  if(ncol(params)<8) stop("Not enough parameters supplied: probable attempt to pass NULL values?")
  if(!is.numeric(params)) stop("Parameters need to be numeric.")
  if (any(is.na(params)) || !all(is.finite(params))) stop("Parameters need to be numeric and finite.")
  
  randRTs    <- vector("numeric",length=n)
  randBounds <- vector("numeric",length=n)
  
  uniques <- unique(params)
  for (i in seq_len(nrow(uniques))) {
    ok_rows <- apply(params, 1, identical, y = uniques[i,])
    
    # Calculate n for this row
    current_n <- sum(ok_rows)
    
    out <- r_fastdm (current_n, uniques[i,1:8], precision)
    
    randRTs[ok_rows]    <- out$rt       # unlist(output[4])
    randBounds[ok_rows] <- out$boundary # unlist(output[5])
  }
  response <- factor(randBounds, levels = 0:1, labels = c("lower", "upper"))
  data.frame(rt = randRTs, response)
}



