# TODO: Which bits of this are worth moving into C++? 
#      Easy but probably not much difference: 
#          recalc_t0
#          scale z and sz by a

recalc_t0 <- function (t0, st0) { t0 <- t0 + st0/2 }


rcpp_fast_pdiffusion <- function (rt, response = "upper",
                 a, v, t0, z = 0.5, d = 0, sz = 0, sv = 0, st0 = 0,
                 precision = 3, maxt = 1e4)
{
  if(any(missing(a), missing(v), missing(t0))) stop("a, v, and/or t0 must be supplied")

  rt[rt>maxt] <- maxt
  if(!all(rt == sort(rt)))  stop("rt needs to be sorted")

  # Convert boundaries to numeric
  nn <- length(rt)
  # Build parameter matrix
  # Convert boundaries to numeric
  if (is.character(response)) {
    response <- match.arg(response, choices=c("upper", "lower"),several.ok = TRUE)
    numeric_bounds <- ifelse(response == "upper", 2L, 1L)
    }
  else {
    response <- as.numeric(response)
    if(any(!(response %in% 1:2))) stop("response needs to be either 'upper', 'lower', or as.numeric(response) %in% 1:2!")
    numeric_bounds <- as.integer(response)
  }
  numeric_bounds <- rep(numeric_bounds, length.out = nn)
  # all parameters brought to length of rt
  a <- rep(a, length.out = nn)
  v <- rep(v, length.out = nn)
  t0 <- rep(t0, length.out = nn)
  z <- rep(z, length.out = nn)
  d <- rep(d, length.out = nn)
  sz <- rep(sz, length.out = nn)
  sv <- rep(sv, length.out = nn)
  st0 <- rep(st0, length.out = nn)
  t0 <- recalc_t0 (t0, st0)

  # bind params to matrix
  params <- cbind (a, v, t0, d, sz, sv, st0, z, numeric_bounds)


  # Check for illegal parameter values
  if(ncol(params)<9) stop("Not enough parameters supplied: probable attempt to pass NULL values?")
  if(!is.numeric(params)) stop("Parameters need to be numeric.")
  if (any(is.na(params)) || !all(is.finite(params))) stop("Parameters need to be numeric and finite.")

  pvalues <- vector("numeric", length=length(rt))
  uniques <- unique(params)
  for (i in seq_len(nrow(uniques))) {
    ok_rows <- apply(params, 1, identical, y = uniques[i,])
    
    pvalues[ok_rows] <- p_fastdm (rt[ok_rows], uniques[i,1:8], precision, uniques[i,9])
  }
  pvalues
}


rcpp_pdiffusion <- function (rt, response = "upper",
                             a, v, t0, z = 0.5, d = 0, sz = 0, sv = 0, st0 = 0,
                             precision = 3, maxt = 1e4)
{
  if(any(missing(a), missing(v), missing(t0))) stop("a, v, and/or t0 must be supplied")
  
  rt[rt>maxt] <- maxt
  if(!all(rt == sort(rt)))  stop("rt needs to be sorted")
  
  # Convert boundaries to numeric
  nn <- length(rt)
  # Build parameter matrix
  # Convert boundaries to numeric
  if (is.character(response)) {
    response <- match.arg(response, choices=c("upper", "lower"),several.ok = TRUE)
    numeric_bounds <- ifelse(response == "upper", 2L, 1L)
  }
  else {
    response <- as.numeric(response)
    if(any(!(response %in% 1:2))) stop("response needs to be either 'upper', 'lower', or as.numeric(response) %in% 1:2!")
    numeric_bounds <- as.integer(response)
  }
  numeric_bounds <- rep(numeric_bounds, length.out = nn)
  # all parameters brought to length of rt
  a <- rep(a, length.out = nn)
  v <- rep(v, length.out = nn)
  t0 <- rep(t0, length.out = nn)
  z <- rep(z, length.out = nn)
  d <- rep(d, length.out = nn)
  sz <- rep(sz, length.out = nn)
  sv <- rep(sv, length.out = nn)
  st0 <- rep(st0, length.out = nn)
  t0 <- recalc_t0 (t0, st0)
  
  # bind params to matrix
  params <- cbind (a, v, t0, d, sz, sv, st0, z, numeric_bounds)
  
  
  # Check for illegal parameter values
  if(ncol(params)<9) stop("Not enough parameters supplied: probable attempt to pass NULL values?")
  if(!is.numeric(params)) stop("Parameters need to be numeric.")
  if (any(is.na(params)) || !all(is.finite(params))) stop("Parameters need to be numeric and finite.")
  
  pvalues <- vector("numeric", length=length(rt))
  uniques <- unique(params)
  for (i in seq_len(nrow(uniques))) {
    ok_rows <- apply(params, 1, identical, y = uniques[i,])
    
    pvalues[ok_rows] <- p_precise_fastdm (rt[ok_rows], uniques[i,1:8], precision, uniques[i,9])
  }
  pvalues
}

