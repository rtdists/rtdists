
#' Single accumulator of the racing diffusion model (RDM)
#' 
#' Density, distribution, and random generation functions for a single Wald 
#' accumulator, that is, a diffusion process with within-trial variability in 
#' drift rate racing towards a single threshold with a uniformly distributed 
#' start point.
#' 
#' @param rt a vector of RTs.
#' @param n desired number of observations (scalar integer).
#' @param A start point interval or evidence in accumulator before beginning of decision process. Start point varies from trial to trial in the interval [0, \code{A}] (uniform distribution). Average amount of evidence before evidence accumulation across trials is \code{A}/2.
#' @param b response threshold. (\code{b} - \code{A}/2) is a measure of "response caution". 
#' @param t0 non-decision time or response time constant (in seconds). Lower bound for the duration of all non-decisional processes (encoding and response execution).
#' @param st0 variability of non-decision time, such that \code{t0} is uniformly distributed between \code{t0} and \code{t0} + \code{st0}. Only available in random number generation function \code{rwald}.
#' @param v drift rate. Mean rate of evidence accumulation within a trial. Needs to be positive, accumulators with a negative rate never finish.
#' @param s within-trial standard deviation of the evidence accumulation process (i.e., the diffusion constant), scales \code{A}, \code{b}, and \code{v}. Needs to be fixed to a constant in most applications. Default is 1.
#'   
#' @details These functions are mainly for internal purposes. We do not recommend to use them. Use the high-level functions described in \code{\link{RDM}} instead.
#' 
#' @return All functions starting with a \code{d} return the density (PDF), all functions starting with \code{p} return the distribution function (CDF), and all functions starting with \code{r} return random response times and responses (in a \code{matrix}).
#' 
#' @note Density (i.e., \code{dwald}) and distribution (i.e., \code{pwald}) functions are vectorized for all parameters (i.e., in case parameters are not of the same length as \code{rt}, parameters are recycled). Furthermore, the random derivative function also accepts a matrix of length \code{n} in which each column corresponds to a accumulator specific value (see \code{\link{rRDM}} for a more user-friendly way).
#'
#' @references 
#' Tillman, G., Van Zandt, T., & Logan, G. D. (2020). Sequential sampling models without random between-trial variability: The racing diffusion model of speeded decision making. \emph{Psychonomic Bulletin & Review}, 27(5), 911-936. doi:10.3758/s13423-020-01719-6
#' 
#' Logan, G. D., Van Zandt, T., Verbruggen, F., & Wagenmakers, E.-J. (2014). On the ability to inhibit thought and action: General and special theories of an act of control. \emph{Psychological Review}, 121(1), 66-95. doi:10.1037/a0035230
#' 
#' @importFrom stats dnorm pnorm qnorm rnorm runif
#' 
#' @name single-RDM
#' 
#' @example examples/examples.srdm.R
#'
NULL

# The first passage time equations below are a reparameterisation of the pigt,
# digt, and rwaldt functions, Copyright (C) 2013 Trisha Van Zandt, distributed
# with Logan, Van Zandt, Verbruggen, and Wagenmakers (2014), with comments and
# changes by Andrew Heathcote and Gabriel Tillman. Their criterion k and half
# width a correspond to b - A/2 and A/2 here.

rlevy <- function(n=1, m=0, c=1) {
  if (any(c < 0)) stop("c must be positive")
  c/qnorm(1-runif(n)/2)^2 + m
}

rwaldt <- function(n, k, l, tiny = 1e-6) {
  flag <- l > tiny
  x <- rep(NA, times = n)
  if (any(!flag)) x[!flag] <- rlevy(sum(!flag), 0, k[!flag]^2)
  mu <- k/l
  lambda <- k^2
  y <- rnorm(sum(flag))^2
  mu0 <- mu[flag]
  lambda0 <- lambda[flag]
  x0 <- mu0 + mu0^2*y/(2*lambda0) - 
    sqrt(4*mu0*lambda0*y + mu0^2*y^2)*mu0/(2*lambda0)
  z <- runif(length(x0))
  test <- mu0/(mu0+x0)
  x0[z > test] <- mu0[z > test]^2/x0[z > test]
  x[flag] <- x0
  x[x < 0] <- max(x)
  x
}

make_r_wald <- function(v, n, b, A, n_v, t0, st0=0) {
  if (is.null(dim(A))) starts <- matrix(runif(min=0,max=A,n=n*n_v),ncol=n_v,byrow=TRUE)
  else starts <- apply(A, c(1,2), function(x) runif(min=0, max = x, 1))
  if (is.null(dim(b))) bs <- t((b-t(starts)))
  else bs <- b-starts
  ttf <- matrix(Inf, nrow = n, ncol = n_v)
  ok <- as.vector(v) >= 0
  if (any(ok)) 
    ttf[ok] <- rwaldt(sum(ok), k = as.vector(bs)[ok], l = as.vector(v)[ok])
  rt <- apply(ttf+t0+runif(min=0,max=st0,n=n),1,min)
  resp <- apply(ttf+t0,1,which.min)
  bad <- !is.finite(rt)
  if (any(bad)) {
    warning(paste(sum(bad),"infinite RTs removed and less than", n, "rts returned"))
    resp <- resp[!bad]
    rt <- rt[!bad]
  }
  cbind(rt=rt,response=resp)
}

####### Wald:

#' @rdname single-RDM
#' @export dwald
dwald <- function(rt, A, b, t0, v, s = 1) {
  check_vector(rt, A, b, t0, v, s)
  nn <- length(rt)
  A <- rep(A, length.out = nn)
  b <- rep(b, length.out = nn)
  t0 <- rep(t0, length.out = nn)
  v <- rep(v, length.out = nn)
  s <- rep(s, length.out = nn)
  if (any(b < A)) stop(error_message_b_smaller_A)
  dwald_core(rt = rt, A = A, b = b, t0 = t0, v = v, s = s, nn = nn)
}

dwald_core <- function(rt, A, b, t0, v, s = 1, nn) {
  A <- rep(A, length.out = nn)
  b <- rep(b, length.out = nn)
  t0 <- rep(t0, length.out = nn)
  v <- rep(v, length.out = nn)
  s <- rep(s, length.out = nn)
  A <- A/s
  b <- b/s
  v <- v/s
  rt <- rem_t0(rt, t0)
  k <- b - A/2
  a <- A/2
  out <- numeric(nn)
  ok <- (rt > 0) & (v >= 0) & is.finite(A) & is.finite(b)
  ok[is.na(ok)] <- FALSE
  if (!any(ok)) return(out)
  A_small <- ok & (a < 1e-10)
  v_small <- ok & !A_small & (v < 1e-10)
  rest <- ok & !A_small & !v_small
  
  if (any(A_small)) {
    tt <- rt[A_small]
    kk <- k[A_small]
    vv <- v[A_small]
    vv[vv < 1e-10] <- 0
    lambda <- kk^2
    e <- numeric(length(tt))
    nz <- vv != 0
    if (any(nz)) {
      mu <- kk[nz]/vv[nz]
      e[nz] <- -(lambda[nz]/(2*tt[nz]))*(tt[nz]^2/mu^2 - 2*tt[nz]/mu + 1)
    }
    if (any(!nz)) e[!nz] <- -.5*lambda[!nz]/tt[!nz]
    out[A_small] <- exp(e + .5*log(lambda) - .5*log(2*tt^3*pi))
  }
  
  if (any(rest)) {
    tt <- rt[rest]
    kk <- k[rest]
    aa <- a[rest]
    vv <- v[rest]
    sqr_t <- sqrt(tt)
    term_1a <- -(aa - kk + tt*vv)^2/(2*tt)
    term_1b <- -(aa + kk - tt*vv)^2/(2*tt)
    term_1 <- (exp(term_1a) - exp(term_1b))/sqrt(2*pi*tt)
    term_2a <- log(.5) + log(vv)
    term_2b <- 2*pnorm((-kk + aa)/sqr_t + sqr_t*vv) - 1
    term_2c <- 2*pnorm((kk + aa)/sqr_t - sqr_t*vv) - 1
    term_3 <- term_1 + exp(term_2a)*(term_2b + term_2c)
    tmp <- numeric(length(tt))
    pos <- !is.na(term_3) & term_3 > 0
    tmp[pos] <- exp(log(term_3[pos]) - log(2) - log(aa[pos]))
    out[rest] <- tmp
  }
  
  if (any(v_small)) {
    tt <- rt[v_small]
    kk <- k[v_small]
    aa <- a[v_small]
    term_1 <- -.5*(log(2) + log(pi) + log(tt))
    term_2 <- (kk - aa)^2/(2*tt)
    term_3 <- (kk + aa)^2/(2*tt)
    term_4 <- exp(-term_2) - exp(-term_3)
    tmp <- numeric(length(tt))
    pos <- !is.na(term_4) & term_4 > 0
    tmp[pos] <- exp(term_1[pos] + log(term_4[pos]) - log(2) - log(aa[pos]))
    out[v_small] <- tmp
  }
  
  out[which(is.na(out) | out < 0)] <- 0
  out
}

#' @rdname single-RDM
#' @export pwald
pwald <- function(rt, A, b, t0, v, s = 1) {
  check_vector(rt, A, b, t0, v, s)
  nn <- length(rt)
  A <- rep(A, length.out = nn)
  b <- rep(b, length.out = nn)
  t0 <- rep(t0, length.out = nn)
  v <- rep(v, length.out = nn)
  s <- rep(s, length.out = nn)
  if (any(b < A)) stop(error_message_b_smaller_A)
  pwald_core(rt = rt, A = A, b = b, t0 = t0, v = v, s = s, nn = nn)
}

pwald_core <- function(rt, A, b, t0, v, s = 1, nn) {
  A <- rep(A, length.out = nn)
  b <- rep(b, length.out = nn)
  t0 <- rep(t0, length.out = nn)
  v <- rep(v, length.out = nn)
  s <- rep(s, length.out = nn)
  A <- A/s
  b <- b/s
  v <- v/s
  rt <- rem_t0(rt, t0)
  k <- b - A/2
  a <- A/2
  out <- numeric(nn)
  ok <- (rt > 0) & (v >= 0) & is.finite(A) & is.finite(b)
  ok[is.na(ok)] <- FALSE
  if (!any(ok)) return(out)
  A_small <- ok & (a < 1e-10)
  v_small <- ok & !A_small & (v < 1e-10)
  rest <- ok & !A_small & !v_small
  
  if (any(A_small)) {
    tt <- rt[A_small]
    kk <- k[A_small]
    vv <- v[A_small]
    vv[vv < 1e-10] <- 0
    lambda <- kk^2
    mu <- kk/vv
    e <- exp(log(2*lambda) - log(mu))
    add <- sqrt(lambda/tt)*(1 + tt/mu)
    sub <- sqrt(lambda/tt)*(1 - tt/mu)
    out[A_small] <- exp(e + log(pnorm(add, lower.tail = FALSE))) + 
      pnorm(sub, lower.tail = FALSE)
  }
  
  if (any(rest)) {
    tt <- rt[rest]
    kk <- k[rest]
    aa <- a[rest]
    vv <- v[rest]
    sqr_t <- sqrt(tt)
    term_1a <- .5*log(tt) - .5*log(2*pi)
    term_1b <- exp(-((kk - aa - tt*vv)^2/tt)/2)
    term_1c <- exp(-((kk + aa - tt*vv)^2/tt)/2)
    term_1 <- exp(term_1a)*(term_1b - term_1c)
    term_2a <- exp(2*vv*(kk - aa) + log(pnorm(-(kk - aa + tt*vv)/sqr_t)))
    term_2b <- exp(2*vv*(kk + aa) + log(pnorm(-(kk + aa + tt*vv)/sqr_t)))
    term_2 <- aa + (term_2b - term_2a)/(2*vv)
    term_4a <- 2*pnorm((kk + aa)/sqr_t - sqr_t*vv) - 1
    term_4b <- 2*pnorm((kk - aa)/sqr_t - sqr_t*vv) - 1
    term_4c <- .5*(tt*vv - aa - kk + .5/vv)
    term_4d <- .5*(kk - aa - tt*vv - .5/vv)
    term_4 <- term_4c*term_4a + term_4d*term_4b
    out[rest] <- (term_4 + term_2 + term_1)/(2*aa)
  }
  
  if (any(v_small)) {
    tt <- rt[v_small]
    kk <- k[v_small]
    aa <- a[v_small]
    sqr_t <- sqrt(tt)
    term_5a <- 2*pnorm((kk + aa)/sqr_t) - 1
    term_5b <- 2*pnorm(-(kk - aa)/sqr_t) - 1
    term_5 <- (-(kk + aa)*term_5a - (kk - aa)*term_5b)/(2*aa)
    term_6a <- -.5*(kk + aa)^2/tt - .5*log(2) - .5*log(pi) + .5*log(tt) - log(aa)
    term_6b <- -.5*(kk - aa)^2/tt - .5*log(2) - .5*log(pi) + .5*log(tt) - log(aa)
    out[v_small] <- term_5 + 1 + exp(term_6b) - exp(term_6a)
  }
  
  out[which(is.na(out) | out < 0)] <- 0
  pmin(out, 1)
}

#' @rdname single-RDM
#' @export rwald
rwald <- function(n, A, b, t0, v, s = 1, st0 = 0) {
  check_single_arg(n)
  if (any(b < A)) stop(error_message_b_smaller_A)
  n_v <- if (is.null(dim(v))) length(v) else ncol(v)
  if (is.null(dim(v))) drifts <- matrix(rep(v, each = n), ncol = n_v) else drifts <- v
  make_r_wald(v = drifts/s, n = n, b = b/s, A = A/s, n_v = n_v, t0 = t0, st0 = st0)
}
