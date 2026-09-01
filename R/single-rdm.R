
#' Wald and shifted Wald distribution (single accumulator of the RDM)
#'
#' Density, distribution, quantile, and random generation functions for a single
#' Wald accumulator, that is, a diffusion process with within-trial variability
#' in drift rate racing towards a single threshold with a uniformly distributed
#' start point. With the default \code{A = 0} (i.e., no start point
#' variability) these are the functions of the \strong{shifted Wald}
#' distribution in the evidence accumulation parameterization.
#'
#' @param rt a vector of RTs.
#' @param p vector of probabilities.
#' @param n desired number of observations (scalar integer).
#' @param A start point interval or evidence in accumulator before beginning of decision process. Start point varies from trial to trial in the interval [0, \code{A}] (uniform distribution). Average amount of evidence before evidence accumulation across trials is \code{A}/2. Default is 0, which gives the shifted Wald distribution.
#' @param b response threshold. (\code{b} - \code{A}/2) is a measure of "response caution".
#' @param t0 non-decision time or response time constant (in seconds). Lower bound for the duration of all non-decisional processes (encoding and response execution).
#' @param st0 variability of non-decision time, such that \code{t0} is uniformly distributed between \code{t0} and \code{t0} + \code{st0}. Only supported for \code{A = 0} in \code{dwald}, \code{pwald}, and \code{qwald} (see Details).
#' @param v drift rate. Mean rate of evidence accumulation within a trial. Needs to be positive, accumulators with a negative rate never finish.
#' @param s within-trial standard deviation of the evidence accumulation process (i.e., the diffusion constant), scales \code{A}, \code{b}, and \code{v}. Needs to be fixed to a constant in most applications. Default is 1.
#' @param log logical. If \code{TRUE}, the log density is returned. Default is \code{FALSE}.
#' @param lower.tail logical. If \code{TRUE} (default), probabilities are \eqn{P[T \le t]}, otherwise \eqn{P[T > t]}.
#' @param log.p logical. If \code{TRUE}, probabilities \code{p} are given as (and returned on) the log scale. Default is \code{FALSE}.
#' @param interval a vector containing the end-points of the interval to be searched for the desired quantiles in \code{qwald} if no analytical solution is available (i.e., if \code{A > 0} or \code{st0 > 0}). Default is \code{c(0, 10)}.
#' @param max_diff numeric. Maximum acceptable difference between desired and observed probability of the quantile function (\code{qwald}) in the numerical case.
#'
#' @details For \code{A = 0} (the default) the first passage time of the
#'   accumulator follows an inverse Gaussian (Wald) distribution with mean
#'   \code{mu = b/v} and shape \code{lambda = (b/s)^2}, shifted by \code{t0}.
#'   This is the shifted Wald distribution as it is commonly used in the
#'   response time literature (e.g., Anders, Alario, & Van Maanen, 2016, who
#'   write \eqn{\alpha} for \code{b}, \eqn{\gamma} for \code{v}, and
#'   \eqn{\theta} for \code{t0}). For \code{A > 0} the start point of the
#'   accumulator varies uniformly in [0, \code{A}] across trials.
#'
#'   The inverse Gaussian density, distribution, and quantile functions are
#'   evaluated on the log scale throughout, following the formulations of Giner
#'   and Smyth (2016) as implemented in the \pkg{statmod} package. This matters
#'   in the tails and whenever \code{2*b*v/s^2} is large, which happens for
#'   example with the common choice \code{s = 0.1}.
#'
#'   Variability in non-decision time, \code{st0}, is only supported for
#'   \code{A = 0}. In that case the density and distribution function are
#'   available in closed form, \eqn{f(t) = [F(t - t_0) - F(t - t_0 -
#'   s_{t0})]/s_{t0}}, where \eqn{F} is the Wald distribution function, and
#'   analogously for the distribution function with \eqn{F} replaced by its
#'   integral. For \code{A > 0} the corresponding integral has no closed form;
#'   use \code{\link{dRDM}}/\code{\link{n1PDF}}, which integrate over \code{st0}
#'   numerically. Note that with \code{st0 > 0} the far upper tail of
#'   \code{pwald} is accurate in absolute, but not in relative terms.
#'
#'   These functions are mainly for internal purposes when used as accumulators
#'   of a race. For a race of Wald accumulators use the high-level functions
#'   described in \code{\link{RDM}} instead.
#'
#' @return \code{dwald} returns the density (PDF), \code{pwald} the distribution
#'   function (CDF), and \code{qwald} the quantile function (i.e., predicted
#'   RTs), all of the same length as \code{rt} or \code{p}, respectively.
#'   \code{rwald} returns random response times: a numeric vector if a single
#'   accumulator is requested (i.e., \code{v} is of length 1) and otherwise a
#'   \code{matrix} with columns \code{rt} and \code{response} giving the winner
#'   of the race.
#'
#' @note Density (i.e., \code{dwald}), distribution (i.e., \code{pwald}), and
#'   quantile (i.e., \code{qwald}) functions are vectorized for all parameters
#'   (i.e., in case parameters are not of the same length as \code{rt}/\code{p},
#'   parameters are recycled). Furthermore, the random derivative function also
#'   accepts a matrix of length \code{n} in which each column corresponds to a
#'   accumulator specific value (see \code{\link{rRDM}} for a more user-friendly
#'   way).
#'
#' @references
#' Tillman, G., Van Zandt, T., & Logan, G. D. (2020). Sequential sampling models without random between-trial variability: The racing diffusion model of speeded decision making. \emph{Psychonomic Bulletin & Review}, 27(5), 911-936. doi:10.3758/s13423-020-01719-6
#'
#' Logan, G. D., Van Zandt, T., Verbruggen, F., & Wagenmakers, E.-J. (2014). On the ability to inhibit thought and action: General and special theories of an act of control. \emph{Psychological Review}, 121(1), 66-95. doi:10.1037/a0035230
#'
#' Anders, R., Alario, F.-X., & Van Maanen, L. (2016). The shifted Wald distribution for response time data analysis. \emph{Psychological Methods}, 21(3), 309-327. doi:10.1037/met0000066
#'
#' Giner, G., & Smyth, G. K. (2016). statmod: Probability calculations for the inverse Gaussian distribution. \emph{The R Journal}, 8(1), 339-351. doi:10.32614/RJ-2016-024
#'
#' Michael, J. R., Schucany, W. R., & Haas, R. W. (1976). Generating random variates using transformations with multiple roots. \emph{The American Statistician}, 30(2), 88-90. doi:10.1080/00031305.1976.10479147
#'
#' @importFrom stats dnorm pnorm qnorm rnorm runif optimize uniroot qchisq qgamma
#'
#' @name single-RDM
#' @aliases Wald ShiftedWald shifted-Wald
#'
#' @example examples/examples.srdm.R
#'
NULL

###############################################################################
## Inverse Gaussian (Wald) primitives                                        ##
###############################################################################

## All primitives are parameterized by mean `mu` and shape `lambda` and assume
## that x/q/p and the parameters have already been recycled to a common length
## and contain only valid values (x > 0, lambda > 0, mu > 0). `mu = Inf` is
## allowed and gives the Levy limit (i.e., a zero drift rate). Everything is
## computed on the log scale following Giner and Smyth (2016), which is what
## statmod::dinvgauss()/pinvgauss()/qinvgauss() do.

## log density.  Note x/mu evaluates to 0 for mu = Inf, so the Levy limit
## -lambda/(2x) falls out of the general expression without a special case.
dig_log <- function(x, mu, lambda) {
  0.5 * (log(lambda) - log(2 * pi) - 3 * log(x) - lambda * (x/mu - 1)^2 / x)
}

## log distribution function.
pig_log <- function(x, mu, lambda, lower.tail = TRUE) {
  ## x = Inf makes sqrt(lambda/x) * (x/mu) an Inf * 0 product below
  inf_x <- !is.finite(x)
  if (any(inf_x)) {
    logp <- rep(if (lower.tail) 0 else -Inf, length(x))
    fin <- !inf_x
    if (any(fin))
      logp[fin] <- pig_log(x[fin], mu[fin], lambda[fin],
                           lower.tail = lower.tail)
    return(logp)
  }
  sqr <- sqrt(lambda / x)
  z1 <- sqr * (x/mu - 1)
  z2 <- -sqr * (x/mu + 1)
  aa <- pnorm(z1, lower.tail = lower.tail, log.p = TRUE)
  bb <- 2 * lambda/mu + pnorm(z2, log.p = TRUE)
  if (lower.tail) bb <- exp(bb - aa) else bb <- -exp(bb - aa)
  logp <- aa + log1p(bb)
  if (!lower.tail) {
    ## far upper tail: the two terms above cancel, use the asymptotic expansion
    ## for log pbar(q; 1, phi) instead (Giner & Smyth, 2016, p. 342, unnumbered;
    ## derived in their appendix). The condition below is statmod's and is
    ## stricter than the paper's phi^(-1/2)*(q - 1) > 1e5.
    y <- x/mu
    w <- x * lambda / (2 * mu^2)
    i <- is.finite(mu) & (y > 1e6) & (w > 5e5)
    i[is.na(i)] <- FALSE
    if (any(i)) {
      logp[i] <- lambda[i]/mu[i] - 0.5 * log(pi) - log(2 * mu[i]/lambda[i]) -
        1.5 * log1p(w[i]) - w[i]
    }
  }
  logp
}

## Integrated distribution function H(x) = int_0^x F(u) du, needed for st0.
##   H(x) = (x - mu) * Phi(z1) + (x + mu) * exp(2*lambda/mu) * Phi(z2)
## For very large mu the two terms are O(mu) while H is O(x), so the general
## expression cancels catastrophically. Its mu -> Inf limit is exact and used
## instead; the two agree to ~1e-8 relative at the switch point.
iig <- function(x, mu, lambda) {
  out <- numeric(length(x))
  levy <- !(is.finite(mu) & (x/mu > 1e-8))
  if (any(!levy)) {
    j <- !levy
    sqr <- sqrt(lambda[j] / x[j])
    z1 <- sqr * (x[j]/mu[j] - 1)
    z2 <- -sqr * (x[j]/mu[j] + 1)
    out[j] <- (x[j] - mu[j]) * pnorm(z1) +
      (x[j] + mu[j]) * exp(2 * lambda[j]/mu[j] + pnorm(z2, log.p = TRUE))
  }
  if (any(levy)) {
    sqr <- sqrt(lambda[levy] / x[levy])
    out[levy] <- 2 * (x[levy] + lambda[levy]) * pnorm(-sqr) -
      2 * sqrt(lambda[levy] * x[levy]) * dnorm(sqr)
  }
  out
}

## Quantile function via Newton iteration on the log-CDF scale, ported from
## statmod::qinvgauss() (Giner & Smyth, 2016). Works on the standardized scale
## (mean 1, dispersion phi = mu/lambda) and rescales at the end.
qig <- function(p, mu, lambda, lower.tail = TRUE, log.p = FALSE,
                maxit = 200L, tol = 1e-14) {
  n <- length(p)
  if (log.p) logp <- p else {
    p[p < 0 | p > 1] <- NA
    logp <- log(p)
  }
  p <- exp(logp)
  q <- rep_len(NA_real_, n)
  bad <- is.na(p) | is.na(mu) | is.na(lambda) | mu <= 0 | lambda <= 0
  if (lower.tail) {
    at_zero <- logp == -Inf
    at_inf <- logp == 0
  } else {
    at_zero <- logp == 0
    at_inf <- logp == -Inf
  }
  q[!bad & at_zero] <- 0
  q[!bad & at_inf] <- Inf
  ok <- !(bad | at_zero | at_inf)

  ## zero drift: Levy, i.e. lambda / chi^2_1
  inf_mu <- ok & !is.finite(mu)
  if (any(inf_mu)) {
    q[inf_mu] <- lambda[inf_mu] /
      qchisq(logp[inf_mu], df = 1, lower.tail = !lower.tail, log.p = TRUE)
    ok[inf_mu] <- FALSE
  }

  phi <- mu / lambda  # standardized dispersion (coefficient of variation^2)
  ## negligible variability: the gamma approximation is accurate and avoids
  ## Newton steps that would immediately underflow.
  small_cv <- ok & (phi < 1e-8)
  if (any(small_cv)) {
    q[small_cv] <- qgamma(logp[small_cv], shape = 1/phi[small_cv],
                          scale = phi[small_cv] * mu[small_cv],
                          lower.tail = lower.tail, log.p = TRUE)
    ok[small_cv] <- FALSE
  }
  if (!any(ok)) return(q)

  ph <- phi[ok]
  lp <- logp[ok]
  pp <- p[ok]
  lam <- 1/ph  # shape of the standardized inverse Gaussian (mean 1)

  ## starting value: the mode of the standardized density
  kappa <- 1.5 * ph
  x <- sqrt(1 + kappa^2) - kappa
  big <- kappa > 1000
  if (any(big)) {
    k1 <- 1/2/kappa[big]
    x[big] <- k1 * (1 - k1^2)
  }
  if (lower.tail) {
    left <- lp < -11.51
    right <- lp > -1e-5
  } else {
    left <- lp > -1e-5
    right <- lp < -11.51
  }
  if (any(left)) {
    ## left tail asymptote F(x) ~ Phi(-1/sqrt(phi*x)). This only holds when the
    ## coefficient of variation is large; for a concentrated distribution it
    ## overshoots to the right of the mode, from where Newton diverges (this is
    ## a bug in statmod <= 1.5.0, which replaces the mode unconditionally).
    x[left] <- pmin(x[left],
                    1/ph[left]/qnorm(lp[left], lower.tail = lower.tail,
                                     log.p = TRUE)^2)
  }
  if (any(right)) {
    alpha <- 1/ph[right]
    x[right] <- pmax(x[right], qgamma(lp[right], shape = alpha, rate = alpha,
                                      lower.tail = lower.tail, log.p = TRUE))
  }

  step <- function(x, p, logp, lam) {
    one <- rep_len(1, length(x))
    logF <- pig_log(x, mu = one, lambda = lam, lower.tail = lower.tail)
    dlogp <- logp - logF
    dp <- dlogp
    ## p - F, computed from the difference of logs when that difference is tiny
    small <- !is.na(dlogp) & abs(dlogp) < 1e-5
    dp[small] <- exp(logp[small] + log1p(-dlogp[small]/2)) * dlogp[small]
    dp[!small] <- p[!small] - exp(logF[!small])
    dp / exp(dig_log(x, mu = one, lambda = lam))
  }
  ## take the step, backtracking if it would leave the support
  advance <- function(x, dx) {
    nx <- if (lower.tail) x + dx else x - dx
    for (k in seq_len(60)) {
      bad <- !(nx > 0)
      bad[is.na(bad)] <- TRUE
      if (!any(bad)) break
      dx[bad] <- dx[bad]/2
      nx <- if (lower.tail) x + dx else x - dx
    }
    nx[!(nx > 0) | is.na(nx)] <- x[!(nx > 0) | is.na(nx)]
    nx
  }

  dx <- step(x, pp, lp, lam)
  dx[is.na(dx)] <- 0
  sdx <- sign(dx)
  x <- advance(x, dx)
  i <- abs(dx) > tol
  iter <- 0L
  while (any(i)) {
    iter <- iter + 1L
    if (iter > maxit) {
      warning("qwald: maximum number of iterations exceeded", call. = FALSE)
      break
    }
    dx <- step(x[i], pp[i], lp[i], lam[i])
    ## never step backwards: the CDF is monotone, so a sign flip means we
    ## overshot into the numerically flat region.
    dx[is.na(dx) | dx * sdx[i] < 0] <- 0
    x[i] <- advance(x[i], dx)
    i[i] <- (abs(dx)/pmax(x[i], 1)) > tol
  }
  q[ok] <- x * mu[ok]
  q
}

## Random generation via the transformation of Michael, Schucany, and Haas
## (1976), with the large-dispersion guard from statmod::rinvgauss().
rig <- function(n, mu, lambda) {
  mu <- rep_len(mu, n)
  lambda <- rep_len(lambda, n)
  out <- rep_len(NA_real_, n)
  ok <- is.finite(mu) & mu > 0 & is.finite(lambda) & lambda > 0
  ok[is.na(ok)] <- FALSE
  ## zero drift (mu = Inf): Levy, i.e. lambda / chi^2_1
  levy <- !is.finite(mu) & is.finite(lambda) & lambda > 0
  levy[is.na(levy)] <- FALSE
  if (any(levy)) out[levy] <- lambda[levy]/rnorm(sum(levy))^2
  if (!any(ok)) return(out)
  m <- sum(ok)
  phi <- mu[ok]/lambda[ok]  # standardized dispersion
  y <- rnorm(m)^2
  yphi <- y * phi
  x1 <- numeric(m)
  big <- yphi > 5e5
  if (any(big)) x1[big] <- 1/yphi[big]
  if (any(!big))
    x1[!big] <- 1 + yphi[!big]/2 * (1 - sqrt(1 + 4/yphi[!big]))
  first <- runif(m) < 1/(1 + x1)
  x1[!first] <- 1/x1[!first]
  out[ok] <- mu[ok] * x1
  out
}

###############################################################################
## Wald accumulator                                                          ##
###############################################################################

# The A > 0 first passage time equations below are a reparameterisation of the
# pigt and digt functions, Copyright (C) 2013 Trisha Van Zandt, distributed
# with Logan, Van Zandt, Verbruggen, and Wagenmakers (2014), with comments and
# changes by Andrew Heathcote and Gabriel Tillman. Their criterion k and half
# width a correspond to b - A/2 and A/2 here.

error_message_st0_A <- paste(
  "st0 > 0 is only supported for A = 0 (i.e., the shifted Wald).",
  "\nUse dRDM()/pRDM()/n1PDF()/n1CDF(), which integrate over st0 numerically."
)

make_r_wald <- function(v, n, b, A, n_v, t0, st0=0) {
  if (is.null(dim(A))) starts <- matrix(runif(min=0,max=A,n=n*n_v),ncol=n_v,byrow=TRUE)
  else starts <- apply(A, c(1,2), function(x) runif(min=0, max = x, 1))
  if (is.null(dim(b))) bs <- t((b-t(starts)))
  else bs <- b-starts
  ttf <- matrix(Inf, nrow = n, ncol = n_v)
  ok <- as.vector(v) >= 0
  if (any(ok)) {
    k <- as.vector(bs)[ok]
    vv <- as.vector(v)[ok]
    ## mu = Inf for a zero drift rate gives the Levy limit inside rig()
    ttf[ok] <- rig(sum(ok), mu = k/vv, lambda = k^2)
  }
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
dwald <- function(rt, A = 0, b, t0, v, s = 1, st0 = 0, log = FALSE) {
  check_vector(rt, A, b, t0, v, s, st0)
  nn <- length(rt)
  A <- rep(A, length.out = nn)
  b <- rep(b, length.out = nn)
  t0 <- rep(t0, length.out = nn)
  v <- rep(v, length.out = nn)
  s <- rep(s, length.out = nn)
  st0 <- rep(st0, length.out = nn)
  if (any(b < A)) stop(error_message_b_smaller_A)
  if (any(st0 > 0 & A > 0)) stop(error_message_st0_A)
  dwald_core(rt = rt, A = A, b = b, t0 = t0, v = v, s = s, nn = nn,
             st0 = st0, log = log)
}

dwald_core <- function(rt, A, b, t0, v, s = 1, nn, st0 = 0, log = FALSE) {
  A <- rep(A, length.out = nn)
  b <- rep(b, length.out = nn)
  t0 <- rep(t0, length.out = nn)
  v <- rep(v, length.out = nn)
  s <- rep(s, length.out = nn)
  st0 <- rep(st0, length.out = nn)
  A <- A/s
  b <- b/s
  v <- v/s
  rt_raw <- rt
  rt <- rem_t0(rt, t0)
  k <- b - A/2
  a <- A/2
  ## the density is accumulated on the log scale, so that `log = TRUE` does not
  ## lose the far tail to underflow.
  out <- rep(-Inf, nn)
  ok <- (rt > 0) & (v >= 0) & is.finite(A) & is.finite(b)
  ok[is.na(ok)] <- FALSE
  if (!any(ok)) return(if (log) out else exp(out))
  A_small <- ok & (a < 1e-10)
  v_small <- ok & !A_small & (v < 1e-10)
  rest <- ok & !A_small & !v_small

  if (any(A_small)) {
    ## shifted Wald: inverse Gaussian with mean k/v and shape k^2
    vv <- v[A_small]
    vv[vv < 1e-10] <- 0
    kk <- k[A_small]
    mu <- kk/vv          # Inf for vv == 0 (Levy limit)
    lambda <- kk^2
    ss <- st0[A_small]
    if (any(ss > 0)) {
      out[A_small] <- dwald_st0(rt_raw[A_small], t0[A_small], ss, mu, lambda)
    } else {
      out[A_small] <- dig_log(rt[A_small], mu, lambda)
    }
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
    tmp <- rep(-Inf, length(tt))
    pos <- !is.na(term_3) & term_3 > 0
    tmp[pos] <- log(term_3[pos]) - log(2) - log(aa[pos])
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
    tmp <- rep(-Inf, length(tt))
    pos <- !is.na(term_4) & term_4 > 0
    tmp[pos] <- term_1[pos] + log(term_4[pos]) - log(2) - log(aa[pos])
    out[v_small] <- tmp
  }

  out[is.na(out)] <- -Inf
  if (log) out else exp(out)
}

## Shifted Wald density with uniform non-decision time variability:
##   f(t) = [F(t - t0) - F(t - t0 - st0)] / st0
## Evaluated with the tail of F that does not cancel.
dwald_st0 <- function(rt, t0, st0, mu, lambda) {
  x1 <- rem_t0(rt, t0)
  x2 <- rem_t0(rt, t0 + st0)
  out <- rep(-Inf, length(rt))
  ok <- x1 > 0
  if (!any(ok)) return(out)
  x1 <- x1[ok]; x2 <- x2[ok]; mu <- mu[ok]; lambda <- lambda[ok]
  lo <- pig_log(x1, mu, lambda)
  ## F(x2) is 0 for rt below t0 + st0, i.e. while the uniform is only partly
  ## covered; pig_log() is undefined at 0 and must not be called there.
  started <- x2 > 0
  ## below the median of F(x1) the lower tails do not cancel; above it, take
  ## the difference of the upper tails instead.
  upper <- lo > log(0.5)
  d <- numeric(length(x1))
  j <- which(!upper)
  if (length(j)) {
    d[j] <- exp(lo[j])
    jj <- j[started[j]]
    if (length(jj)) d[jj] <- d[jj] - exp(pig_log(x2[jj], mu[jj], lambda[jj]))
  }
  j <- which(upper)
  if (length(j)) {
    ## the upper tail of F(x2) is 1 whenever x2 == 0
    hi2 <- numeric(length(j))
    jj <- started[j]
    if (any(jj))
      hi2[jj] <- pig_log(x2[j][jj], mu[j][jj], lambda[j][jj],
                         lower.tail = FALSE)
    d[j] <- exp(hi2) - exp(pig_log(x1[j], mu[j], lambda[j], lower.tail = FALSE))
  }
  d[!is.na(d) & d < 0] <- 0
  out[ok] <- log(d) - log(st0[ok])
  out
}

#' @rdname single-RDM
#' @export pwald
pwald <- function(rt, A = 0, b, t0, v, s = 1, st0 = 0, lower.tail = TRUE,
                  log.p = FALSE) {
  check_vector(rt, A, b, t0, v, s, st0)
  nn <- length(rt)
  A <- rep(A, length.out = nn)
  b <- rep(b, length.out = nn)
  t0 <- rep(t0, length.out = nn)
  v <- rep(v, length.out = nn)
  s <- rep(s, length.out = nn)
  st0 <- rep(st0, length.out = nn)
  if (any(b < A)) stop(error_message_b_smaller_A)
  if (any(st0 > 0 & A > 0)) stop(error_message_st0_A)
  pwald_core(rt = rt, A = A, b = b, t0 = t0, v = v, s = s, nn = nn,
             st0 = st0, lower.tail = lower.tail, log.p = log.p)
}

pwald_core <- function(rt, A, b, t0, v, s = 1, nn, st0 = 0, lower.tail = TRUE,
                       log.p = FALSE) {
  A <- rep(A, length.out = nn)
  b <- rep(b, length.out = nn)
  t0 <- rep(t0, length.out = nn)
  v <- rep(v, length.out = nn)
  s <- rep(s, length.out = nn)
  st0 <- rep(st0, length.out = nn)
  A <- A/s
  b <- b/s
  v <- v/s
  rt_raw <- rt
  rt <- rem_t0(rt, t0)
  k <- b - A/2
  a <- A/2
  ## `out` holds lower-tail probabilities on the natural scale for the A > 0
  ## branches (unchanged from the original implementation); the shifted Wald
  ## branch writes log probabilities (already tail-adjusted) into `logout`.
  out <- numeric(nn)
  logout <- rep(-Inf, nn)
  ok <- (rt > 0) & (v >= 0) & is.finite(A) & is.finite(b)
  ok[is.na(ok)] <- FALSE
  A_small <- ok & (a < 1e-10)
  v_small <- ok & !A_small & (v < 1e-10)
  rest <- ok & !A_small & !v_small

  if (any(A_small)) {
    vv <- v[A_small]
    vv[vv < 1e-10] <- 0
    kk <- k[A_small]
    mu <- kk/vv
    lambda <- kk^2
    ss <- st0[A_small]
    if (any(ss > 0)) {
      logout[A_small] <- pwald_st0(rt_raw[A_small], t0[A_small], ss, mu, lambda,
                                   lower.tail = lower.tail)
    } else {
      logout[A_small] <- pig_log(rt[A_small], mu, lambda,
                                 lower.tail = lower.tail)
    }
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
  out <- pmin(out, 1)
  ## an accumulator that drifts towards its threshold finishes eventually; the
  ## branches above all evaluate to Inf * 0 at rt = Inf.
  inf_rt <- ok & !is.finite(rt)
  if (any(inf_rt)) out[inf_rt] <- 1
  if (!lower.tail) out[!A_small] <- 1 - out[!A_small]
  if (log.p) {
    out <- log(out)
    out[A_small] <- logout[A_small]
  } else {
    out[A_small] <- exp(logout[A_small])
  }
  out
}

## Shifted Wald distribution function with uniform non-decision time
## variability, using H(x) = int_0^x F(u) du:
##   F(t) = [H(t - t0) - H(t - t0 - st0)] / st0
pwald_st0 <- function(rt, t0, st0, mu, lambda, lower.tail = TRUE) {
  x1 <- rem_t0(rt, t0)
  x2 <- rem_t0(rt, t0 + st0)
  out <- rep(-Inf, length(rt))
  ok <- x1 > 0
  if (!lower.tail) out[!ok] <- 0
  ## H(Inf) - H(Inf) is not defined; the distribution is proper
  inf_rt <- !is.finite(rt)
  ok <- ok & !inf_rt
  if (any(inf_rt)) out[inf_rt] <- if (lower.tail) 0 else -Inf
  if (!any(ok)) return(out)
  p <- (iig(x1[ok], mu[ok], lambda[ok]) - iig(x2[ok], mu[ok], lambda[ok])) /
    st0[ok]
  p[is.na(p)] <- 0
  p <- pmin(pmax(p, 0), 1)
  out[ok] <- log(if (lower.tail) p else 1 - p)
  out
}

#' @rdname single-RDM
#' @export qwald
qwald <- function(p, A = 0, b, t0, v, s = 1, st0 = 0, lower.tail = TRUE,
                  log.p = FALSE, interval = c(0, 10), max_diff = 1e-4) {
  check_vector(p, A, b, t0, v, s, st0)
  nn <- length(p)
  A <- rep(A, length.out = nn)
  b <- rep(b, length.out = nn)
  t0 <- rep(t0, length.out = nn)
  v <- rep(v, length.out = nn)
  s <- rep(s, length.out = nn)
  st0 <- rep(st0, length.out = nn)
  if (any(b < A)) stop(error_message_b_smaller_A)
  if (any(st0 > 0 & A > 0)) stop(error_message_st0_A)

  out <- rep(NA_real_, nn)
  ## analytical inversion whenever the accumulator is a plain shifted Wald.
  ## The condition mirrors the A_small branch of pwald_core(), which tests
  ## (A/s)/2 against 1e-10.
  ana <- (A/(2*s) < 1e-10) & (st0 == 0)
  ana[is.na(ana)] <- FALSE
  if (any(ana)) {
    vv <- v[ana]/s[ana]
    vv[vv < 1e-10] <- 0
    kk <- b[ana]/s[ana]
    out[ana] <- t0[ana] + qig(p[ana], mu = kk/vv, lambda = kk^2,
                              lower.tail = lower.tail, log.p = log.p)
  }
  if (any(!ana)) {
    ## no closed form: numerically invert pwald_core, following qLBA()
    idx <- which(!ana)
    pp <- p[idx]
    if (log.p) pp <- exp(pp)
    if (!lower.tail) pp <- 1 - pp
    for (i in seq_along(idx)) {
      j <- idx[i]
      out[j] <- q_invert_wald(pp[i], A = A[j], b = b[j], t0 = t0[j], v = v[j],
                              s = s[j], st0 = st0[j], interval = interval,
                              max_diff = max_diff)
    }
  }
  out
}

inv_cdf_wald <- function(x, A, b, t0, v, s, st0, value, abs = TRUE) {
  cur <- pwald_core(rt = x, A = A, b = b, t0 = t0, v = v, s = s, nn = length(x),
                    st0 = st0)
  if (abs) abs(value - cur) else (value - cur)
}

q_invert_wald <- function(value, A, b, t0, v, s, st0, interval, max_diff) {
  if (is.na(value) || value < 0 || value > 1) return(NA_real_)
  if (value == 0) return(t0)
  if (value == 1) return(Inf)
  args <- list(A = A, b = b, t0 = t0, v = v, s = s, st0 = st0, value = value)
  lo <- max(interval[1], t0)
  tmp <- do.call(optimize, c(list(f = inv_cdf_wald,
                                  interval = c(lo, interval[2]),
                                  tol = .Machine$double.eps^0.5), args))
  if (tmp$objective > max_diff)
    tmp <- do.call(optimize, c(list(f = inv_cdf_wald,
                                    interval = c(lo, max(interval)/2),
                                    tol = .Machine$double.eps^0.5), args))
  if (tmp$objective > max_diff) {
    try({
      uni <- do.call(uniroot, c(list(f = inv_cdf_wald,
                                     interval = c(lo, interval[2]),
                                     tol = .Machine$double.eps^0.5,
                                     abs = FALSE), args))
      tmp$objective <- abs(uni$f.root)
      tmp$minimum <- uni$root
    }, silent = TRUE)
  }
  if (tmp$objective > max_diff) {
    warning("Cannot obtain RT that is less than ", max_diff,
            " away from desired p = ", value,
            ".\nIncrease/decrease interval.", call. = FALSE)
    return(NA_real_)
  }
  tmp$minimum
}

#' @rdname single-RDM
#' @export rwald
rwald <- function(n, A = 0, b, t0, v, s = 1, st0 = 0) {
  check_single_arg(n)
  if (any(b < A)) stop(error_message_b_smaller_A)
  n_v <- if (is.null(dim(v))) length(v) else ncol(v)
  if (is.null(dim(v))) drifts <- matrix(rep(v, each = n), ncol = n_v) else drifts <- v
  out <- make_r_wald(v = drifts/s, n = n, b = b/s, A = A/s, n_v = n_v, t0 = t0,
                     st0 = st0)
  ## a single accumulator is a distribution, not a race: return the RTs
  if (n_v == 1) unname(out[, "rt"]) else out
}
