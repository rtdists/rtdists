
#' Censored shifted Wald distribution
#'
#' Density, distribution, and random generation functions for the censored
#' shifted Wald model of choice response times (Miller, Scherbaum,
#' Heck, Goschke, & Enge, 2018): a competing-risks race (Prentice et al., 1978)
#' between two shifted Wald accumulators that share the non-decision time
#' \code{t0}. The accumulator of the upper response drifts towards its
#' threshold at rate \code{v}, the accumulator of the lower response at rate
#' \code{-v}; whichever finishes first determines response and RT and censors
#' the other process.
#'
#' @param rt a vector of RTs. Or for convenience also a \code{data.frame} with
#'   columns \code{rt} and \code{response} such as returned by \code{rcswald}.
#' @param n desired number of observations (scalar integer).
#' @param response character vector with values in \code{c("upper", "lower")}
#'   (possibly abbreviated), \code{"upper"} being the default. Alternatively, a
#'   numeric vector with values 1=lower and 2=upper. For convenience,
#'   \code{response} is converted via \code{as.numeric} also allowing factors.
#'   Ignored if the first argument is a \code{data.frame}.
#' @param b single response bound, as in the racing diffusion model and the
#'   LBA: each accumulator races towards its own threshold, and with an
#'   unbiased start point (\code{w = 0.5}) both thresholds equal \code{b}.
#'   Supply either \code{b} or \code{a}, not both.
#' @param a boundary separation, as in the Wiener diffusion model (see
#'   \code{\link{ddiffusion}}): the distance between the assumed bounds for
#'   lower and upper responses of a single diffusion process. Identical to the
#'   \code{b} parameterization with \code{a = 2*b}. Supply either \code{b} or
#'   \code{a}, not both.
#' @param t0 non-decision time or response time constant (in seconds), shared
#'   by both accumulators.
#' @param v drift rate of the upper accumulator; the lower accumulator drifts
#'   at \code{-v}. The sign of \code{v} therefore determines which response
#'   dominates: for \code{v > 0} most responses terminate at the upper bound,
#'   for \code{v < 0} at the lower bound.
#' @param w relative start point (bias). Must lie strictly between 0 and 1;
#'   the default 0.5 is an unbiased start. As in the Wiener diffusion model,
#'   \code{w > 0.5} moves the start point towards the upper bound: the
#'   thresholds of the upper and lower accumulator are \code{a*(1 - w)} and
#'   \code{a*w} (equivalently, \code{2*b*(1 - w)} and \code{2*b*w}).
#' @param s within-trial standard deviation of the evidence accumulation
#'   process (i.e., the diffusion constant), scales \code{b}/\code{a} and
#'   \code{v}. Needs to be fixed to a constant in most applications. Default is
#'   1.
#' @param log logical. If \code{TRUE}, the log density is returned. Default is
#'   \code{FALSE}.
#' @param lower.tail logical. If \code{TRUE} (the default), \code{pcswald}
#'   returns the defective probability of the response occurring at or before
#'   \code{rt}; if \code{FALSE}, the defective probability of the response
#'   occurring after \code{rt}. Both sum to the response probability.
#'
#' @details Each accumulator is a shifted Wald (inverse Gaussian) first passage
#'   time process, that is, a single Wald accumulator without start point
#'   variability (see \code{\link{dwald}} with \code{A = 0}). The accumulator
#'   that drifts away from its threshold still finishes on some trials because
#'   of the within-trial noise, and since the other accumulator finishes with
#'   probability 1, every trial produces either an upper or a lower response
#'   (commission errors, not omissions; censoring at an external deadline is
#'   not part of this model).
#'
#'   The functions are built on \code{\link{dwald}} and \code{\link{pwald}}:
#'   an accumulator with threshold \eqn{k} that drifts away from it at rate
#'   \eqn{-v} finishes with total probability \eqn{\exp(-2kv/s^2)}, and its
#'   density and distribution function are this constant times those of the
#'   mirrored accumulator with rate \eqn{v}. The defective density of a
#'   response is the density of its accumulator times the survival function of
#'   the other accumulator (Eq. 5 of Miller et al., 2018). For an unbiased
#'   start point (\code{w = 0.5}) the defective distribution function has a
#'   closed form; for \code{w != 0.5}, \code{pcswald} integrates the density
#'   numerically.
#'
#'   The density and distribution functions are defective, that is, they are
#'   scaled by the probability of the corresponding response, exactly
#'   as for \code{\link{dRDM}} and \code{\link{ddiffusion}}. Consequently,
#'   \code{sum(dcswald(rt, response, ..., log = TRUE))} is the log-likelihood of
#'   a full data set of RTs and responses, and \code{pcswald(Inf, response,
#'   ...)} is the probability of that response. Fitting only the
#'   correct-response RTs with a plain shifted Wald instead biases parameter
#'   estimates, because the censoring by errors is informative (Miller et al.,
#'   2018).
#'
#'   When drift rate and the disfavored bound are high, errors become
#'   negligible: the finishing probability of the losing accumulator drops
#'   below machine precision, its survival factor equals 1 in double precision,
#'   and the model reduces to the plain shifted Wald (Miller et al.'s "simple"
#'   censored shifted Wald, their Eq. 4). \code{dcswald} and \code{pcswald} are
#'   then equal to \code{\link{dwald}} and \code{\link{pwald}} with
#'   \code{A = 0}, so no version needs to be chosen.
#'
#'   Start point variability (\code{A}) and non-decision time variability
#'   (\code{st0}) are not part of this model.
#'
#' @return \code{dcswald} returns the defective density (PDF) and
#'   \code{pcswald} the defective distribution function (CDF), both of the
#'   same length as \code{rt}. \code{rcswald} returns a \code{data.frame} with
#'   columns \code{rt} (numeric) and \code{response} (factor with levels
#'   \code{"lower"} and \code{"upper"}).
#'
#' @note The density and distribution functions are vectorized for all
#'   parameters (i.e., parameters are recycled to the length of \code{rt}). \code{rcswald} recycles parameters to \code{n}.
#'
#' @references
#' Miller, R., Scherbaum, S., Heck, D. W., Goschke, T., & Enge, S. (2018). On the relation between the (censored) shifted Wald and the Wiener distribution as measurement models for choice response times. \emph{Applied Psychological Measurement}, 42(2), 116-135. doi:10.1177/0146621617710465
#'
#' Prentice, R. L., Kalbfleisch, J. D., Peterson, A. V., Flournoy, N., Farewell, V. T., & Breslow, N. E. (1978). The analysis of failure times in the presence of competing risks. \emph{Biometrics}, 34(4), 541-554. doi:10.2307/2530374
#'
#' @importFrom stats integrate runif
#'
#' @name CensoredShiftedWald
#' @aliases cswald censored-shifted-Wald
#'
#' @seealso \code{\link{dwald}} for a single (shifted) Wald accumulator,
#'   \code{\link{dRDM}} for a race of Wald accumulators towards a common bound,
#'   and \code{\link{ddiffusion}} for the Wiener diffusion model this model
#'   approximates.
#'
#' @example examples/examples.cswald.R
#'
NULL

###############################################################################
## internal helpers                                                          ##
###############################################################################

error_message_b_a <- paste(
  "Supply exactly one of b (single response bound, as in the RDM/LBA)",
  "or a (boundary separation, as in the Wiener diffusion model; a = 2*b)."
)

## check_vector() derives argument names from its own call, so it is invoked
## directly in the exported functions; this only settles the b/a choice.
csw_check <- function(b, a) {
  if (is.null(b) == is.null(a)) stop(error_message_b_a, call. = FALSE)
}

## resolve the two parameterizations into per-accumulator thresholds and
## recycle everything to a common length nn; all evidence-scale parameters are
## divided by s. Returns list(ku, kl, vv, t0) with ku/kl the thresholds of the
## upper/lower accumulator and vv the drift rate of the upper accumulator.
csw_prep <- function(b, a, w, t0, v, s, nn) {
  bound <- if (is.null(b)) a else 2 * b
  bound <- rep(bound, length.out = nn)
  t0 <- rep(t0, length.out = nn)
  v <- rep(v, length.out = nn)
  w <- rep(w, length.out = nn)
  s <- rep(s, length.out = nn)
  if (any(!is.na(w) & (w <= 0 | w >= 1)))
    stop("w must be strictly between 0 and 1!", call. = FALSE)
  if (any(!is.na(bound) & bound <= 0))
    stop("b/a must be positive!", call. = FALSE)
  if (any(!is.na(s) & s <= 0))
    stop("s must be positive!", call. = FALSE)
  list(ku = bound * (1 - w) / s, kl = bound * w / s, vv = v / s, t0 = t0)
}

## response -> integer (1 = lower, 2 = upper), following diffusion.R
csw_response <- function(response, nn) {
  if (is.character(response) || is.factor(response)) {
    response <- match.arg(as.character(response),
                          choices = c("upper", "lower"), several.ok = TRUE)
    response <- ifelse(response == "upper", 2L, 1L)
  } else {
    response <- as.numeric(response)
    if (any(!(response %in% 1:2)))
      stop("response needs to be either 'upper', 'lower', or as.numeric(response) %in% 1:2 (1 = lower, 2 = upper)!")
    response <- as.integer(response)
  }
  rep(response, length.out = nn)
}

## winner/loser parameters for a vector of responses (2 = upper, 1 = lower):
## kw/vw are threshold and drift rate of the accumulator of the observed
## response, kl/vl those of the other accumulator (vl = -vw)
csw_race <- function(pars, resp) {
  up <- resp == 2L
  vw <- ifelse(up, pars$vv, -pars$vv)
  list(kw = ifelse(up, pars$ku, pars$kl), vw = vw,
       kl = ifelse(up, pars$kl, pars$ku), vl = -vw)
}

## An accumulator that drifts away from its threshold k (v < 0) finishes with
## total probability exp(2*k*v) < 1, and its first passage density and CDF are
## this constant times those of the mirrored accumulator with drift -v > 0:
## the exponents (k + |v|t)^2 and (k - |v|t)^2 of the two densities differ by
## 4*k*|v|*t, which the 2t in the denominator reduces to the constant 2*k*|v|
## (Miller et al., 2018, Eqs. 2-3). For v >= 0 the constant is 1.
csw_mass <- function(k, v) exp(2 * k * pmin(v, 0))

## first passage density and CDF of one accumulator with threshold k and drift
## rate v of either sign (both on the s = 1 scale), via the single Wald
## accumulator without start point variability
csw_d1 <- function(rt, k, v, t0, nn) {
  csw_mass(k, v) *
    dwald_core(rt = rt, A = 0, b = k, t0 = t0, v = abs(v), s = 1, nn = nn)
}
csw_p1 <- function(rt, k, v, t0, nn) {
  csw_mass(k, v) *
    pwald_core(rt = rt, A = 0, b = k, t0 = t0, v = abs(v), s = 1, nn = nn)
}

## defective density of the response whose accumulator has threshold kw and
## drift vw, racing against the accumulator with threshold kl and drift vl:
## density of the winner times survival of the loser (Eq. 5). Once the loser's
## finishing probability is below one ulp its survival factor is exactly 1 and
## this is dwald(rt, A = 0, ...) (the "simple" censored shifted Wald, Eq. 4).
csw_d_core <- function(rt, kw, vw, kl, vl, t0, nn) {
  csw_d1(rt, kw, vw, t0, nn) * (1 - csw_p1(rt, kl, vl, t0, nn))
}

## defective CDF (lower.tail = TRUE) or defective survival function of the
## response. With equal thresholds (w = 0.5) both accumulators share the CDF F
## of the favored one, up to the constants cw/cl of csw_mass(), and
## integrating the defective density gives the closed forms below. Otherwise
## the defective density is integrated numerically, as in pLBA/pRDM.
csw_p_core <- function(rt, kw, vw, kl, vl, t0, nn, lower.tail = TRUE) {
  out <- rep(NA_real_, nn)
  same <- kw == kl
  same[is.na(same)] <- FALSE
  if (any(same)) {
    m <- sum(same)
    F_ <- pwald_core(rt = rt[same], A = 0, b = kw[same], t0 = t0[same],
                     v = abs(vw[same]), s = 1, nn = m)
    cw <- csw_mass(kw[same], vw[same])
    cl <- csw_mass(kl[same], vl[same])
    out[same] <- if (lower.tail) cw * (F_ - cl * F_^2 / 2)
                 else cw * ((1 - F_) - cl * (1 - F_^2) / 2)
  }
  for (i in which(!same)) {
    x <- rem_t0(rt[i], t0[i])
    if (is.na(x)) next
    if ((lower.tail && x <= 0) || (!lower.tail && !is.finite(x))) {
      out[i] <- 0
      next
    }
    dens <- function(u) csw_d_core(u, kw[i], vw[i], kl[i], vl[i],
                                   t0 = 0, nn = length(u))
    lo <- if (lower.tail) 0 else x
    hi <- if (lower.tail) x else Inf
    val <- tryCatch(
      integrate(dens, lower = lo, upper = hi,
                rel.tol = 1e-10, abs.tol = 1e-300)$value,
      error = function(e) NA_real_)
    if (is.na(val))
      val <- tryCatch(integrate(dens, lower = lo, upper = hi)$value,
                      error = function(e) NA_real_)
    out[i] <- val
  }
  pmin(pmax(out, 0), 1)
}

###############################################################################
## exported functions                                                        ##
###############################################################################

#' @rdname CensoredShiftedWald
#' @export dcswald
dcswald <- function(rt, response = "upper", b = NULL, t0, v, w = 0.5, s = 1,
                    a = NULL, log = FALSE) {
  ## for convenience accept data.frame as first argument
  if (is.data.frame(rt)) {
    response <- rt$response
    rt <- rt$rt
  }
  csw_check(b, a)
  if (is.null(a)) check_vector(rt, b, t0, v, w, s)
  else check_vector(rt, a, t0, v, w, s)
  nn <- length(rt)
  pars <- csw_prep(b = b, a = a, w = w, t0 = t0, v = v, s = s, nn = nn)
  race <- csw_race(pars, csw_response(response, nn))
  out <- csw_d_core(rt, race$kw, race$vw, race$kl, race$vl, pars$t0, nn)
  if (log) log(out) else out
}

#' @rdname CensoredShiftedWald
#' @export pcswald
pcswald <- function(rt, response = "upper", b = NULL, t0, v, w = 0.5, s = 1,
                    a = NULL, lower.tail = TRUE) {
  if (is.data.frame(rt)) {
    response <- rt$response
    rt <- rt$rt
  }
  csw_check(b, a)
  if (is.null(a)) check_vector(rt, b, t0, v, w, s)
  else check_vector(rt, a, t0, v, w, s)
  nn <- length(rt)
  pars <- csw_prep(b = b, a = a, w = w, t0 = t0, v = v, s = s, nn = nn)
  race <- csw_race(pars, csw_response(response, nn))
  csw_p_core(rt, race$kw, race$vw, race$kl, race$vl, pars$t0, nn,
             lower.tail = lower.tail)
}

#' @rdname CensoredShiftedWald
#' @export rcswald
rcswald <- function(n, b = NULL, t0, v, w = 0.5, s = 1, a = NULL) {
  check_single_arg(n)
  csw_check(b, a)
  if (is.null(a)) check_vector(b, t0, v, w, s)
  else check_vector(a, t0, v, w, s)
  pars <- csw_prep(b = b, a = a, w = w, t0 = t0, v = v, s = s, nn = n)
  av <- abs(pars$vv)
  ## The accumulator drifting towards its threshold (the upper one for v >= 0,
  ## the lower one for v < 0) always finishes. The other one finishes with
  ## probability exp(-2*k*|v|) and, conditional on finishing, its first
  ## passage time is that of the mirrored accumulator (see csw_mass), so it is
  ## a coin flip followed by an ordinary Wald draw.
  up <- !is.na(pars$vv) & pars$vv >= 0
  k_fav <- ifelse(up, pars$ku, pars$kl)
  k_dis <- ifelse(up, pars$kl, pars$ku)
  t_fav <- rwaldt(n, k = k_fav, l = av)
  t_dis <- rep(Inf, n)
  hit <- runif(n) < exp(-2 * k_dis * av)
  hit[is.na(hit)] <- FALSE
  if (any(hit)) t_dis[hit] <- rwaldt(sum(hit), k = k_dis[hit], l = av[hit])
  tu <- ifelse(up, t_fav, t_dis)
  tl <- ifelse(up, t_dis, t_fav)
  ## the race is always proper: at least one accumulator drifts towards its
  ## threshold (or both are in the zero-drift Levy limit), so min(tu, tl) is
  ## finite and the response proportions follow the sign of v.
  rt <- pmin(tu, tl) + pars$t0
  response <- factor(ifelse(tu <= tl, 1L, 0L), levels = 0:1,
                     labels = c("lower", "upper"))
  data.frame(rt = rt, response = response)
}
