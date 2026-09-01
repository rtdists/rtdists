
#' Censored shifted Wald distribution
#'
#' Density, distribution, quantile, and random generation functions for the
#' censored shifted Wald model of choice response times (Miller, Scherbaum,
#' Heck, Goschke, & Enge, 2018): a competing-risks race (Prentice et al., 1978)
#' between two shifted Wald accumulators that share the non-decision time
#' \code{t0}. The accumulator of the upper response drifts towards its
#' threshold at rate \code{v}, the accumulator of the lower response at rate
#' \code{-v}; whichever finishes first determines response and RT and censors
#' the other process.
#'
#' @param rt a vector of RTs. Or for convenience also a \code{data.frame} with
#'   columns \code{rt} and \code{response} such as returned by \code{rcswald}.
#' @param p vector of probabilities. Or for convenience also a
#'   \code{data.frame} with columns \code{p} and \code{response}.
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
#' @param interval a vector containing the end-points of the interval to be
#'   searched for the desired quantiles in \code{qcswald} when no analytical
#'   solution is available. Default is \code{c(0, 10)}.
#' @param max_diff numeric. Maximum acceptable difference between desired and
#'   observed probability of the quantile function (\code{qcswald}) in the
#'   numerical case.
#' @param scale_p logical. Should entered probabilities automatically be scaled
#'   by maximally predicted probability? Default is \code{FALSE}.
#' @param scale_max numerical scalar. Value at which maximally predicted RT
#'   should be calculated if \code{scale_p} is \code{TRUE}.
#'
#' @details Each accumulator is a shifted Wald (inverse Gaussian) first passage
#'   time process (see \code{\link{dwald}}). The accumulator that drifts away
#'   from its threshold still finishes on some trials because of the
#'   within-trial noise, and since the other accumulator finishes with
#'   probability 1, every trial produces either an upper or a lower response
#'   (commission errors, not omissions; censoring at an external deadline is
#'   not part of this model).
#'
#'   The density, distribution, and quantile functions are defective, that is,
#'   they are scaled by the probability of the corresponding response (Eq. 5 of
#'   Miller et al., 2018), exactly as for \code{\link{dRDM}} and
#'   \code{\link{ddiffusion}}. Consequently, \code{sum(dcswald(rt, response,
#'   ..., log = TRUE))} is the log-likelihood of a full data set of RTs and
#'   responses, and \code{pcswald(Inf, response, ...)} is the probability of
#'   that response. Fitting only the correct-response RTs with a plain shifted
#'   Wald instead biases parameter estimates, because the censoring by errors
#'   is informative (Miller et al., 2018).
#'
#'   When drift rate and the disfavored bound are high, errors become
#'   negligible and the model reduces to the plain shifted Wald (Miller et
#'   al.'s "simple" censored shifted Wald, their Eq. 4). The functions detect
#'   this internally and then coincide with
#'   \code{\link{dwald}}/\code{\link{pwald}}/\code{\link{qwald}}, so no
#'   version needs to be chosen.
#'
#'   Start point variability (\code{A}) and non-decision time variability
#'   (\code{st0}) are not part of this model.
#'
#' @return \code{dcswald} returns the defective density (PDF), \code{pcswald}
#'   the defective distribution function (CDF), and \code{qcswald} the
#'   quantile function (i.e., predicted RTs for the defective probabilities
#'   \code{p}), all of the same length as \code{rt}/\code{p}. \code{qcswald}
#'   returns \code{NA} (with a warning) for probabilities above the predicted
#'   response probability. \code{rcswald} returns a \code{data.frame} with
#'   columns \code{rt} (numeric) and \code{response} (factor with levels
#'   \code{"lower"} and \code{"upper"}).
#'
#' @note The density, distribution, and quantile functions are vectorized for
#'   all parameters (i.e., parameters are recycled to the length of
#'   \code{rt}/\code{p}). \code{rcswald} recycles parameters to \code{n}.
#'
#' @references
#' Miller, R., Scherbaum, S., Heck, D. W., Goschke, T., & Enge, S. (2018). On the relation between the (censored) shifted Wald and the Wiener distribution as measurement models for choice response times. \emph{Applied Psychological Measurement}, 42(2), 116-135. doi:10.1177/0146621617710465
#'
#' Prentice, R. L., Kalbfleisch, J. D., Peterson, A. V., Flournoy, N., Farewell, V. T., & Breslow, N. E. (1978). The analysis of failure times in the presence of competing risks. \emph{Biometrics}, 34(4), 541-554. doi:10.2307/2530374
#'
#' @importFrom stats integrate
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

## An accumulator that drifts away from its threshold (drift v < 0) finishes
## with total probability exp(2*k*v/s^2) < 1. Once that mass is below one ulp,
## its survival function is 1 to machine precision for every t, and the
## competing-risks expressions for the other response reduce to the plain
## shifted Wald ("simple" censored shifted Wald, Miller et al., 2018, Eq. 4).
## The relative error of dropping the survival factor is bounded by the total
## mass itself (verified numerically, see issues/cswald-simple-vs-crisk.md),
## so below this cutoff the switch is exact in double precision.
csw_log_mass_cutoff <- log(.Machine$double.eps / 2)  # about -36.7

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

## mean of the inverse Gaussian first passage time; negative for a negative
## drift rate (the defective case, handled by analytic continuation in
## dig_log/pig_log), Inf for a zero drift rate (the Levy limit).
csw_acc_mu <- function(k, v) {
  v[abs(v) < 1e-10] <- 0
  k/v
}

## log first passage density / log survival of one accumulator, any sign of v
csw_dlog <- function(x, k, v) dig_log(x, csw_acc_mu(k, v), k^2)
csw_slog <- function(x, k, v) pig_log(x, csw_acc_mu(k, v), k^2,
                                      lower.tail = FALSE)

## log defective density of the response whose accumulator has threshold kw
## and drift vw, racing against the accumulator with threshold kl, drift vl
## (= -vw). x is the decision time (rt - t0).
csw_dlog_core <- function(x, kw, vw, kl, vl) {
  out <- rep(-Inf, length(x))
  ok <- (x > 0) & is.finite(x)
  ok[is.na(ok)] <- FALSE
  if (!any(ok)) return(out)
  lf <- csw_dlog(x[ok], kw[ok], vw[ok])
  ls <- numeric(sum(ok))
  ## survival factor of the losing accumulator, dropped when its total
  ## finishing probability is below one ulp (see csw_log_mass_cutoff)
  need <- !(2 * kl[ok] * vl[ok] < csw_log_mass_cutoff)
  need[is.na(need)] <- TRUE
  if (any(need)) ls[need] <- csw_slog(x[ok][need], kl[ok][need], vl[ok][need])
  out[ok] <- lf + ls
  out
}

## defective CDF at decision time x (natural scale). Closed form (shifted Wald
## CDF) when the losing accumulator's mass is negligible, otherwise numerical
## integration of the defective density.
csw_p_core <- function(x, kw, vw, kl, vl) {
  out <- numeric(length(x))
  ok <- (x > 0)
  ok[is.na(ok)] <- FALSE
  if (!any(ok)) return(out)
  simple <- ok & (2 * kl * vl < csw_log_mass_cutoff)
  simple[is.na(simple)] <- FALSE
  if (any(simple))
    out[simple] <- exp(pig_log(x[simple], csw_acc_mu(kw[simple], vw[simple]),
                               kw[simple]^2))
  for (i in which(ok & !simple)) {
    intgr <- function(u) {
      m <- length(u)
      exp(csw_dlog_core(u, rep_len(kw[i], m), rep_len(vw[i], m),
                        rep_len(kl[i], m), rep_len(vl[i], m)))
    }
    val <- tryCatch(
      integrate(intgr, lower = 0, upper = x[i],
                rel.tol = 1e-10, abs.tol = 1e-300)$value,
      error = function(e) NA_real_)
    if (is.na(val))
      val <- tryCatch(integrate(intgr, lower = 0, upper = x[i])$value,
                      error = function(e) NA_real_)
    out[i] <- val
  }
  pmin(pmax(out, 0), 1)
}

## winner/loser parameters for a vector of responses (2 = upper, 1 = lower)
csw_race <- function(pars, resp) {
  up <- resp == 2L
  vw <- ifelse(up, pars$vv, -pars$vv)
  list(kw = ifelse(up, pars$ku, pars$kl), vw = vw,
       kl = ifelse(up, pars$kl, pars$ku), vl = -vw)
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
  out <- csw_dlog_core(rem_t0(rt, pars$t0),
                       race$kw, race$vw, race$kl, race$vl)
  if (log) out else exp(out)
}

#' @rdname CensoredShiftedWald
#' @export pcswald
pcswald <- function(rt, response = "upper", b = NULL, t0, v, w = 0.5, s = 1,
                    a = NULL) {
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
  ## rem_t0 turns rt = Inf into Inf, giving the total response probability
  csw_p_core(rem_t0(rt, pars$t0), race$kw, race$vw, race$kl, race$vl)
}

#' @rdname CensoredShiftedWald
#' @export qcswald
qcswald <- function(p, response = "upper", b = NULL, t0, v, w = 0.5, s = 1,
                    a = NULL, interval = c(0, 10), max_diff = 1e-4,
                    scale_p = FALSE, scale_max = Inf) {
  if (is.data.frame(p)) {
    response <- p$response
    p <- p$p
  }
  csw_check(b, a)
  if (is.null(a)) check_vector(p, b, t0, v, w, s)
  else check_vector(p, a, t0, v, w, s)
  nn <- length(p)
  pars <- csw_prep(b = b, a = a, w = w, t0 = t0, v = v, s = s, nn = nn)
  race <- csw_race(pars, csw_response(response, nn))
  if (scale_p)
    p <- p * csw_p_core(rem_t0(rep(scale_max, nn), pars$t0),
                        race$kw, race$vw, race$kl, race$vl)
  out <- rep(NA_real_, nn)
  bad <- is.na(p) | p < 0 | p > 1
  ## closed form: losing accumulator negligible, so the defective CDF is the
  ## plain shifted Wald CDF of the winning accumulator (response probability
  ## is 1 up to machine precision)
  ana <- !bad & (2 * race$kl * race$vl < csw_log_mass_cutoff)
  ana[is.na(ana)] <- FALSE
  if (any(ana))
    out[ana] <- pars$t0[ana] +
      qig(p[ana], mu = csw_acc_mu(race$kw[ana], race$vw[ana]),
          lambda = race$kw[ana]^2)
  for (i in which(!bad & !ana)) {
    out[i] <- q_invert_cswald(p[i], kw = race$kw[i], vw = race$vw[i],
                              kl = race$kl[i], vl = race$vl[i],
                              t0 = pars$t0[i], interval = interval,
                              max_diff = max_diff)
  }
  out
}

q_invert_cswald <- function(value, kw, vw, kl, vl, t0, interval, max_diff) {
  if (value == 0) return(t0)
  pfun <- function(rt) {
    m <- length(rt)
    csw_p_core(rem_t0(rt, rep_len(t0, m)), rep_len(kw, m), rep_len(vw, m),
               rep_len(kl, m), rep_len(vl, m))
  }
  mass <- pfun(Inf)
  if (is.na(mass) || value > mass) {
    warning("p = ", value, " is above the predicted response probability (",
            format(mass), "); NA returned. Use pcswald(Inf, ...) to obtain ",
            "response probabilities or scale_p = TRUE.", call. = FALSE)
    return(NA_real_)
  }
  lo <- max(interval[1], t0)
  root <- tryCatch(
    uniroot(function(x) pfun(x) - value, interval = c(lo, interval[2]),
            extendInt = "upX", tol = .Machine$double.eps^0.5),
    error = function(e) NULL)
  if (is.null(root) || abs(root$f.root) > max_diff) {
    warning("Cannot obtain RT that is less than ", max_diff,
            " away from desired p = ", value,
            ".\nIncrease/decrease interval.", call. = FALSE)
    return(NA_real_)
  }
  root$root
}

#' @rdname CensoredShiftedWald
#' @export rcswald
rcswald <- function(n, b = NULL, t0, v, w = 0.5, s = 1, a = NULL) {
  check_single_arg(n)
  csw_check(b, a)
  if (is.null(a)) check_vector(b, t0, v, w, s)
  else check_vector(a, t0, v, w, s)
  pars <- csw_prep(b = b, a = a, w = w, t0 = t0, v = v, s = s, nn = n)
  tu <- csw_r_acc(n, pars$ku, pars$vv)
  tl <- csw_r_acc(n, pars$kl, -pars$vv)
  ## the race is always proper: at least one accumulator drifts towards its
  ## threshold (or both are in the zero-drift Levy limit), so min(tu, tl) is
  ## finite and the response proportions follow the sign of v.
  rt <- pmin(tu, tl) + pars$t0
  response <- factor(ifelse(tu <= tl, 1L, 0L), levels = 0:1,
                     labels = c("lower", "upper"))
  data.frame(rt = rt, response = response)
}

## first passage time of a single accumulator, any sign of the drift rate.
## For v < 0 the distribution is defective: the accumulator finishes at all
## only with probability exp(2*k*v), and conditional on finishing its first
## passage time is the ordinary inverse Gaussian of the mirrored (positive)
## drift, because F_defective(t | k, v) = exp(2*k*v) * F_IG(t | k/|v|, k^2).
## Sampling by inversion or rejection would be numerically unstable or
## arbitrarily slow; this two-stage scheme costs one runif per defective draw
## and one rig per hit.
csw_r_acc <- function(n, k, v) {
  v[abs(v) < 1e-10] <- 0
  out <- rep(Inf, n)
  hit <- rep(TRUE, n)
  neg <- !is.na(v) & v < 0
  if (any(neg)) hit[neg] <- runif(sum(neg)) < exp(2 * k[neg] * v[neg])
  hit[is.na(k) | is.na(v)] <- NA
  ok <- !is.na(hit) & hit
  ## mu = k/|v| is Inf for a zero drift rate, which rig() handles (Levy)
  if (any(ok)) out[ok] <- rig(sum(ok), mu = k[ok]/abs(v[ok]), lambda = k[ok]^2)
  out[is.na(hit)] <- NA_real_
  out
}
