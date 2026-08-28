
#' The Racing Diffusion Model (RDM)
#' 
#' Density, distribution, quantile, and random generation functions for the 
#' racing diffusion model (RDM) with any number of response alternatives.
#' 
#' @param rt vector of RTs. Or for convenience also a \code{data.frame} with
#'   columns \code{rt} and \code{response} (see details).
#' @param response integer vector of winning accumulators/responses
#'   corresponding to the vector of RTs/p (i.e., used for specifying the
#'   response for a given RT/probability). Will be recycled if necessary. Cannot
#'   contain values larger than the number of accumulators. First
#'   response/accumulator must receive value 1, second 2, and so forth.
#' @param p vector of probabilities. Or for convenience also a \code{data.frame}
#'   with columns \code{p} and \code{response}.
#' @param n desired number of observations (scalar integer).
#' @param A start point interval or evidence in accumulator before beginning of
#'   decision process. Start point varies from trial to trial in the interval
#'   [0, \code{A}] (uniform distribution). Average amount of evidence before
#'   evidence accumulation across trials is \code{A}/2.
#' @param b response threshold. (\code{b} - \code{A}/2) is a measure of
#'   "response caution".
#' @param t0 non-decision time or response time constant (in seconds). Lower
#'   bound for the duration of all non-decisional processes (encoding and
#'   response execution).
#' @param v drift rate. Mean rate of evidence accumulation within a trial. Needs
#'   to be positive, accumulators with a negative rate never finish.
#' @param s within-trial standard deviation of the evidence accumulation process
#'   (i.e., the diffusion constant), scales \code{A}, \code{b}, and \code{v}.
#'   Needs to be fixed to a constant in most applications. Default is 1.
#' @param st0 variability of non-decision time, such that \code{t0} is uniformly
#'   distributed between \code{t0} and \code{t0} + \code{st0}. Only available in
#'   \code{rRDM}.
#' @param silent logical. Should the number of accumulators used be suppressed?
#'   Default is \code{FALSE} which prints the number of accumulators.
#' @param interval a vector containing the end-points of the interval to be
#'   searched for the desired quantiles in \code{qRDM}. Default is \code{c(0,
#'   10)}.
#' @param scale_p logical. Should entered probabilities automatically be scaled
#'   by maximally predicted probability? Default is \code{FALSE}.
#' @param scale_max numerical scalar. Value at which maximally predicted RT
#'   should be calculated if \code{scale_p} is \code{TRUE}.
#' 
#' @details The RDM is a race between independent diffusion processes, one per
#'   response alternative, that all start at a point drawn uniformly from [0,
#'   \code{A}] and drift towards their own threshold \code{b} at rate \code{v}
#'   with within-trial noise \code{s}. The first accumulator to reach its
#'   threshold determines both the response and the decision time, to which the
#'   non-decision time \code{t0} is added. Unlike the LBA, the RDM has no
#'   between-trial variability in drift rate, the noise that produces both fast
#'   and slow errors is within-trial.
#'   
#'   \code{A}, \code{b}, \code{t0}, and \code{v} can either be a single numeric
#'   vector (which will be recycled across the accumulators) or a \code{list} of
#'   numeric vectors with one element per accumulator. In the latter case, each
#'   element is recycled to the number of trials. This is the same convention as
#'   for \code{\link{dLBA}}.
#'   
#'   The RDM is implemented as a race of Wald accumulators and shares the
#'   machinery of the LBA. \code{dRDM}, \code{pRDM}, \code{qRDM}, and
#'   \code{rRDM} are wrappers around \code{\link{dLBA}}, \code{\link{pLBA}},
#'   \code{\link{qLBA}}, and \code{\link{rLBA}} called with
#'   \code{distribution = "wald"}.
#' 
#' @return \code{dRDM} returns the density (PDF), \code{pRDM} returns the
#'   distribution function (CDF), \code{qRDM} returns the quantile function, and
#'   \code{rRDM} returns random response times and responses (in a
#'   \code{data.frame}).
#'   
#'   The density, distribution, and quantile functions are defective, that is,
#'   they are scaled by the probability of the corresponding response.
#' 
#' @note These are the top-level functions intended for end-users. To obtain the
#'   density and distribution of a single accumulator, see
#'   \code{\link{single-RDM}}.
#' 
#' @references 
#' Tillman, G., Van Zandt, T., & Logan, G. D. (2020). Sequential sampling models without random between-trial variability: The racing diffusion model of speeded decision making. \emph{Psychonomic Bulletin & Review}, 27(5), 911-936. \doi{10.3758/s13423-020-01719-6}
#' 
#' @seealso \code{\link{LBA}} for the linear ballistic accumulator and
#'   \code{\link{Diffusion}} for the Ratcliff diffusion model.
#' 
#' @name RDM
#' 
#' @example examples/examples.rdm.R
#'
NULL

#' @rdname RDM
#' @export 
dRDM <- function(rt, response, A, b, t0, v, s = 1, st0 = 0, silent = FALSE) {
  dLBA(rt = rt, response = response, A = A, b = b, t0 = t0, v = v, st0 = st0, 
       distribution = "wald", args.dist = list(s = s), silent = silent)
}

#' @rdname RDM
#' @export 
pRDM <- function(rt, response, A, b, t0, v, s = 1, st0 = 0, silent = FALSE) {
  pLBA(rt = rt, response = response, A = A, b = b, t0 = t0, v = v, st0 = st0, 
       distribution = "wald", args.dist = list(s = s), silent = silent)
}

#' @rdname RDM
#' @export 
qRDM <- function(p, response, A, b, t0, v, s = 1, st0 = 0, silent = FALSE, 
                 interval = c(0, 10), scale_p = FALSE, scale_max = Inf) {
  qLBA(p = p, response = response, A = A, b = b, t0 = t0, v = v, st0 = st0, 
       distribution = "wald", args.dist = list(s = s), silent = silent, 
       interval = interval, scale_p = scale_p, scale_max = scale_max)
}

#' @rdname RDM
#' @export
rRDM <- function(n, A, b, t0, v, s = 1, st0 = 0, silent = FALSE) {
  rLBA(n = n, A = A, b = b, t0 = t0, v = v, st0 = st0, 
       distribution = "wald", args.dist = list(s = s), silent = silent)
}
