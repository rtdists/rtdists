#' Deprecated functions
#'
#' These functions have been renamed and deprecated in \pkg{rtdists}:
#' \code{drd()} (use \code{\link{ddiffusion}()}), 
#' \code{prd()} (use \code{\link{pdiffusion}()}),
#' \code{rrd()} (use \code{\link{rdiffusion}()}).
#' @rdname deprecated
#' @keywords internal
#' @aliases rtdists-deprecated
#' @param ... arguments passed from the old functions to the new functions
#' @export
drd <- function(...) {
  .Deprecated("ddiffusion", "rtdists", "drd was renamed to ddiffusion and is now deprecated.")
  ddiffusion(...)
}
#' @rdname deprecated
#' @export
prd <- function(...) {
  .Deprecated("pdiffusion", "rtdists", "prd was renamed to pdiffusion and is now deprecated.")
  pdiffusion(...)
}
#' @rdname deprecated
#' @export
rrd <- function(...) {
  .Deprecated("rdiffusion", "rtdists", "rrd was renamed to rdiffusion and is now deprecated.")
  rdiffusion(...)
}
